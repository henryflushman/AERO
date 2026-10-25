"""Orbit propagation helpers for two-body, J2, and custom dynamics.

``propagate_two_body`` supports two propagation styles:

1. SciPy ``solve_ivp`` methods such as ``DOP853`` (the default), ``RK45``,
   ``Radau``, etc. These use ``rtol`` and ``atol``.
2. ``method="universal_anomaly"`` for exact two-body Kepler propagation with
   universal variables. This uses ``time_step_s`` for the output spacing
   instead of ``rtol``/``atol``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Iterable

import numpy as np
from scipy.integrate import solve_ivp

try:
    from .Constants import J2_EARTH, MU_EARTH_KM3_S2, R_EARTH_KM
    from .OrbitalElements import StateVector
except ImportError:
    from Constants import J2_EARTH, MU_EARTH_KM3_S2, R_EARTH_KM
    from OrbitalElements import StateVector


DEFAULT_RTOL = 1e-10
DEFAULT_ATOL = 1e-12
UNIVERSAL_SOLVER_TOLERANCE = 1e-12
UNIVERSAL_MAX_ITERATIONS = 100


@dataclass(frozen=True)
class PropagationResult:
    """Common propagation result for numerical and universal propagation.

    Attributes
    ----------
    times_s : ndarray
        Output times in seconds.
    states : ndarray
        N x 6 array containing position and velocity histories.
    raw_solution : object
        For ``solve_ivp`` propagation this is the SciPy solution object.
        For universal-anomaly propagation this is a small metadata dictionary.
    """

    times_s: np.ndarray
    states: np.ndarray
    raw_solution: Any

    @property
    def positions_km(self) -> np.ndarray:
        return self.states[:, :3]

    @property
    def velocities_km_s(self) -> np.ndarray:
        return self.states[:, 3:]

    @property
    def final_state(self) -> StateVector:
        return StateVector.from_vector6(self.states[-1])


def two_body_acceleration(
    r_km,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
) -> np.ndarray:
    """Point-mass gravitational acceleration in km/s^2."""
    r = np.asarray(r_km, dtype=float).reshape(3)
    rmag = np.linalg.norm(r)
    if rmag <= 0.0:
        raise ValueError("Position magnitude must be positive.")
    return -mu_km3_s2 * r / rmag**3


def j2_acceleration(
    r_km,
    *,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    radius_km: float = R_EARTH_KM,
    j2: float = J2_EARTH,
) -> np.ndarray:
    """J2 perturbing acceleration only, in km/s^2."""
    r = np.asarray(r_km, dtype=float).reshape(3)
    x, y, z = r
    rmag = np.linalg.norm(r)
    if rmag <= 0.0:
        raise ValueError("Position magnitude must be positive.")

    factor = 1.5 * j2 * mu_km3_s2 * radius_km**2 / rmag**5
    z2_over_r2 = z**2 / rmag**2
    return factor * np.array(
        [
            x * (5.0 * z2_over_r2 - 1.0),
            y * (5.0 * z2_over_r2 - 1.0),
            z * (5.0 * z2_over_r2 - 3.0),
        ]
    )


def two_body_rhs(
    t_s: float,
    state,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
) -> np.ndarray:
    """Six-state ODE for two-body motion."""
    state = np.asarray(state, dtype=float)
    return np.hstack(
        (
            state[3:],
            two_body_acceleration(state[:3], mu_km3_s2),
        )
    )


def j2_rhs(
    t_s: float,
    state,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    radius_km: float = R_EARTH_KM,
    j2: float = J2_EARTH,
) -> np.ndarray:
    """Six-state ODE for two-body gravity plus J2."""
    state = np.asarray(state, dtype=float)
    acceleration = two_body_acceleration(
        state[:3],
        mu_km3_s2,
    ) + j2_acceleration(
        state[:3],
        mu_km3_s2=mu_km3_s2,
        radius_km=radius_km,
        j2=j2,
    )
    return np.hstack((state[3:], acceleration))


def _state6(initial_state) -> np.ndarray:
    """Convert StateVector or six-value collection to a numeric state array."""
    if isinstance(initial_state, StateVector):
        return initial_state.vector6

    state = np.asarray(initial_state, dtype=float).reshape(-1)
    if state.size != 6:
        raise ValueError(
            "initial_state must contain six Cartesian state values."
        )
    return state


def _stumpff_c(z: float) -> float:
    """Universal-variable Stumpff C(z)."""
    if z > 1e-8:
        root_z = np.sqrt(z)
        return float((1.0 - np.cos(root_z)) / z)

    if z < -1e-8:
        root_neg_z = np.sqrt(-z)
        return float((np.cosh(root_neg_z) - 1.0) / (-z))

    # Series expansion near z = 0 avoids cancellation.
    return float(
        0.5
        - z / 24.0
        + z**2 / 720.0
        - z**3 / 40320.0
        + z**4 / 3628800.0
    )


def _stumpff_s(z: float) -> float:
    """Universal-variable Stumpff S(z)."""
    if z > 1e-8:
        root_z = np.sqrt(z)
        return float((root_z - np.sin(root_z)) / root_z**3)

    if z < -1e-8:
        root_neg_z = np.sqrt(-z)
        return float((np.sinh(root_neg_z) - root_neg_z) / root_neg_z**3)

    # Series expansion near z = 0 avoids cancellation.
    return float(
        1.0 / 6.0
        - z / 120.0
        + z**2 / 5040.0
        - z**3 / 362880.0
        + z**4 / 39916800.0
    )


def _initial_universal_anomaly_guess(
    r0_km: np.ndarray,
    v0_km_s: np.ndarray,
    delta_t_s: float,
    mu_km3_s2: float,
    alpha: float,
) -> float:
    """Choose a practical starting value for the universal anomaly."""
    if delta_t_s == 0.0:
        return 0.0

    sqrt_mu = np.sqrt(mu_km3_s2)
    r0_mag = np.linalg.norm(r0_km)

    # Elliptic orbit.
    if alpha > 1e-10:
        return float(sqrt_mu * alpha * delta_t_s)

    # Hyperbolic orbit. Use the standard logarithmic initial estimate when it
    # is well-defined, then fall back to a simple scale estimate if necessary.
    if alpha < -1e-10:
        sign_dt = np.sign(delta_t_s)
        rv_dot = float(np.dot(r0_km, v0_km_s))
        denominator = (
            rv_dot
            + sign_dt
            * np.sqrt(-mu_km3_s2 / alpha)
            * (1.0 - r0_mag * alpha)
        )
        numerator = -2.0 * mu_km3_s2 * alpha * delta_t_s

        if denominator != 0.0:
            log_argument = numerator / denominator
            if log_argument > 0.0:
                return float(
                    sign_dt
                    * np.sqrt(-1.0 / alpha)
                    * np.log(log_argument)
                )

        return float(sign_dt * sqrt_mu * abs(alpha) * abs(delta_t_s))

    # Near-parabolic orbit.
    return float(sqrt_mu * delta_t_s / r0_mag)


def universal_anomaly_state(
    initial_state,
    delta_t_s: float,
    *,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
) -> StateVector:
    """Propagate a two-body state by ``delta_t_s`` using universal variables.

    This solves the universal Kepler equation and evaluates the Lagrange
    ``f`` and ``g`` coefficients. It does not perform numerical time stepping,
    so ``rtol`` and ``atol`` are not part of the public interface.

    Parameters
    ----------
    initial_state : StateVector or six-value collection
        Cartesian state at the initial epoch.
    delta_t_s : float
        Propagation time relative to the initial epoch, in seconds. Negative
        values propagate backward.
    mu_km3_s2 : float, optional
        Gravitational parameter.
    """
    y0 = _state6(initial_state)
    r0 = y0[:3]
    v0 = y0[3:]

    dt = float(delta_t_s)
    mu = float(mu_km3_s2)
    if mu <= 0.0:
        raise ValueError("mu_km3_s2 must be positive.")

    if dt == 0.0:
        return StateVector(r0.copy(), v0.copy())

    r0_mag = np.linalg.norm(r0)
    if r0_mag <= 0.0:
        raise ValueError("Initial position magnitude must be positive.")

    v0_sq = float(np.dot(v0, v0))
    radial_velocity = float(np.dot(r0, v0) / r0_mag)
    sqrt_mu = np.sqrt(mu)

    # Reciprocal semi-major axis. Works for elliptic, parabolic, and
    # hyperbolic conics.
    alpha = 2.0 / r0_mag - v0_sq / mu

    chi = _initial_universal_anomaly_guess(
        r0,
        v0,
        dt,
        mu,
        alpha,
    )

    for _ in range(UNIVERSAL_MAX_ITERATIONS):
        z = alpha * chi**2
        c = _stumpff_c(z)
        s = _stumpff_s(z)

        residual = (
            (r0_mag * radial_velocity / sqrt_mu) * chi**2 * c
            + (1.0 - alpha * r0_mag) * chi**3 * s
            + r0_mag * chi
            - sqrt_mu * dt
        )

        derivative = (
            (r0_mag * radial_velocity / sqrt_mu)
            * chi
            * (1.0 - z * s)
            + (1.0 - alpha * r0_mag) * chi**2 * c
            + r0_mag
        )

        if derivative == 0.0:
            raise RuntimeError(
                "Universal-anomaly solver encountered a zero derivative."
            )

        correction = residual / derivative
        chi -= correction

        if abs(correction) <= UNIVERSAL_SOLVER_TOLERANCE * max(1.0, abs(chi)):
            break
    else:
        raise RuntimeError(
            "Universal-anomaly solver did not converge within "
            f"{UNIVERSAL_MAX_ITERATIONS} iterations."
        )

    z = alpha * chi**2
    c = _stumpff_c(z)
    s = _stumpff_s(z)

    f = 1.0 - (chi**2 / r0_mag) * c
    g = dt - (chi**3 / sqrt_mu) * s

    r = f * r0 + g * v0
    r_mag = np.linalg.norm(r)
    if r_mag <= 0.0:
        raise RuntimeError(
            "Universal-anomaly propagation produced zero position magnitude."
        )

    fdot = (
        sqrt_mu
        / (r_mag * r0_mag)
        * chi
        * (z * s - 1.0)
    )
    gdot = 1.0 - (chi**2 / r_mag) * c

    v = fdot * r0 + gdot * v0
    return StateVector(r, v)


def _fixed_time_grid(
    t_span_s: tuple[float, float],
    time_step_s: float,
) -> np.ndarray:
    """Create an inclusive fixed-spacing time grid for universal propagation."""
    t0, tf = map(float, t_span_s)
    step_mag = float(time_step_s)

    if step_mag <= 0.0:
        raise ValueError("time_step_s must be positive.")

    if t0 == tf:
        return np.array([t0], dtype=float)

    direction = 1.0 if tf > t0 else -1.0
    step = direction * step_mag
    duration = abs(tf - t0)
    full_steps = int(np.floor(duration / step_mag))

    times = t0 + step * np.arange(full_steps + 1, dtype=float)

    if not np.isclose(times[-1], tf, rtol=0.0, atol=1e-12):
        times = np.append(times, tf)
    else:
        times[-1] = tf

    return times


def propagate_universal_anomaly(
    initial_state,
    t_span_s: tuple[float, float],
    *,
    time_step_s: float,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
) -> PropagationResult:
    """Propagate a two-body orbit with the universal-anomaly formulation.

    ``time_step_s`` controls the returned sample spacing. Each requested state
    is propagated directly from the initial state, so error is not accumulated
    by repeatedly marching from one sample to the next.
    """
    y0 = _state6(initial_state)
    t0, tf = map(float, t_span_s)
    times = _fixed_time_grid((t0, tf), time_step_s)

    states = np.empty((len(times), 6), dtype=float)
    initial = StateVector.from_vector6(y0)

    for index, time_s in enumerate(times):
        propagated = universal_anomaly_state(
            initial,
            time_s - t0,
            mu_km3_s2=mu_km3_s2,
        )
        states[index] = propagated.vector6

    metadata = {
        "method": "universal_anomaly",
        "time_step_s": float(time_step_s),
        "mu_km3_s2": float(mu_km3_s2),
        "t_span_s": (t0, tf),
    }

    return PropagationResult(
        times_s=times,
        states=states,
        raw_solution=metadata,
    )


def altitude_event(
    threshold_km: float,
    *,
    direction: str = "below",
    radius_km: float = R_EARTH_KM,
) -> Callable[[float, np.ndarray, np.ndarray], bool]:
    """Build a stop event that triggers when altitude crosses ``threshold_km``.

    ``direction`` is ``"below"`` (e.g. reentry) or ``"above"``. Altitude is
    measured as ``|r| - radius_km`` (spherical Earth).
    """
    if direction not in ("below", "above"):
        raise ValueError("direction must be 'below' or 'above'.")
    threshold = float(threshold_km)

    def event(t_s: float, r_km: np.ndarray, v_km_s: np.ndarray) -> bool:
        altitude = float(np.linalg.norm(r_km)) - radius_km
        return altitude <= threshold if direction == "below" else altitude >= threshold

    return event


def enckes_method(
    initial_state,
    t_span_s: tuple[float, float],
    time_step_s: float,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    perturbations: Iterable[Callable[[float, np.ndarray, np.ndarray], np.ndarray]] | None = None,
    *,
    include_perturbations: bool = True,
    rectification_threshold: float | None = 0.1,
    events: Iterable[Callable[[float, np.ndarray, np.ndarray], bool]] | None = None,
) -> PropagationResult:
    """Propagate an orbit with fixed-step Encke updates.

    Deviation is actual state minus the two-body reference state. Each
    perturbation is a callable ``acceleration(t_s, r_km, v_km_s)`` that
    returns an acceleration vector in km/s^2.

    Parameters
    ----------
    initial_state : StateVector or six-value collection
        Initial Cartesian state.
    t_span_s : (float, float)
        Initial and final propagation times in seconds.
    time_step_s : float
        Positive output and integration step, in seconds.
    mu_km3_s2 : float, optional
        Gravitational parameter for the reference orbit and Encke equations.
    perturbations : iterable of callable, optional
        Perturbing acceleration models. Each receives time in seconds,
        position in km, and velocity in km/s, and returns acceleration in
        km/s^2.
    include_perturbations : bool, optional
        Set false to skip all perturbation acceleration models.
    rectification_threshold : float or None, optional
        Reset the reference orbit when the deviation magnitude divided by
        reference radius exceeds this value. Use None to disable rectification.
    events : iterable of callable, optional
        Stop events ``event(t_s, r_km, v_km_s) -> bool`` checked after each
        step on the actual state (see :func:`altitude_event`). Propagation
        ends at the first step where any returns True; that step is the last
        row of the result and ``raw_solution["terminated_by_event"]`` is True.

    Returns
    -------
    PropagationResult
        Cartesian state history on the fixed time grid.
    """
    mu = float(mu_km3_s2)
    stop_events = tuple(events) if events is not None else ()
    if any(not callable(event) for event in stop_events):
        raise TypeError("Every event must be callable.")
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu_km3_s2 must be finite and positive.")
    if rectification_threshold is not None and (
        not np.isfinite(rectification_threshold) or rectification_threshold <= 0.0
    ):
        raise ValueError("rectification_threshold must be positive or None.")

    acceleration_models = (
        tuple(perturbations) if perturbations is not None else ()
    )
    if any(not callable(model) for model in acceleration_models):
        raise TypeError("Every perturbation must be callable.")
    if not include_perturbations:
        acceleration_models = ()

    y0 = _state6(initial_state)
    t0, tf = map(float, t_span_s)
    times = _fixed_time_grid((t0, tf), time_step_s)
    states = np.empty((len(times), 6), dtype=float)
    reference = StateVector.from_vector6(y0)
    reference_epoch_s = t0
    delta_r = np.zeros(3, dtype=float)
    delta_v = np.zeros(3, dtype=float)

    def reference_at(time_s: float) -> StateVector:
        return universal_anomaly_state(
            reference,
            time_s - reference_epoch_s,
            mu_km3_s2=mu,
        )

    states[0] = y0
    last_index = len(times) - 1
    terminated = False
    start_reference = reference_at(t0)
    for index in range(1, len(times)):
        step_start_s = times[index - 1]
        step_end_s = times[index]
        step_s = step_end_s - step_start_s
        r_osc = start_reference.r_km
        r_osc_mag = np.linalg.norm(r_osc)
        actual_r = r_osc + delta_r
        actual_v = start_reference.v_km_s + delta_v

        q_value = float(
            np.dot(delta_r, delta_r + 2.0 * r_osc) / r_osc_mag**2
        )
        if q_value <= -1.0:
            raise RuntimeError(
                "Encke deviation places the object at or through Earth's center."
            )
        f_q = q_value * (q_value**2 + 3.0 * q_value + 3.0) / (
            1.0 + (1.0 + q_value) ** 1.5
        )
        delta_acceleration = -mu / r_osc_mag**3 * (
            delta_r - f_q * r_osc
        )

        for acceleration_model in acceleration_models:
            contribution = np.asarray(
                acceleration_model(step_start_s, actual_r.copy(), actual_v.copy()),
                dtype=float,
            ).reshape(-1)
            if contribution.size != 3 or not np.all(np.isfinite(contribution)):
                raise ValueError(
                    "Each perturbation must return three finite "
                    "acceleration components in km/s^2."
                )
            delta_acceleration += contribution

        delta_r = (
            delta_r
            + delta_v * step_s
            + 0.5 * delta_acceleration * step_s**2
        )
        delta_v = delta_v + delta_acceleration * step_s
        propagated_reference = reference_at(step_end_s)
        actual_r = propagated_reference.r_km + delta_r
        actual_v = propagated_reference.v_km_s + delta_v
        states[index, :3] = actual_r
        states[index, 3:] = actual_v
        start_reference = propagated_reference

        if any(event(step_end_s, actual_r, actual_v) for event in stop_events):
            last_index = index
            terminated = True
            break

        if (
            rectification_threshold is not None
            and np.linalg.norm(delta_r) / np.linalg.norm(propagated_reference.r_km)
            >= rectification_threshold
        ):
            reference = StateVector(actual_r.copy(), actual_v.copy())
            reference_epoch_s = step_end_s
            start_reference = reference
            delta_r = np.zeros(3, dtype=float)
            delta_v = np.zeros(3, dtype=float)

    metadata = {
        "method": "encke",
        "time_step_s": float(time_step_s),
        "mu_km3_s2": mu,
        "t_span_s": (t0, tf),
        "perturbations_enabled": bool(acceleration_models),
        "rectification_threshold": rectification_threshold,
        "terminated_by_event": terminated,
    }
    return PropagationResult(
        times_s=times[: last_index + 1],
        states=states[: last_index + 1],
        raw_solution=metadata,
    )


def altitude_crossing_event(
    threshold_km: float,
    *,
    direction: str = "below",
    radius_km: float = R_EARTH_KM,
) -> Callable[[float, np.ndarray], float]:
    """Build a terminal ``solve_ivp`` event for an altitude crossing.

    For use with :func:`cowells_method`. ``direction`` is ``"below"`` or
    ``"above"``. Altitude is ``|r| - radius_km`` (spherical Earth).
    """
    if direction not in ("below", "above"):
        raise ValueError("direction must be 'below' or 'above'.")
    threshold = float(threshold_km)

    def event(t_s: float, y: np.ndarray, *args) -> float:
        return float(np.linalg.norm(y[:3])) - radius_km - threshold

    event.terminal = True
    event.direction = -1.0 if direction == "below" else 1.0
    return event


def cowells_method(
    initial_state,
    t_span_s: tuple[float, float],
    time_step_s: float,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    perturbations: Iterable[Callable[[float, np.ndarray, np.ndarray], np.ndarray]] | None = None,
    *,
    include_perturbations: bool = True,
    events: Iterable[Callable] | None = None,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    method: str = "DOP853",
) -> PropagationResult:
    """Propagate an orbit by direct integration of Cowell's equation.

    Integrates ``r_ddot = -mu / r^3 * r + a_perts`` in ECI. Each perturbation
    is a callable ``acceleration(t_s, r_km, v_km_s)`` returning km/s^2, the
    same interface used by :func:`enckes_method`.

    Parameters
    ----------
    initial_state : StateVector or six-value collection
        Initial Cartesian state.
    t_span_s : (float, float)
        Initial and final propagation times in seconds.
    time_step_s : float
        Positive output spacing in seconds. The integrator step is adaptive.
    perturbations : iterable of callable, optional
        Perturbing acceleration models.
    include_perturbations : bool, optional
        Set false to skip all perturbation models.
    events : iterable of callable, optional
        Native ``solve_ivp`` events ``event(t, y)``, e.g.
        :func:`altitude_crossing_event`. Terminal events end propagation at
        the root; the event state is appended as the last row.
        ``raw_solution.status == 1`` indicates a terminal event fired.
    rtol, atol, method : optional
        ``solve_ivp`` settings.

    Returns
    -------
    PropagationResult
        State history on the output grid; ``raw_solution`` is the SciPy solution.
    """
    mu = float(mu_km3_s2)
    if not np.isfinite(mu) or mu <= 0.0:
        raise ValueError("mu_km3_s2 must be finite and positive.")

    acceleration_models = tuple(perturbations) if perturbations is not None else ()
    if any(not callable(model) for model in acceleration_models):
        raise TypeError("Every perturbation must be callable.")
    if not include_perturbations:
        acceleration_models = ()

    def rhs(t_s: float, y: np.ndarray) -> np.ndarray:
        r = y[:3]
        v = y[3:]
        acceleration = -mu / np.linalg.norm(r) ** 3 * r
        for model in acceleration_models:
            acceleration = acceleration + np.asarray(
                model(t_s, r.copy(), v.copy()), dtype=float
            ).reshape(3)
        return np.concatenate((v, acceleration))

    y0 = _state6(initial_state)
    t0, tf = map(float, t_span_s)
    times = _fixed_time_grid((t0, tf), time_step_s)
    event_list = list(events) if events is not None else None

    result = propagate_custom(
        y0,
        (t0, tf),
        rhs,
        times_s=times,
        rtol=rtol,
        atol=atol,
        method=method,
        events=event_list,
    )
    solution = result.raw_solution
    if solution.status == 1:
        for t_event, y_event in zip(solution.t_events, solution.y_events):
            if len(t_event):
                result = PropagationResult(
                    times_s=np.append(result.times_s, t_event[0]),
                    states=np.vstack((result.states, y_event[0])),
                    raw_solution=solution,
                )
                break
    return result


def propagate_custom(
    initial_state,
    t_span_s: tuple[float, float],
    rhs: Callable,
    *,
    times_s: Iterable[float] | None = None,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    method: str = "DOP853",
    args: tuple = (),
    allow_any_size: bool = False,
    **solve_ivp_kwargs,
) -> PropagationResult:
    """Generic ``solve_ivp`` wrapper for a six-state orbital ODE."""
    y0 = np.asarray(initial_state.vector6 if isinstance(initial_state, StateVector) else initial_state, dtype=float).reshape(-1)
    if y0.size != 6 and not allow_any_size:
        raise ValueError(
            "initial_state must contain six Cartesian state values."
        )
    t_eval = None if times_s is None else np.asarray(list(times_s), dtype=float)

    solution = solve_ivp(
        rhs,
        tuple(map(float, t_span_s)),
        y0,
        t_eval=t_eval,
        rtol=rtol,
        atol=atol,
        method=method,
        args=args,
        **solve_ivp_kwargs,
    )
    if not solution.success:
        raise RuntimeError(f"solve_ivp failed: {solution.message}")

    return PropagationResult(
        times_s=solution.t.copy(),
        states=solution.y.T.copy(),
        raw_solution=solution,
    )


def propagate_two_body(
    initial_state,
    t_span_s: tuple[float, float],
    *,
    times_s: Iterable[float] | None = None,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    method: str = "DOP853",
    time_step_s: float | None = None,
    **solve_ivp_kwargs,
) -> PropagationResult:
    """Propagate an ECI Cartesian state using point-mass gravity.

    Parameters
    ----------
    initial_state : StateVector or six-value collection
        Initial Cartesian state.
    t_span_s : (float, float)
        Initial and final propagation times in seconds.
    times_s : iterable of float, optional
        Requested output times for ``solve_ivp`` methods.
    mu_km3_s2 : float, optional
        Gravitational parameter.
    rtol, atol : float, optional
        Numerical integration tolerances for SciPy ``solve_ivp`` methods.
        They are not used by the universal-anomaly method.
    method : str, optional
        Any SciPy ``solve_ivp`` method, such as ``DOP853``, ``RK45``,
        ``RK23``, ``Radau``, ``BDF``, or ``LSODA``. Also accepts
        ``"universal_anomaly"`` (and the aliases ``"universal"`` and
        ``"universal_variable"``).
    time_step_s : float, optional
        Required when ``method="universal_anomaly"``. Sets the returned time
        spacing and replaces the need to specify ``rtol`` and ``atol``.
    **solve_ivp_kwargs
        Additional SciPy options for numerical integration methods only.
    """
    normalized_method = (
        str(method)
        .strip()
        .lower()
        .replace("-", "_")
        .replace(" ", "_")
    )

    universal_methods = {
        "universal",
        "universal_anomaly",
        "universal_variable",
        "universal_variables",
    }

    if normalized_method in universal_methods:
        if time_step_s is None:
            raise ValueError(
                "method='universal_anomaly' requires time_step_s."
            )
        if times_s is not None:
            raise ValueError(
                "Universal-anomaly propagation uses time_step_s instead of "
                "times_s."
            )
        if solve_ivp_kwargs:
            names = ", ".join(sorted(solve_ivp_kwargs))
            raise ValueError(
                "solve_ivp keyword arguments are not used by the "
                f"universal-anomaly method: {names}"
            )

        return propagate_universal_anomaly(
            initial_state,
            t_span_s,
            time_step_s=time_step_s,
            mu_km3_s2=mu_km3_s2,
        )

    if time_step_s is not None:
        raise ValueError(
            "time_step_s is only used when method='universal_anomaly'."
        )

    return propagate_custom(
        initial_state,
        t_span_s,
        two_body_rhs,
        times_s=times_s,
        rtol=rtol,
        atol=atol,
        method=method,
        args=(mu_km3_s2,),
        **solve_ivp_kwargs,
    )


def propagate_j2(
    initial_state,
    t_span_s: tuple[float, float],
    *,
    times_s: Iterable[float] | None = None,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    radius_km: float = R_EARTH_KM,
    j2: float = J2_EARTH,
    rtol: float = DEFAULT_RTOL,
    atol: float = DEFAULT_ATOL,
    method: str = "DOP853",
    **solve_ivp_kwargs,
) -> PropagationResult:
    """Propagate an ECI Cartesian state using point-mass gravity plus J2."""
    return propagate_custom(
        initial_state,
        t_span_s,
        j2_rhs,
        times_s=times_s,
        rtol=rtol,
        atol=atol,
        method=method,
        args=(mu_km3_s2, radius_km, j2),
        **solve_ivp_kwargs,
    )


__all__ = [
    "PropagationResult",
    "two_body_acceleration",
    "j2_acceleration",
    "two_body_rhs",
    "j2_rhs",
    "universal_anomaly_state",
    "propagate_universal_anomaly",
    "altitude_event",
    "altitude_crossing_event",
    "cowells_method",
    "enckes_method",
    "propagate_custom",
    "propagate_two_body",
    "propagate_j2",
]
