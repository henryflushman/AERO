"""Relative-motion utilities for chief/deputy spacecraft rendezvous analysis.

This module is designed for the Intermediate Orbits course and follows the
standard LVLH/RTN framing used in the rest of the toolkit:

- +x/radial: away from the central body
- +y/along-track: tangent to the chief orbit
- +z/cross-track: orbit normal

The key distinction implemented here is that a relative velocity must include the
frame-rotation term when converting a relative state into the rotating LVLH
frame. A pure DCM rotation is not sufficient.
"""

from __future__ import annotations

import warnings

import numpy as np

try:
    from .Constants import MU_EARTH_KM3_S2
    from .OrbitalFrames import eci_to_lvlh_dcm
except ImportError:  # pragma: no cover
    from Constants import MU_EARTH_KM3_S2
    from OrbitalFrames import eci_to_lvlh_dcm


def _vec3(value, name: str) -> np.ndarray:
    """Return a 3-vector as a float array."""
    vector = np.asarray(value, dtype=float).reshape(-1)
    if vector.size != 3:
        raise ValueError(f"{name} must contain exactly three values.")
    return vector


def mean_motion_from_orbit_radius(radius_km: float) -> float:
    """Return mean motion n in rad/s for a circular orbit of radius km."""
    radius = float(radius_km)
    if radius <= 0.0:
        raise ValueError("radius_km must be positive.")
    return np.sqrt(MU_EARTH_KM3_S2 / radius**3)


def cwhill_state_transition_matrix(elapsed_s: float, mean_motion_rad_s: float) -> np.ndarray:
    """Return the 6x6 CW state-transition matrix for a circular chief orbit.

    The state vector is ordered as:
    [x, y, z, xdot, ydot, zdot]^T
    where x/y/z are in the chief LVLH frame.
    """
    t = float(elapsed_s)
    n = float(mean_motion_rad_s)
    if n <= 0.0:
        raise ValueError("mean_motion_rad_s must be positive.")

    c = np.cos(n * t)
    s = np.sin(n * t)

    Phi = np.array(
        [
            [4.0 - 3.0 * c, 0.0, 0.0, s / n, 2.0 * (1.0 - c) / n, 0.0],
            [6.0 * (s - n * t), 1.0, 0.0,
             2.0 * (c - 1.0) / n, (4.0 * s - 3.0 * n * t) / n, 0.0],
            [0.0, 0.0, c, 0.0, 0.0, s / n],
            [3.0 * n * s, 0.0, 0.0, c, 2.0 * s, 0.0],
            [6.0 * n * (c - 1.0), 0.0, 0.0,
             -2.0 * s, 4.0 * c - 3.0, 0.0],
            [0.0, 0.0, -n * s, 0.0, 0.0, c],
        ],
        dtype=float,
    )
    return Phi

def cwhill_transition_matrix(elapsed_s: float, mean_motion_rad_s: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    s = np.sin(mean_motion_rad_s * elapsed_s)
    c = np.cos(mean_motion_rad_s * elapsed_s)
    phi_rr = np.array(
        [
            [4.0 - 3.0*c, 0.0, 0.0],
            [6.0 * (s - mean_motion_rad_s * elapsed_s), 1.0, 0.0],
            [0.0, 0.0, c],
        ]
    )
    phi_rv = np.array(
        [
            [1.0/(mean_motion_rad_s)*s, 2.0/(mean_motion_rad_s)*(1.0 - c), 0.0],
            [2.0/(mean_motion_rad_s)*(c-1.0), (4.0*s - 3.0*mean_motion_rad_s*elapsed_s)/(mean_motion_rad_s), 0.0],
            [0.0, 0.0, 1.0/(mean_motion_rad_s)*s]
        ]
    )
    phi_vr = np.array(
        [
            [3.0*mean_motion_rad_s*s, 0.0, 0.0],
            [6.0*mean_motion_rad_s*(c-1.0), 0.0, 0.0],
            [0.0, 0.0, -mean_motion_rad_s*s]
        ]
    )
    phi_vv = np.array(
        [
            [c, 2.0*s, 0.0],
            [-2.0*s, 4.0*c - 3.0, 0.0],
            [0.0, 0.0, c]
        ]
    )
    return phi_rr, phi_rv, phi_vr, phi_vv


def _warn_if_cw_assumptions_are_stretched(delta_r0_km, delta_v0_km_s, eccentricity: float | None = None) -> None:
    """Warn when the CW linearization is likely outside its valid regime."""
    r_mag = float(np.linalg.norm(_vec3(delta_r0_km, "delta_r0_km")))
    if r_mag > 1.0:
        warnings.warn(
            "The Clohessy-Wiltshire equations assume a small relative separation. "
            f"Initial separation magnitude is {r_mag:.3f} km, which may be outside the linear regime.",
            RuntimeWarning,
            stacklevel=2,
        )
    if eccentricity is not None and abs(float(eccentricity)) > 0.05:
        warnings.warn(
            "The CW equations are derived for a near-circular chief orbit. "
            f"Received eccentricity e={eccentricity:.3f}; the approximation may be noticeably inaccurate.",
            RuntimeWarning,
            stacklevel=2,
        )


def cwhill_delr_delv(
    delta_r0_km,
    delta_v0_km_s,
    elapsed_s: float,
    mean_motion_rad_s: float,
    *,
    eccentricity: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate a chief-relative state using the CW transition matrix.

    This formulation writes the initial state as

        x0 = [delr0; delv0]

    and propagates it by

        x(t) = Phi(t) @ x0

    where ``Phi`` is the Clohessy-Wiltshire state-transition matrix for a
    circular chief orbit.

    A warning is raised when the initial relative position is large or the chief
    orbit eccentricity is large enough that the linearized model may be invalid.
    """
    _warn_if_cw_assumptions_are_stretched(delta_r0_km, delta_v0_km_s, eccentricity=eccentricity)
    r0 = _vec3(delta_r0_km, "delta_r0_km")
    v0 = _vec3(delta_v0_km_s, "delta_v0_km_s")
    x0 = np.concatenate((r0, v0))

    Phi = cwhill_state_transition_matrix(elapsed_s, mean_motion_rad_s)
    x1 = Phi @ x0
    return x1[:3], x1[3:]


def cwhill_propagate(
    delta_r0_km,
    delta_v0_km_s,
    elapsed_s: float,
    mean_motion_rad_s: float,
    *,
    eccentricity: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Propagate a small relative state using the CW equations.

    Parameters
    ----------
    delta_r0_km : array-like of shape (3,)
        Initial relative position in the chief LVLH frame, km.
    delta_v0_km_s : array-like of shape (3,)
        Initial relative velocity in the chief LVLH frame, km/s.
    elapsed_s : float
        Propagation time, s.
    mean_motion_rad_s : float
        Chief orbital mean motion, rad/s.
    eccentricity : float, optional
        Chief orbit eccentricity. Used only to emit a warning when the chief orbit
        is not sufficiently circular for the CW approximation.

    Returns
    -------
    delta_r : np.ndarray, shape (3,)
        Relative position after the elapsed time, km.
    delta_v : np.ndarray, shape (3,)
        Relative velocity after the elapsed time, km/s.
    """
    return cwhill_delr_delv(
        delta_r0_km,
        delta_v0_km_s,
        elapsed_s,
        mean_motion_rad_s,
        eccentricity=eccentricity,
    )


def exact_relative_motion(
    chief_r_eci_km,
    chief_v_eci_km_s,
    deputy_r_eci_km,
    deputy_v_eci_km_s,
) -> tuple[np.ndarray, np.ndarray]:
    """Return exact relative motion in LVLH without the CW approximation."""
    rc = _vec3(chief_r_eci_km, "chief_r_eci_km")
    vc = _vec3(chief_v_eci_km_s, "chief_v_eci_km_s")
    rho_eci = _vec3(deputy_r_eci_km, "deputy_r_eci_km") - rc
    rho_dot_eci = _vec3(deputy_v_eci_km_s, "deputy_v_eci_km_s") - vc

    C = eci_to_lvlh_dcm(rc, vc)
    rho_lvlh = C @ rho_eci
    h_mag = np.linalg.norm(np.cross(rc, vc))
    omega_lvlh = np.array([0.0, 0.0, h_mag / np.dot(rc, rc)])
    rho_dot_lvlh = C @ rho_dot_eci - np.cross(omega_lvlh, rho_lvlh)
    return rho_lvlh, rho_dot_lvlh


def relative_separation_norm(delta_r_km) -> float:
    """Return the Euclidean norm of the relative position vector."""
    return float(np.linalg.norm(_vec3(delta_r_km, "delta_r_km")))


__all__ = [
    "mean_motion_from_orbit_radius",
    "cwhill_state_transition_matrix",
    "cwhill_delr_delv",
    "cwhill_propagate",
    "exact_relative_motion",
    "relative_separation_norm",
]
