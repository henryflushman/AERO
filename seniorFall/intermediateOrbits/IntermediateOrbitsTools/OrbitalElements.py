"""Classical orbital elements and Cartesian state-vector conversions.

Conventions
-----------
Distance : km
Time     : s
Velocity : km/s
Angles   : radians internally

``ClassicalOrbitalElements`` accepts exactly six independent orbital quantities,
converts them to the canonical set ``a, e, i, RAAN, argp, nu``, and exposes
other equivalent quantities as automatically calculated attributes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field

import numpy as np

try:
    from .Constants import MU_EARTH_KM3_S2
    from .OrbitalFrames import perifocal_to_eci_dcm
except ImportError:
    from Constants import MU_EARTH_KM3_S2
    from OrbitalFrames import perifocal_to_eci_dcm


TWOPI = 2.0 * np.pi
DEFAULT_COE_ORDER = ("a", "e", "i", "raan", "argp", "nu")


def wrap_2pi(angle_rad: float) -> float:
    return float(angle_rad % TWOPI)


def _vec3(value, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float).reshape(-1)
    if vector.size != 3:
        raise ValueError(f"{name} must contain exactly three values.")
    return vector


def _angle_to_rad(value: float, unit: str) -> float:
    unit = unit.strip().lower()
    if unit in {"rad", "radian", "radians"}:
        return float(value)
    if unit in {"deg", "degree", "degrees"}:
        return float(np.deg2rad(value))
    raise ValueError("angle_unit must be 'rad' or 'deg'.")


def _parse_element_name(name: str) -> tuple[str, str | None]:
    raw = str(name).strip()
    unit = None

    if raw.lower().endswith("_deg"):
        raw, unit = raw[:-4], "deg"
    elif raw.lower().endswith("_rad"):
        raw, unit = raw[:-4], "rad"

    if raw == "E":
        return "E", unit
    if raw == "M":
        return "M", unit
    if raw == "Ω":
        return "raan", unit
    if raw == "ω":
        return "argp", unit
    if raw in {"ν", "θ"}:
        return "nu", unit

    key = raw.lower().replace("_", "").replace("-", "").replace(" ", "")
    for suffix in ("km2s", "km"):
        if key.endswith(suffix):
            key = key[: -len(suffix)]
            break

    aliases = {
        "a": {"a", "sma", "semimajoraxis"},
        "h": {"h", "specificangularmomentum"},
        "p": {"p", "semilatusrectum"},
        "rp": {"rp", "periapsis", "periapsisradius", "perigee", "perigeeradius"},
        "ra": {"ra", "apoapsis", "apoapsisradius", "apogee", "apogeeradius"},
        "e": {"e", "ecc", "eccentricity"},
        "i": {"i", "inc", "inclination"},
        "raan": {"raan", "rightascensionofascendingnode"},
        "argp": {"argp", "omega", "argumentofperiapsis", "argumentofperigee"},
        "nu": {"nu", "theta", "trueanomaly"},
        "E": {"eccentricanomaly"},
        "M": {"meananomaly"},
    }

    for canonical, choices in aliases.items():
        if key in choices:
            return canonical, unit
    raise KeyError(f"Unsupported classical orbital element name: {name!r}")


def eccentric_anomaly_from_true(nu_rad: float, e: float) -> float:
    denominator = 1.0 + e * np.cos(nu_rad)
    sin_E = np.sqrt(1.0 - e**2) * np.sin(nu_rad) / denominator
    cos_E = (e + np.cos(nu_rad)) / denominator
    return wrap_2pi(np.arctan2(sin_E, cos_E))


def true_anomaly_from_eccentric(E_rad: float, e: float) -> float:
    denominator = 1.0 - e * np.cos(E_rad)
    sin_nu = np.sqrt(1.0 - e**2) * np.sin(E_rad) / denominator
    cos_nu = (np.cos(E_rad) - e) / denominator
    return wrap_2pi(np.arctan2(sin_nu, cos_nu))


def mean_anomaly_from_eccentric(E_rad: float, e: float) -> float:
    return wrap_2pi(E_rad - e * np.sin(E_rad))


def eccentric_anomaly_from_mean(
    M_rad: float,
    e: float,
    *,
    tolerance: float = 1e-13,
    max_iterations: int = 50,
) -> float:
    M = wrap_2pi(M_rad)
    E = M if e < 0.8 else np.pi

    for _ in range(max_iterations):
        step = (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
        E -= step
        if abs(step) <= tolerance:
            return wrap_2pi(E)
    raise RuntimeError("Kepler equation did not converge.")


@dataclass(frozen=True)
class StateVector:
    """Cartesian orbital state using km and km/s."""

    r_km: np.ndarray
    v_km_s: np.ndarray

    def __post_init__(self) -> None:
        object.__setattr__(self, "r_km", _vec3(self.r_km, "r_km"))
        object.__setattr__(self, "v_km_s", _vec3(self.v_km_s, "v_km_s"))

    @property
    def vector6(self) -> np.ndarray:
        return np.hstack((self.r_km, self.v_km_s))

    @classmethod
    def from_vector6(cls, state) -> "StateVector":
        state = np.asarray(state, dtype=float).reshape(-1)
        if state.size != 6:
            raise ValueError("state must contain six values [rx, ry, rz, vx, vy, vz].")
        return cls(state[:3], state[3:])


@dataclass(frozen=True, init=False)
class ClassicalOrbitalElements:
    """Complete elliptical classical orbital elements.

    Input may be a six-value sequence, mapping, six positional values, or six
    keyword values. Use ``names=...`` when a sequence is not in the default
    order ``[a, e, i, raan, argp, nu]``.

    Supported equivalent inputs:
        size : a+e, h+e, p+e, or rp+ra
        phase: nu/theta, E/eccentric_anomaly, or M/mean_anomaly
    """

    a_km: float
    e: float
    i_rad: float
    raan_rad: float
    argp_rad: float
    nu_rad: float
    _mu_km3_s2: float = field(repr=False, compare=False)

    def __init__(
        self,
        *elements,
        names: Sequence[str] | None = None,
        angle_unit: str = "rad",
        mu_km3_s2: float = MU_EARTH_KM3_S2,
        **named_elements,
    ) -> None:
        raw = self._prepare_input(elements, named_elements, names)
        canonical = self._canonicalize(raw, angle_unit, float(mu_km3_s2))

        for name, value in canonical.items():
            object.__setattr__(self, name, value)
        object.__setattr__(self, "_mu_km3_s2", float(mu_km3_s2))

        if not 0.0 <= self.i_rad <= np.pi:
            raise ValueError("Inclination must be between 0 and 180 degrees.")
        if self._mu_km3_s2 <= 0.0:
            raise ValueError("mu_km3_s2 must be positive.")

    @staticmethod
    def _sequence_to_mapping(values, names: Sequence[str] | None) -> dict:
        values = list(values)
        if len(values) != 6:
            raise ValueError(f"A COE collection must contain exactly six values; received {len(values)}.")

        input_names = DEFAULT_COE_ORDER if names is None else tuple(names)
        if len(input_names) != 6:
            raise ValueError("names must contain exactly six element names.")
        if len(set(input_names)) != 6:
            raise ValueError("names must contain six unique element names.")
        return dict(zip(input_names, values))

    @classmethod
    def _prepare_input(cls, positional, named_elements, names) -> dict:
        if positional and named_elements:
            raise TypeError("Provide COEs either positionally/as a collection or as keywords, not both.")

        if named_elements:
            if names is not None:
                raise TypeError("names cannot be used with keyword orbital elements.")
            raw = dict(named_elements)
        elif len(positional) == 1 and isinstance(positional[0], Mapping):
            if names is not None:
                raise TypeError("names is unnecessary when a mapping is supplied.")
            raw = dict(positional[0])
        elif len(positional) == 1 and not isinstance(positional[0], (str, bytes)):
            raw = cls._sequence_to_mapping(positional[0], names)
        elif len(positional) == 6:
            raw = cls._sequence_to_mapping(positional, names)
        elif len(positional) == 0:
            raise ValueError("ClassicalOrbitalElements requires exactly six orbital quantities.")
        else:
            raise ValueError("Provide six orbital values, or one mapping/collection containing six values.")

        if len(raw) != 6:
            raise ValueError(
                f"ClassicalOrbitalElements requires exactly six independent quantities; received {len(raw)}."
            )
        return raw

    @staticmethod
    def _get_angle(values, name: str, default_unit: str) -> float:
        if name not in values:
            raise ValueError(f"Missing required orbital angle: {name}.")
        value, explicit_unit = values[name]
        return _angle_to_rad(value, explicit_unit or default_unit)

    @classmethod
    def _canonicalize(cls, raw, angle_unit: str, mu: float) -> dict[str, float]:
        values: dict[str, tuple[float, str | None]] = {}

        for raw_name, raw_value in raw.items():
            name, explicit_unit = _parse_element_name(raw_name)
            if name in values:
                raise ValueError(f"{raw_name!r} duplicates another supplied orbital element.")
            value = float(raw_value)
            if not np.isfinite(value):
                raise ValueError(f"{raw_name!r} must be finite.")
            values[name] = (value, explicit_unit)

        if "rp" in values or "ra" in values:
            if "rp" not in values or "ra" not in values:
                raise ValueError("rp and ra must be supplied together.")
            if any(name in values for name in ("a", "h", "p", "e")):
                raise ValueError("Do not combine rp/ra with a, h, p, or e.")
            rp, ra = values["rp"][0], values["ra"][0]
            if rp <= 0.0 or ra < rp:
                raise ValueError("Require 0 < rp <= ra.")
            a = 0.5 * (rp + ra)
            e = (ra - rp) / (ra + rp)
        else:
            if "e" not in values:
                raise ValueError("e is required unless rp and ra are supplied.")
            e = values["e"][0]
            if not 0.0 <= e < 1.0:
                raise ValueError("This helper supports elliptical orbits with 0 <= e < 1.")

            size_names = [name for name in ("a", "h", "p") if name in values]
            if len(size_names) != 1:
                raise ValueError("Provide exactly one of a, h, or p.")

            name = size_names[0]
            value = values[name][0]
            if value <= 0.0:
                raise ValueError(f"{name} must be positive.")

            if name == "a":
                a = value
            elif name == "h":
                a = (value**2 / mu) / (1.0 - e**2)
            else:
                a = value / (1.0 - e**2)

        i = cls._get_angle(values, "i", angle_unit)
        raan = cls._get_angle(values, "raan", angle_unit)
        argp = cls._get_angle(values, "argp", angle_unit)

        phase_names = [name for name in ("nu", "E", "M") if name in values]
        if len(phase_names) != 1:
            raise ValueError("Provide exactly one of nu/theta, E, or M.")

        phase_name = phase_names[0]
        phase = cls._get_angle(values, phase_name, angle_unit)
        if phase_name == "nu":
            nu = wrap_2pi(phase)
        elif phase_name == "E":
            nu = true_anomaly_from_eccentric(phase, e)
        else:
            nu = true_anomaly_from_eccentric(eccentric_anomaly_from_mean(phase, e), e)

        return {
            "a_km": float(a),
            "e": float(e),
            "i_rad": float(i),
            "raan_rad": wrap_2pi(raan),
            "argp_rad": wrap_2pi(argp),
            "nu_rad": wrap_2pi(nu),
        }

    @classmethod
    def from_degrees(
        cls,
        a_km: float,
        e: float,
        i_deg: float,
        raan_deg: float,
        argp_deg: float,
        nu_deg: float,
        *,
        mu_km3_s2: float = MU_EARTH_KM3_S2,
    ) -> "ClassicalOrbitalElements":
        return cls(
            a_km=a_km,
            e=e,
            i_deg=i_deg,
            raan_deg=raan_deg,
            argp_deg=argp_deg,
            nu_deg=nu_deg,
            mu_km3_s2=mu_km3_s2,
        )

    @property
    def mu_km3_s2(self) -> float:
        return self._mu_km3_s2

    @property
    def inclination_deg(self) -> float:
        return float(np.rad2deg(self.i_rad))

    @property
    def raan_deg(self) -> float:
        return float(np.rad2deg(self.raan_rad))

    @property
    def argument_of_periapsis_deg(self) -> float:
        return float(np.rad2deg(self.argp_rad))

    @property
    def true_anomaly_deg(self) -> float:
        return float(np.rad2deg(self.nu_rad))

    @property
    def eccentric_anomaly_rad(self) -> float:
        return eccentric_anomaly_from_true(self.nu_rad, self.e)

    @property
    def eccentric_anomaly_deg(self) -> float:
        return float(np.rad2deg(self.eccentric_anomaly_rad))

    @property
    def mean_anomaly_rad(self) -> float:
        return mean_anomaly_from_eccentric(self.eccentric_anomaly_rad, self.e)

    @property
    def mean_anomaly_deg(self) -> float:
        return float(np.rad2deg(self.mean_anomaly_rad))

    @property
    def p_km(self) -> float:
        return self.a_km * (1.0 - self.e**2)

    @property
    def h_km2_s(self) -> float:
        return float(np.sqrt(self.mu_km3_s2 * self.p_km))

    @property
    def periapsis_radius_km(self) -> float:
        return self.a_km * (1.0 - self.e)

    @property
    def apoapsis_radius_km(self) -> float:
        return self.a_km * (1.0 + self.e)

    @property
    def radius_km(self) -> float:
        return self.p_km / (1.0 + self.e * np.cos(self.nu_rad))

    @property
    def mean_motion_rad_s(self) -> float:
        return float(np.sqrt(self.mu_km3_s2 / self.a_km**3))

    @property
    def mean_motion_rev_day(self) -> float:
        return self.mean_motion_rad_s * 86400.0 / TWOPI

    @property
    def period_s(self) -> float:
        return TWOPI / self.mean_motion_rad_s

    @property
    def period_min(self) -> float:
        return self.period_s / 60.0

    @property
    def specific_orbital_energy_km2_s2(self) -> float:
        return -self.mu_km3_s2 / (2.0 * self.a_km)

    # Common symbolic aliases. Angles are radians.
    a = property(lambda self: self.a_km)
    h = property(lambda self: self.h_km2_s)
    p = property(lambda self: self.p_km)
    rp = property(lambda self: self.periapsis_radius_km)
    ra = property(lambda self: self.apoapsis_radius_km)
    i = property(lambda self: self.i_rad)
    raan = property(lambda self: self.raan_rad)
    argp = property(lambda self: self.argp_rad)
    nu = property(lambda self: self.nu_rad)
    theta = property(lambda self: self.nu_rad)
    E = property(lambda self: self.eccentric_anomaly_rad)
    M = property(lambda self: self.mean_anomaly_rad)
    eccentric_anomaly = property(lambda self: self.eccentric_anomaly_rad)
    mean_anomaly = property(lambda self: self.mean_anomaly_rad)

    @property
    def degrees(self) -> dict[str, float]:
        return {
            "a_km": self.a_km,
            "e": self.e,
            "i_deg": self.inclination_deg,
            "raan_deg": self.raan_deg,
            "argp_deg": self.argument_of_periapsis_deg,
            "nu_deg": self.true_anomaly_deg,
        }

    @property
    def derived(self) -> dict[str, float]:
        return {
            "p_km": self.p_km,
            "h_km2_s": self.h_km2_s,
            "periapsis_radius_km": self.periapsis_radius_km,
            "apoapsis_radius_km": self.apoapsis_radius_km,
            "radius_km": self.radius_km,
            "eccentric_anomaly_deg": self.eccentric_anomaly_deg,
            "mean_anomaly_deg": self.mean_anomaly_deg,
            "mean_motion_rad_s": self.mean_motion_rad_s,
            "mean_motion_rev_day": self.mean_motion_rev_day,
            "period_s": self.period_s,
            "specific_orbital_energy_km2_s2": self.specific_orbital_energy_km2_s2,
        }

    @property
    def vector6(self) -> np.ndarray:
        return np.array(
            [self.a_km, self.e, self.i_rad, self.raan_rad, self.argp_rad, self.nu_rad],
            dtype=float,
        )

    def to_state(self, mu_km3_s2: float | None = None) -> StateVector:
        return coe_to_rv(self, mu_km3_s2=mu_km3_s2)


def coe_to_rv(
    coes: ClassicalOrbitalElements,
    *,
    mu_km3_s2: float | None = None,
) -> StateVector:
    """Convert classical orbital elements to an ECI Cartesian state."""
    if not isinstance(coes, ClassicalOrbitalElements):
        raise TypeError("coes must be a ClassicalOrbitalElements instance.")

    mu = coes.mu_km3_s2 if mu_km3_s2 is None else float(mu_km3_s2)
    p = coes.p_km
    cnu, snu = np.cos(coes.nu_rad), np.sin(coes.nu_rad)

    r_pqw = p / (1.0 + coes.e * cnu) * np.array([cnu, snu, 0.0])
    v_pqw = np.sqrt(mu / p) * np.array([-snu, coes.e + cnu, 0.0])
    Q = perifocal_to_eci_dcm(coes.i_rad, coes.raan_rad, coes.argp_rad)
    return StateVector(Q @ r_pqw, Q @ v_pqw)


def rv_to_coe(
    r_eci_km,
    v_eci_km_s,
    *,
    mu_km3_s2: float = MU_EARTH_KM3_S2,
    tolerance: float = 1e-10,
) -> ClassicalOrbitalElements:
    """Convert an ECI Cartesian state to elliptical classical orbital elements."""
    r = _vec3(r_eci_km, "r_eci_km")
    v = _vec3(v_eci_km_s, "v_eci_km_s")
    rmag, vmag = np.linalg.norm(r), np.linalg.norm(v)

    if rmag <= 0.0:
        raise ValueError("Position magnitude must be positive.")

    h_vec = np.cross(r, v)
    hmag = np.linalg.norm(h_vec)
    if hmag <= tolerance:
        raise ValueError("State has near-zero angular momentum.")
    h_hat = h_vec / hmag

    n_vec = np.cross([0.0, 0.0, 1.0], h_vec)
    nmag = np.linalg.norm(n_vec)
    e_vec = np.cross(v, h_vec) / mu_km3_s2 - r / rmag
    e = float(np.linalg.norm(e_vec))

    energy = vmag**2 / 2.0 - mu_km3_s2 / rmag
    if abs(energy) <= tolerance:
        raise ValueError("Parabolic states are not supported by this helper.")
    a = -mu_km3_s2 / (2.0 * energy)
    if a <= 0.0 or e >= 1.0:
        raise ValueError("This helper currently supports elliptical orbits only.")

    i = float(np.arccos(np.clip(h_vec[2] / hmag, -1.0, 1.0)))
    circular, equatorial = e <= tolerance, nmag <= tolerance
    raan = wrap_2pi(np.arctan2(n_vec[1], n_vec[0])) if not equatorial else 0.0

    if not circular and not equatorial:
        argp = wrap_2pi(
            np.arctan2(
                np.dot(np.cross(n_vec, e_vec), h_hat) / nmag,
                np.dot(n_vec, e_vec) / nmag,
            )
        )
        nu = wrap_2pi(
            np.arctan2(
                np.dot(np.cross(e_vec, r), h_hat) / (e * rmag),
                np.dot(e_vec, r) / (e * rmag),
            )
        )
    elif circular and not equatorial:
        argp = 0.0
        nu = wrap_2pi(
            np.arctan2(
                np.dot(np.cross(n_vec, r), h_hat) / (nmag * rmag),
                np.dot(n_vec, r) / (nmag * rmag),
            )
        )
    elif not circular and equatorial:
        argp = wrap_2pi(np.arctan2(e_vec[1], e_vec[0]))
        nu = wrap_2pi(
            np.arctan2(
                np.dot(np.cross(e_vec, r), h_hat) / (e * rmag),
                np.dot(e_vec, r) / (e * rmag),
            )
        )
    else:
        argp, nu = 0.0, wrap_2pi(np.arctan2(r[1], r[0]))

    return ClassicalOrbitalElements(
        a_km=a,
        e=e,
        i_rad=i,
        raan_rad=raan,
        argp_rad=argp,
        nu_rad=nu,
        mu_km3_s2=mu_km3_s2,
    )


def mean_motion_rad_s(a_km: float, mu_km3_s2: float = MU_EARTH_KM3_S2) -> float:
    return float(np.sqrt(mu_km3_s2 / float(a_km) ** 3))


def orbital_period_s(a_km: float, mu_km3_s2: float = MU_EARTH_KM3_S2) -> float:
    return float(TWOPI / mean_motion_rad_s(a_km, mu_km3_s2))


__all__ = [
    "StateVector",
    "ClassicalOrbitalElements",
    "coe_to_rv",
    "rv_to_coe",
    "mean_motion_rad_s",
    "orbital_period_s",
    "wrap_2pi",
    "eccentric_anomaly_from_true",
    "true_anomaly_from_eccentric",
    "mean_anomaly_from_eccentric",
    "eccentric_anomaly_from_mean",
]
