"""Direction-cosine matrices and common orbital coordinate transformations.

Conventions
-----------
All DCMs in this file are used as::

    vector_in_new_frame = C_new_from_old @ vector_in_old_frame

The LVLH frame used here is the RTN convention:
- +R / +x: radial, away from Earth
- +T / +y: transverse / along-track
- +N / +z: orbit normal, along angular momentum

The orbital-plane frame is the classical perifocal PQW frame.
"""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np


def _vec3(value, name: str) -> np.ndarray:
    vector = np.asarray(value, dtype=float).reshape(-1)
    if vector.size != 3:
        raise ValueError(f"{name} must contain exactly three values.")
    return vector


def unit(vector, *, name: str = "vector") -> np.ndarray:
    """Return a unit vector and reject zero-length inputs."""
    vector = _vec3(vector, name)
    norm = np.linalg.norm(vector)
    if norm <= 0.0:
        raise ValueError(f"{name} must have nonzero magnitude.")
    return vector / norm


def rot1(angle_rad: float) -> np.ndarray:
    """Active right-handed rotation about +x."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, c, -s],
            [0.0, s, c],
        ]
    )


def rot2(angle_rad: float) -> np.ndarray:
    """Active right-handed rotation about +y."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array(
        [
            [c, 0.0, s],
            [0.0, 1.0, 0.0],
            [-s, 0.0, c],
        ]
    )


def rot3(angle_rad: float) -> np.ndarray:
    """Active right-handed rotation about +z."""
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    return np.array(
        [
            [c, -s, 0.0],
            [s, c, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )


def dcm_from_axes(x_axis, y_axis, z_axis) -> np.ndarray:
    """Build a DCM whose rows are the new frame axes written in the old frame.

    The supplied axes must already describe one orthogonal frame. Small floating
    point errors are tolerated, but the function rejects clearly non-orthogonal
    or left-handed inputs.
    """
    x_hat = unit(x_axis, name="x_axis")
    y_hat = unit(y_axis, name="y_axis")
    z_hat = unit(z_axis, name="z_axis")

    C = np.vstack((x_hat, y_hat, z_hat))
    if not np.allclose(C @ C.T, np.eye(3), atol=1e-10):
        raise ValueError("Axes must be mutually orthogonal.")
    if np.linalg.det(C) < 0.0:
        raise ValueError("Axes must form a right-handed coordinate system.")
    return C


def eci_to_lvlh_dcm(r_eci_km, v_eci_km_s) -> np.ndarray:
    """Return the ECI -> LVLH/RTN DCM for the supplied spacecraft state."""
    r = _vec3(r_eci_km, "r_eci_km")
    v = _vec3(v_eci_km_s, "v_eci_km_s")

    r_hat = unit(r, name="r_eci_km")
    h_hat = unit(np.cross(r, v), name="angular momentum")
    t_hat = unit(np.cross(h_hat, r_hat), name="along-track axis")
    return dcm_from_axes(r_hat, t_hat, h_hat)


def lvlh_to_eci_dcm(r_eci_km, v_eci_km_s) -> np.ndarray:
    """Return the LVLH/RTN -> ECI DCM."""
    return eci_to_lvlh_dcm(r_eci_km, v_eci_km_s).T


def eci_to_lvlh(vector_eci, r_eci_km, v_eci_km_s) -> np.ndarray:
    """Rotate a vector from ECI components into LVLH/RTN components."""
    return eci_to_lvlh_dcm(r_eci_km, v_eci_km_s) @ _vec3(vector_eci, "vector_eci")


def lvlh_to_eci(vector_lvlh, r_eci_km, v_eci_km_s) -> np.ndarray:
    """Rotate a vector from LVLH/RTN components into ECI components."""
    return lvlh_to_eci_dcm(r_eci_km, v_eci_km_s) @ _vec3(vector_lvlh, "vector_lvlh")


def relative_state_eci_to_lvlh(
    chief_r_eci_km,
    chief_v_eci_km_s,
    deputy_r_eci_km,
    deputy_v_eci_km_s,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert an ECI relative state into the rotating chief LVLH frame.

    Unlike rotating a velocity vector with the DCM alone, this includes the
    ``omega x rho`` term caused by rotation of the LVLH frame.
    """
    rc = _vec3(chief_r_eci_km, "chief_r_eci_km")
    vc = _vec3(chief_v_eci_km_s, "chief_v_eci_km_s")
    rd = _vec3(deputy_r_eci_km, "deputy_r_eci_km")
    vd = _vec3(deputy_v_eci_km_s, "deputy_v_eci_km_s")

    C = eci_to_lvlh_dcm(rc, vc)
    rho_eci = rd - rc
    rho_dot_eci = vd - vc

    rho_lvlh = C @ rho_eci
    h_mag = np.linalg.norm(np.cross(rc, vc))
    omega_lvlh = np.array([0.0, 0.0, h_mag / np.dot(rc, rc)])
    rho_dot_lvlh = C @ rho_dot_eci - np.cross(omega_lvlh, rho_lvlh)
    return rho_lvlh, rho_dot_lvlh


def relative_state_lvlh_to_eci(
    chief_r_eci_km,
    chief_v_eci_km_s,
    rho_lvlh,
    rho_dot_lvlh,
) -> tuple[np.ndarray, np.ndarray]:
    """Convert a rotating LVLH relative state into a deputy ECI state."""
    rc = _vec3(chief_r_eci_km, "chief_r_eci_km")
    vc = _vec3(chief_v_eci_km_s, "chief_v_eci_km_s")
    rho = _vec3(rho_lvlh, "rho_lvlh")
    rho_dot = _vec3(rho_dot_lvlh, "rho_dot_lvlh")

    C_lvlh_to_eci = lvlh_to_eci_dcm(rc, vc)

    h_mag = np.linalg.norm(np.cross(rc, vc))
    r_mag = np.linalg.norm(rc)

    omega_lvlh = np.array([
        0.0,
        0.0,
        h_mag / r_mag**2,
    ])

    rho_eci = C_lvlh_to_eci @ rho

    rho_dot_eci = C_lvlh_to_eci @ (
        rho_dot + np.cross(omega_lvlh, rho)
    )

    deputy_r_eci = rc + rho_eci
    deputy_v_eci = vc + rho_dot_eci

    return deputy_r_eci, deputy_v_eci


def perifocal_to_eci_dcm(i_rad: float, raan_rad: float, argp_rad: float) -> np.ndarray:
    """Return the PQW/perifocal -> ECI DCM."""
    return rot3(raan_rad) @ rot1(i_rad) @ rot3(argp_rad)


def eci_to_perifocal_dcm(i_rad: float, raan_rad: float, argp_rad: float) -> np.ndarray:
    """Return the ECI -> PQW/perifocal DCM."""
    return perifocal_to_eci_dcm(i_rad, raan_rad, argp_rad).T


def eci_to_orbital_plane(vector_eci, i_rad: float, raan_rad: float, argp_rad: float) -> np.ndarray:
    """Rotate an ECI vector into the perifocal/orbital-plane PQW frame."""
    return eci_to_perifocal_dcm(i_rad, raan_rad, argp_rad) @ _vec3(vector_eci, "vector_eci")


def orbital_plane_to_eci(vector_pqw, i_rad: float, raan_rad: float, argp_rad: float) -> np.ndarray:
    """Rotate a perifocal/PQW vector into ECI."""
    return perifocal_to_eci_dcm(i_rad, raan_rad, argp_rad) @ _vec3(vector_pqw, "vector_pqw")


def _as_utc(value: datetime) -> datetime:
    if not isinstance(value, datetime):
        raise TypeError("time must be a datetime.")
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def julian_date(time: datetime) -> float:
    """Convert a UTC datetime to Julian Date."""
    t = _as_utc(time)
    year = t.year
    month = t.month
    day_fraction = (
        t.day
        + (t.hour + (t.minute + (t.second + t.microsecond / 1e6) / 60.0) / 60.0) / 24.0
    )

    if month <= 2:
        year -= 1
        month += 12

    A = year // 100
    B = 2 - A + A // 4
    return (
        int(365.25 * (year + 4716))
        + int(30.6001 * (month + 1))
        + day_fraction
        + B
        - 1524.5
    )


def gmst_rad(time: datetime) -> float:
    """Greenwich mean sidereal time in radians.

    This is appropriate for classroom-level ECI/ECEF visualization work. It does
    not include polar motion or the full high-precision Earth orientation model.
    """
    jd = julian_date(time)
    T = (jd - 2451545.0) / 36525.0
    gmst_deg = (
        280.46061837
        + 360.98564736629 * (jd - 2451545.0)
        + 0.000387933 * T**2
        - T**3 / 38710000.0
    ) % 360.0
    return np.deg2rad(gmst_deg)


def eci_to_ecef(vector_eci, time: datetime) -> np.ndarray:
    """Rotate a position/vector from a simple ECI frame into Earth-fixed ECEF."""
    return rot3(-gmst_rad(time)) @ _vec3(vector_eci, "vector_eci")


def ecef_to_eci(vector_ecef, time: datetime) -> np.ndarray:
    """Rotate a position/vector from Earth-fixed ECEF into a simple ECI frame."""
    return rot3(gmst_rad(time)) @ _vec3(vector_ecef, "vector_ecef")


def teme_to_ecef(r_teme_km, time: datetime) -> np.ndarray:
    """Approximate TEME -> ECEF position conversion using GMST rotation.

    This is convenient for TLE/Cesium visualization. For precision orbit
    determination, use a full Earth-orientation transformation library.
    """
    return eci_to_ecef(r_teme_km, time)


__all__ = [
    "unit",
    "rot1",
    "rot2",
    "rot3",
    "dcm_from_axes",
    "eci_to_lvlh_dcm",
    "lvlh_to_eci_dcm",
    "eci_to_lvlh",
    "lvlh_to_eci",
    "relative_state_eci_to_lvlh",
    "relative_state_lvlh_to_eci",
    "perifocal_to_eci_dcm",
    "eci_to_perifocal_dcm",
    "eci_to_orbital_plane",
    "orbital_plane_to_eci",
    "julian_date",
    "gmst_rad",
    "eci_to_ecef",
    "ecef_to_eci",
    "teme_to_ecef",
]
