"""
Physical and mathematical constants used by IntermediateOrbitsTools.

Package conventions
-------------------
Distance : km
Time     : s
Velocity : km/s
Angles   : radians unless otherwise specified

Constants are grouped into dictionaries so new bodies or categories can be
added without changing existing code.
"""

from __future__ import annotations

import math


# ============================================================================
# Mathematical constants
# ============================================================================

MATH = {
    "pi": math.pi,
    "two_pi": 2.0 * math.pi,
    "deg_to_rad": math.pi / 180.0,
    "rad_to_deg": 180.0 / math.pi,
}


# ============================================================================
# Time constants
# ============================================================================

TIME = {
    "seconds_per_minute": 60.0,
    "seconds_per_hour": 3600.0,
    "seconds_per_day": 86400.0,
    "minutes_per_day": 1440.0,
    "julian_day_j2000": 2451545.0,
}


# ============================================================================
# Earth
# ============================================================================

EARTH = {
    "name": "Earth",

    # Gravity
    "mu_km3_s2": 398600.4418,
    "j2": 1.08262668e-3,

    # Geometry
    "radius_equatorial_km": 6378.137,
    "radius_polar_km": 6356.7523142,

    # Rotation
    "rotation_rate_rad_s": 7.2921150e-5,

    # Miscellaneous
    "mass_kg": 5.97219e24,
}


# ============================================================================
# Moon
# ============================================================================

MOON = {
    "name": "Moon",

    "mu_km3_s2": 4902.800066,
    "radius_equatorial_km": 1737.4,

    "mass_kg": 7.342e22,
}


# ============================================================================
# Sun
# ============================================================================

SUN = {
    "name": "Sun",

    "mu_km3_s2": 132712440041.9394,
    "radius_equatorial_km": 695700.0,

    "mass_kg": 1.98847e30,
}


# ============================================================================
# Mars
# ============================================================================

MARS = {
    "name": "Mars",

    "mu_km3_s2": 42828.375214,
    "j2": 1.96045e-3,

    "radius_equatorial_km": 3396.19,
    "radius_polar_km": 3376.20,

    "rotation_rate_rad_s": 7.0882181e-5,

    "mass_kg": 6.4171e23,
}


# ============================================================================
# Planet lookup
# ============================================================================

BODIES = {
    "earth": EARTH,
    "moon": MOON,
    "sun": SUN,
    "mars": MARS,
}


def get_body(name: str) -> dict:
    """
    Return the constants dictionary for a celestial body.

    Parameters
    ----------
    name : str
        Body name, such as "Earth", "Moon", "Mars", or "Sun".

    Returns
    -------
    dict
        Dictionary containing the body's constants.
    """
    key = name.strip().lower()

    if key not in BODIES:
        available = ", ".join(sorted(BODIES))
        raise KeyError(
            f"Unknown body {name!r}. "
            f"Available bodies: {available}"
        )

    return BODIES[key]


# ============================================================================
# Backward-compatible Earth aliases
# ============================================================================
#
# These names allow the existing orbital mechanics files to keep working
# without modification.
#

MU_EARTH_KM3_S2 = EARTH["mu_km3_s2"]

R_EARTH_KM = EARTH["radius_equatorial_km"]

J2_EARTH = EARTH["j2"]

OMEGA_EARTH_RAD_S = EARTH["rotation_rate_rad_s"]

SECONDS_PER_DAY = TIME["seconds_per_day"]