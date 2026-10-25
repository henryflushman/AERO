"""Perturbing accelerations for orbital propagation."""

from __future__ import annotations

import numpy as np

from .Atmosphere import ExponentialAtmosphere
from .Constants import OMEGA_EARTH_RAD_S, R_EARTH_KM


def relative_earth_velocity(
    r_km: np.ndarray,
    v_km_s: np.ndarray,
    earth_angular_velocity_rad_s: float = OMEGA_EARTH_RAD_S,
) -> np.ndarray:
    """Return velocity relative to the co-rotating atmosphere, in km/s."""
    position = np.asarray(r_km, dtype=float).reshape(-1)
    velocity = np.asarray(v_km_s, dtype=float).reshape(-1)
    if position.size != 3 or velocity.size != 3:
        raise ValueError("r_km and v_km_s must each contain three components.")
    if not np.all(np.isfinite(position)) or not np.all(np.isfinite(velocity)):
        raise ValueError("Position and velocity components must be finite.")
    if not np.isfinite(earth_angular_velocity_rad_s):
        raise ValueError("earth_angular_velocity_rad_s must be finite.")

    omega_vector = np.array([0.0, 0.0, earth_angular_velocity_rad_s])
    return velocity - np.cross(omega_vector, position)


def atmospheric_drag(
    t_s: float,
    r_km: np.ndarray,
    v_km_s: np.ndarray,
    *,
    atmosphere: ExponentialAtmosphere,
    mass_kg: float,
    area_m2: float,
    coefficient_of_drag: float = 2.2,
    earth_angular_velocity_rad_s: float = OMEGA_EARTH_RAD_S,
) -> np.ndarray:
    """Return atmospheric drag acceleration in km/s^2.

    This callable has the signature expected by ``enckes_method``. It assumes
    spherical Earth altitude and a co-rotating atmosphere. ``area_m2`` is the
    reference/projected area used by the drag model.
    """
    del t_s

    if not isinstance(atmosphere, ExponentialAtmosphere):
        raise TypeError("atmosphere must be an ExponentialAtmosphere instance.")
    if not np.isfinite(mass_kg) or mass_kg <= 0.0:
        raise ValueError("mass_kg must be finite and positive.")
    if not np.isfinite(area_m2) or area_m2 <= 0.0:
        raise ValueError("area_m2 must be finite and positive.")
    if not np.isfinite(coefficient_of_drag) or coefficient_of_drag <= 0.0:
        raise ValueError("coefficient_of_drag must be finite and positive.")

    position = np.asarray(r_km, dtype=float).reshape(-1)
    if position.size != 3 or not np.all(np.isfinite(position)):
        raise ValueError("r_km must contain three finite components.")
    altitude_km = float(np.linalg.norm(position) - R_EARTH_KM)
    density_kg_m3 = atmosphere.density_kg_m3(altitude_km)

    relative_velocity_km_s = relative_earth_velocity(
        position,
        v_km_s,
        earth_angular_velocity_rad_s,
    )
    relative_velocity_m_s = 1000.0 * relative_velocity_km_s
    acceleration_m_s2 = (
        -0.5
        * coefficient_of_drag
        * area_m2
        / mass_kg
        * density_kg_m3
        * np.linalg.norm(relative_velocity_m_s)
        * relative_velocity_m_s
    )
    return acceleration_m_s2 / 1000.0


__all__ = ["atmospheric_drag", "relative_earth_velocity"]
