"""Atmospheric density models used by IntermediateOrbitsTools."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


# Standard sea-level air properties: rho = p / (R_specific * T).
SEA_LEVEL_PRESSURE_PA = 101325.0
SEA_LEVEL_TEMPERATURE_K = 288.15
AIR_SPECIFIC_GAS_CONSTANT_J_KG_K = 287.05287
SEA_LEVEL_DENSITY_KG_M3 = (
    SEA_LEVEL_PRESSURE_PA
    / (AIR_SPECIFIC_GAS_CONSTANT_J_KG_K * SEA_LEVEL_TEMPERATURE_K)
)


@dataclass(frozen=True, init=False)
class ExponentialAtmosphere:
    """Single-layer exponential density model.

    Density follows ``rho = rho_ref * exp(-(h - h_ref) / H)``.
    Altitude and scale height are in km; density is in kg/m^3. By default,
    the reference is standard sea-level air at 101325 Pa and 288.15 K.
    Scale height must be provided for the selected atmosphere model.
    """

    reference_density_kg_m3: float
    reference_altitude_km: float
    scale_height_km: float

    def __init__(
        self,
        reference_density_kg_m3: float = SEA_LEVEL_DENSITY_KG_M3,
        reference_altitude_km: float = 0.0,
        scale_height_km: float | None = None,
    ) -> None:
        if scale_height_km is None:
            raise ValueError(
                "scale_height_km is required for the selected atmosphere model."
            )
        object.__setattr__(self, "reference_density_kg_m3", reference_density_kg_m3)
        object.__setattr__(self, "reference_altitude_km", reference_altitude_km)
        object.__setattr__(self, "scale_height_km", scale_height_km)
        self.__post_init__()

    def __post_init__(self) -> None:
        parameters = (
            self.reference_density_kg_m3,
            self.reference_altitude_km,
            self.scale_height_km,
        )
        if not all(np.isfinite(value) for value in parameters):
            raise ValueError("Atmosphere parameters must be finite.")
        if self.reference_density_kg_m3 <= 0.0:
            raise ValueError("reference_density_kg_m3 must be positive.")
        if self.scale_height_km <= 0.0:
            raise ValueError("scale_height_km must be positive.")

    def density_kg_m3(self, altitude_km: float | np.ndarray) -> float | np.ndarray:
        """Return density at one or more altitudes above Earth's surface."""
        altitude = np.asarray(altitude_km, dtype=float)
        if not np.all(np.isfinite(altitude)):
            raise ValueError("altitude_km values must be finite.")

        density = self.reference_density_kg_m3 * np.exp(
            -(altitude - self.reference_altitude_km) / self.scale_height_km
        )
        if altitude.ndim == 0:
            return float(density)
        return density


__all__ = [
    "AIR_SPECIFIC_GAS_CONSTANT_J_KG_K",
    "ExponentialAtmosphere",
    "SEA_LEVEL_DENSITY_KG_M3",
    "SEA_LEVEL_PRESSURE_PA",
    "SEA_LEVEL_TEMPERATURE_K",
]
