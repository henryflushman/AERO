"""Spacecraft definitions combining an orbit and physical properties."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from functools import partial

import numpy as np

from .Atmosphere import ExponentialAtmosphere
from .OrbitalElements import ClassicalOrbitalElements, StateVector
from .Perturbations import atmospheric_drag
from .TLETools import TLERecord


@dataclass(frozen=True)
class OrbitingSpacecraft:
    """An orbit together with optional spacecraft physical properties.

    Orbital elements use the toolkit's km, km/s, and radians conventions.
    Mass is in kg, drag area is in m^2, and ``coefficient_of_drag`` is
    dimensionless. Physical properties are optional for non-drag applications.
    """

    orbital_elements: ClassicalOrbitalElements
    mass_kg: float | None = None
    ram_area_m2: float | None = None
    coefficient_of_drag: float = 2.2
    name: str | None = None
    epoch: datetime | None = None
    tle_record: TLERecord | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.orbital_elements, ClassicalOrbitalElements):
            raise TypeError(
                "orbital_elements must be a ClassicalOrbitalElements instance."
            )
        if self.mass_kg is not None and (
            not np.isfinite(self.mass_kg) or self.mass_kg <= 0.0
        ):
            raise ValueError("mass_kg must be finite and positive when provided.")
        if self.ram_area_m2 is not None and (
            not np.isfinite(self.ram_area_m2) or self.ram_area_m2 <= 0.0
        ):
            raise ValueError(
                "ram_area_m2 must be finite and positive when provided."
            )
        if not np.isfinite(self.coefficient_of_drag) or self.coefficient_of_drag <= 0.0:
            raise ValueError("coefficient_of_drag must be finite and positive.")
        if self.name is not None:
            cleaned_name = self.name.strip()
            object.__setattr__(self, "name", cleaned_name or None)

    @property
    def state(self) -> StateVector:
        """Return the Cartesian state represented by ``orbital_elements``."""
        return self.orbital_elements.to_state()

    @classmethod
    def from_tle(
        cls,
        tle: TLERecord | str,
        *,
        mass_kg: float | None = None,
        ram_area_m2: float | None = None,
        coefficient_of_drag: float = 2.2,
        name: str | None = None,
    ) -> "OrbitingSpacecraft":
        """Create a spacecraft from a TLE and optional physical properties.

        TLE mean elements are not the same as osculating classical elements.
        This method obtains the TLE's SGP4 state at its epoch and converts that
        state to osculating two-body elements for use with this toolkit's
        numerical propagators. Use ``TLERecord`` directly for SGP4 propagation.
        """
        if isinstance(tle, str):
            tle_record = TLERecord.from_text(tle)
        elif isinstance(tle, TLERecord):
            tle_record = tle
        else:
            raise TypeError("tle must be a TLERecord or TLE text string.")

        return cls(
            orbital_elements=tle_record.coes_at_epoch(),
            mass_kg=mass_kg,
            ram_area_m2=ram_area_m2,
            coefficient_of_drag=coefficient_of_drag,
            name=name or tle_record.name or f"NORAD {tle_record.norad_id}",
            epoch=tle_record.epoch,
            tle_record=tle_record,
        )

    def drag_acceleration_model(
        self,
        atmosphere: ExponentialAtmosphere,
    ):
        """Return a drag callback configured for this spacecraft.

        The result can be passed directly as an item in
        ``enckes_method(..., perturbations=[...])``.
        """
        if self.mass_kg is None:
            raise ValueError("mass_kg is required to create a drag model.")
        if self.ram_area_m2 is None:
            raise ValueError("ram_area_m2 is required to create a drag model.")
        if not isinstance(atmosphere, ExponentialAtmosphere):
            raise TypeError("atmosphere must be an ExponentialAtmosphere instance.")

        return partial(
            atmospheric_drag,
            atmosphere=atmosphere,
            mass_kg=self.mass_kg,
            area_m2=self.ram_area_m2,
            coefficient_of_drag=self.coefficient_of_drag,
        )


__all__ = ["OrbitingSpacecraft"]
