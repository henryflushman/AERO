"""
mli_blanket.py
==============
A Python class for defining, building, and analysing Multi-Layer Insulation (MLI)
blanket configurations used in spacecraft thermal control.

Construction routes
-------------------
1. MLIBlanket(layers)                          – from an explicit list of Layer objects
2. MLIBlanket.from_json(path_or_dict)          – from a JSON file path or dict
3. MLIBlanket.from_effective_emissivity(...)   – back-calculate layer count from ε_eff
4. MLIBlanket.from_layer_count(n, material)    – uniform blanket of n identical layers
5. MLIBlanket.from_layer_types(layer_spec)     – mixed blanket from {material: count} dict
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Union


# ---------------------------------------------------------------------------
# Material library
# ---------------------------------------------------------------------------

#: Built-in material library keyed by canonical name.
#: Each entry stores the *single-surface* emissivity (ε) and a description.
MATERIAL_LIBRARY: Dict[str, Dict] = {
    "double_aluminized_mylar": {
        "emissivity": 0.03,
        "description": "Double-aluminized Mylar (DAM), 6 µm",
        "alias": ["dam", "mylar"],
    },
    "single_aluminized_mylar": {
        "emissivity": 0.04,
        "description": "Single-aluminized Mylar, 6 µm",
        "alias": ["sam"],
    },
    "aluminized_kapton": {
        "emissivity": 0.035,
        "description": "Aluminized Kapton, 25 µm",
        "alias": ["kapton"],
    },
    "aluminized_pet": {
        "emissivity": 0.04,
        "description": "Aluminized PET (polyethylene terephthalate)",
        "alias": ["pet"],
    },
    "silverized_teflon": {
        "emissivity": 0.025,
        "description": "Silverized Teflon (FEP)",
        "alias": ["fep", "teflon"],
    },
    "dacron_net": {
        "emissivity": 0.85,
        "description": "Dacron net spacer",
        "alias": ["dacron"],
    },
    "silk_net": {
        "emissivity": 0.78,
        "description": "Silk net spacer",
        "alias": ["silk"],
    },
    "custom": {
        "emissivity": None,
        "description": "User-defined custom material",
        "alias": [],
    },
}

STEFAN_BOLTZMANN = 5.670374419e-8  # W m⁻² K⁻⁴


def _resolve_material(name: str) -> str:
    """Return the canonical material key for *name* (case-insensitive alias lookup)."""
    key = name.lower().strip()
    if key in MATERIAL_LIBRARY:
        return key
    for canonical, info in MATERIAL_LIBRARY.items():
        if key in [a.lower() for a in info.get("alias", [])]:
            return canonical
    raise ValueError(
        f"Unknown material '{name}'. "
        f"Available materials: {list(MATERIAL_LIBRARY.keys())}"
    )


# ---------------------------------------------------------------------------
# Layer dataclass
# ---------------------------------------------------------------------------

@dataclass
class Layer:
    """Represents a single reflective or spacer layer in an MLI blanket.

    Parameters
    ----------
    material : str
        Canonical material name or alias (resolved against MATERIAL_LIBRARY).
    emissivity : float, optional
        Override the library emissivity (required when material='custom').
    thickness_um : float, optional
        Physical thickness in micrometres (informational only).
    is_spacer : bool
        True if this layer is a spacer (not a reflective shield).
    name : str, optional
        Human-readable label for the layer.
    """

    material: str
    emissivity: Optional[float] = None
    thickness_um: Optional[float] = None
    is_spacer: bool = False
    name: str = ""

    def __post_init__(self) -> None:
        self.material = _resolve_material(self.material)
        lib_emissivity = MATERIAL_LIBRARY[self.material]["emissivity"]

        if self.emissivity is None:
            if lib_emissivity is None:
                raise ValueError(
                    "Emissivity must be provided explicitly for 'custom' material."
                )
            self.emissivity = lib_emissivity
        else:
            if not (0.0 < self.emissivity <= 1.0):
                raise ValueError(
                    f"Emissivity must be in (0, 1], got {self.emissivity}."
                )

        if not self.name:
            self.name = MATERIAL_LIBRARY[self.material]["description"]

    # ------------------------------------------------------------------
    # Convenience constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_dict(cls, data: dict) -> "Layer":
        """Construct a Layer from a plain dictionary (e.g. parsed from JSON)."""
        return cls(
            material=data["material"],
            emissivity=data.get("emissivity"),
            thickness_um=data.get("thickness_um"),
            is_spacer=data.get("is_spacer", False),
            name=data.get("name", ""),
        )

    def to_dict(self) -> dict:
        """Serialise the layer to a plain dictionary."""
        return {
            "material": self.material,
            "emissivity": self.emissivity,
            "thickness_um": self.thickness_um,
            "is_spacer": self.is_spacer,
            "name": self.name,
        }

    def __repr__(self) -> str:
        return (
            f"Layer(material={self.material!r}, ε={self.emissivity:.4f}, "
            f"spacer={self.is_spacer}, name={self.name!r})"
        )


# ---------------------------------------------------------------------------
# MLIBlanket class
# ---------------------------------------------------------------------------

class MLIBlanket:
    """Multi-Layer Insulation (MLI) blanket model.

    The blanket is stored internally as an ordered list of :class:`Layer` objects.
    All thermal calculations assume pure radiation between layers (no conduction)
    and the standard parallel-plate radiative resistance model.

    The **first** and **last** layers in the stack are the two boundary surfaces
    (ε₁ and ε₂). Every layer in between is an interior shield counted by n in
    the formula. For a stack of N total layers, n = N − 2.

    Effective emissivity formula:

        ε* = ( 2n/ε_mylar − n − 1 + 1/ε₁ + 1/ε₂ )⁻¹

    Which in general series-resistance form (supports mixed materials) is:

        1/ε* = 1/ε₁ + 1/ε₂ − 1 + Σᵢ (2/εᵢ − 1)

    where ε₁ = layers[0].emissivity, ε₂ = layers[-1].emissivity, and the sum
    runs over every interior *reflective* (non-spacer) layer.

    Parameters
    ----------
    layers : list of Layer
        Full ordered stack including boundary surfaces as first and last entries.
        Must contain at least 2 layers (the two boundary surfaces).
    name : str, optional
        Blanket identifier / description.
    notes : str, optional
        Free-text notes stored alongside the blanket.
    """

    def __init__(
        self,
        layers: List[Layer],
        name: str = "MLI Blanket",
        notes: str = "",
    ) -> None:
        if len(layers) < 2:
            raise ValueError(
                "layers must contain at least 2 entries "
                "(the outer and inner boundary surfaces)."
            )
        self._layers: List[Layer] = list(layers)
        self.name = name
        self.notes = notes

    # ==================================================================
    # Alternative constructors
    # ==================================================================

    @classmethod
    def from_json(cls, source: Union[str, Path, dict]) -> "MLIBlanket":
        """Build an MLIBlanket from a JSON file path, Path object, or dict.

        Expected JSON schema::

            {
                "name": "My Blanket",          // optional
                "notes": "...",                // optional
                "layers": [
                    {
                        "material": "double_aluminized_mylar",
                        "emissivity": 0.03,    // optional override
                        "thickness_um": 6,     // optional
                        "is_spacer": false,
                        "name": "Outer shield" // optional
                    },
                    ...
                ]
            }

        The first and last layer entries are treated as the boundary surfaces
        (ε₁ and ε₂). All layers in between are interior shields.
        """
        if isinstance(source, dict):
            data = source
        else:
            path = Path(source)
            with path.open("r", encoding="utf-8") as fh:
                data = json.load(fh)

        layers = [Layer.from_dict(ld) for ld in data["layers"]]
        return cls(
            layers=layers,
            name=data.get("name", "MLI Blanket"),
            notes=data.get("notes", ""),
        )

    @classmethod
    def from_effective_emissivity(
        cls,
        effective_emissivity: float,
        layer_material: str = "double_aluminized_mylar",
        layer_emissivity: Optional[float] = None,
        name: str = "MLI Blanket",
    ) -> "MLIBlanket":
        """Back-calculate the total number of layers required to achieve
        *effective_emissivity* for a uniform single-material blanket.

        The first and last layers are the boundary surfaces (ε₁ = ε₂ = ε_mylar).
        The number of interior shields n is solved from:

            n = (1/ε* − 2/ε + 1) / (2/ε − 1)

        Total layers returned = n + 2.

        Parameters
        ----------
        effective_emissivity : float
            Target effective emissivity (0 < ε* < 1).
        layer_material : str
            Material for all layers (boundary surfaces and interior shields).
        layer_emissivity : float, optional
            Override the library emissivity for the chosen material.
        name : str
            Blanket label.
        """
        if not (0.0 < effective_emissivity < 1.0):
            raise ValueError("effective_emissivity must be in (0, 1).")

        canonical = _resolve_material(layer_material)
        eps = layer_emissivity or MATERIAL_LIBRARY[canonical]["emissivity"]
        if eps is None:
            raise ValueError("Provide layer_emissivity for custom materials.")

        # With ε₁ = ε₂ = ε: 1/ε* = 2/ε - 1 + n*(2/ε - 1) = (n+1)*(2/ε - 1)
        # So n = 1/ε* / (2/ε - 1) - 1
        two_over_eps_minus_one = 2.0 / eps - 1.0
        n_float = (1.0 / effective_emissivity) / two_over_eps_minus_one - 1.0
        n_interior = max(0, math.ceil(n_float))
        total = n_interior + 2  # add the two boundary surface layers

        layers = [
            Layer(
                material=layer_material,
                emissivity=layer_emissivity,
                name=("Outer boundary surface" if i == 0
                      else "Inner boundary surface" if i == total - 1
                      else f"Interior shield {i}"),
            )
            for i in range(total)
        ]
        return cls(
            layers=layers,
            name=name,
            notes=(
                f"Auto-generated from target ε*={effective_emissivity:.4f}; "
                f"n_interior={n_float:.2f} → {n_interior}, total layers={total}."
            ),
        )

    @classmethod
    def from_layer_count(
        cls,
        n_layers: int,
        material: str = "double_aluminized_mylar",
        layer_emissivity: Optional[float] = None,
        name: str = "MLI Blanket",
    ) -> "MLIBlanket":
        """Build a uniform blanket of *n_layers* identical layers.

        The first and last layers are the boundary surfaces (ε₁ and ε₂).
        The remaining n_layers − 2 are interior shields. Minimum 2 layers.

        Parameters
        ----------
        n_layers : int
            Total number of layers including the two boundary surfaces (≥ 2).
        material : str
            Material applied to every layer.
        layer_emissivity : float, optional
            Override the library emissivity.
        name : str
            Blanket label.
        """
        if n_layers < 2:
            raise ValueError("n_layers must be ≥ 2 (at least the two boundary surfaces).")

        layers = [
            Layer(
                material=material,
                emissivity=layer_emissivity,
                name=("Outer boundary surface" if i == 0
                      else "Inner boundary surface" if i == n_layers - 1
                      else f"Interior shield {i}"),
            )
            for i in range(n_layers)
        ]
        return cls(layers=layers, name=name)

    @classmethod
    def from_layer_types(
        cls,
        layer_spec: Dict[str, Union[int, dict]],
        name: str = "MLI Blanket",
        notes: str = "",
    ) -> "MLIBlanket":
        """Build a blanket from a dict mapping material names to layer counts
        (or detailed configuration dicts).

        Simple form — equal distribution of each material::

            MLIBlanket.from_layer_types({
                "double_aluminized_mylar": 15,
                "aluminized_kapton": 5,
            })

        Extended form — per-material overrides::

            MLIBlanket.from_layer_types({
                "double_aluminized_mylar": {
                    "count": 10,
                    "emissivity": 0.03,
                    "thickness_um": 6,
                },
                "aluminized_kapton": {
                    "count": 5,
                    "emissivity": 0.035,
                },
            })

        Layers are stacked in the dict iteration order (Python 3.7+ preserves
        insertion order). To interleave materials, pass an explicit layer list to
        the normal constructor instead.

        Parameters
        ----------
        layer_spec : dict
            Mapping of ``{material: count}`` or ``{material: {count, ...}}``.
            The first layer of the first material and last layer of the last
            material become the boundary surfaces (ε₁ and ε₂).
        name : str
            Blanket label.
        notes : str
            Free-text notes.
        """
        layers: List[Layer] = []
        layer_index = 1
        for material, spec in layer_spec.items():
            if isinstance(spec, int):
                count = spec
                emissivity = None
                thickness_um = None
                is_spacer = False
            elif isinstance(spec, dict):
                count = int(spec.get("count", 1))
                emissivity = spec.get("emissivity")
                thickness_um = spec.get("thickness_um")
                is_spacer = bool(spec.get("is_spacer", False))
            else:
                raise TypeError(
                    f"layer_spec values must be int or dict, got {type(spec)} "
                    f"for material '{material}'."
                )

            for i in range(count):
                layers.append(
                    Layer(
                        material=material,
                        emissivity=emissivity,
                        thickness_um=thickness_um,
                        is_spacer=is_spacer,
                        name=f"Shield {layer_index}",
                    )
                )
                layer_index += 1

        return cls(layers=layers, name=name, notes=notes)

    @classmethod
    def from_emissivity_list(
        cls,
        emissivities: List[float],
        material: str = "custom",
        name: str = "MLI Blanket",
        notes: str = "",
    ) -> "MLIBlanket":
        """Build a blanket by specifying each layer's emissivity individually.

        The first and last values in *emissivities* are treated as the boundary
        surfaces (ε₁ and ε₂). All values in between are interior shields.

        Parameters
        ----------
        emissivities : list of float
            Ordered emissivity values for every layer including boundary surfaces.
            Must have at least 2 entries.
        material : str
            Material label applied to every layer.
        name : str
            Blanket label.
        notes : str
            Free-text notes.
        """
        if len(emissivities) < 2:
            raise ValueError("emissivities must have at least 2 values.")
        layers = [
            Layer(
                material=material,
                emissivity=eps,
                name=("Outer boundary surface" if i == 0
                      else "Inner boundary surface" if i == len(emissivities) - 1
                      else f"Interior shield {i}"),
            )
            for i, eps in enumerate(emissivities)
        ]
        return cls(layers=layers, name=name, notes=notes)

    # ==================================================================
    # Core properties
    # ==================================================================

    @property
    def layers(self) -> List[Layer]:
        """Full ordered layer list."""
        return list(self._layers)

    @property
    def reflective_layers(self) -> List[Layer]:
        """Layers that are *not* spacers."""
        return [l for l in self._layers if not l.is_spacer]

    @property
    def spacer_layers(self) -> List[Layer]:
        """Layers that are spacers."""
        return [l for l in self._layers if l.is_spacer]

    @property
    def n_layers(self) -> int:
        """Total number of layers (reflective + spacer)."""
        return len(self._layers)

    @property
    def n_reflective_layers(self) -> int:
        """Number of reflective (non-spacer) layers."""
        return len(self.reflective_layers)

    @property
    def n_spacer_layers(self) -> int:
        """Number of spacer layers."""
        return len(self.spacer_layers)

    @property
    def outer_surface_emissivity(self) -> float:
        """Emissivity of the outer boundary surface (ε₁ = layers[0].emissivity)."""
        return self._layers[0].emissivity

    @property
    def inner_surface_emissivity(self) -> float:
        """Emissivity of the inner boundary surface (ε₂ = layers[-1].emissivity)."""
        return self._layers[-1].emissivity

    @property
    def interior_reflective_layers(self) -> List[Layer]:
        """Reflective layers that are NOT the boundary surfaces (layers[1:-1])."""
        return [l for l in self._layers[1:-1] if not l.is_spacer]

    @property
    def n_interior_layers(self) -> int:
        """Number of interior reflective shields (n = n_layers − 2)."""
        return len(self.interior_reflective_layers)

    @property
    def effective_emissivity(self) -> float:
        """Effective emissivity using the correct boundary-surface formula.

        The first and last layers are ε₁ and ε₂. Interior non-spacer layers
        (layers[1:-1]) are the n shields summed over:

            1/ε* = 1/ε₁ + 1/ε₂ − 1 + Σᵢ (2/εᵢ − 1)
        """
        eps1 = self.outer_surface_emissivity
        eps2 = self.inner_surface_emissivity
        interior = self.interior_reflective_layers
        inv_eps_star = (
            1.0 / eps1
            + 1.0 / eps2
            - 1.0
            + sum(2.0 / layer.emissivity - 1.0 for layer in interior)
        )
        return 1.0 / inv_eps_star

    @property
    def layer_emissivities(self) -> List[float]:
        """Emissivity of every layer (in stack order)."""
        return [layer.emissivity for layer in self._layers]

    @property
    def mean_layer_emissivity(self) -> float:
        """Arithmetic mean emissivity of reflective layers."""
        eps_list = [l.emissivity for l in self.reflective_layers]
        return sum(eps_list) / len(eps_list)

    @property
    def layer_materials(self) -> List[str]:
        """Material name for every layer (in stack order)."""
        return [layer.material for layer in self._layers]

    @property
    def material_counts(self) -> Dict[str, int]:
        """Number of layers for each unique material present."""
        counts: Dict[str, int] = {}
        for layer in self._layers:
            counts[layer.material] = counts.get(layer.material, 0) + 1
        return counts

    # ==================================================================
    # Thermal calculations
    # ==================================================================

    def heat_flux(self, T_hot: float, T_cold: float) -> float:
        """Radiative heat flux through the blanket (W m⁻²).

        Parameters
        ----------
        T_hot : float
            Hot-side temperature in Kelvin.
        T_cold : float
            Cold-side temperature in Kelvin.

        Returns
        -------
        float
            Net heat flux q = ε_eff · σ · (T_hot⁴ − T_cold⁴)  [W m⁻²]
        """
        if T_hot < 0 or T_cold < 0:
            raise ValueError("Temperatures must be non-negative (Kelvin).")
        return (
            self.effective_emissivity
            * STEFAN_BOLTZMANN
            * (T_hot**4 - T_cold**4)
        )

    def heat_transfer_rate(
        self, T_hot: float, T_cold: float, area: float
    ) -> float:
        """Total heat transfer rate through the blanket (W).

        Parameters
        ----------
        T_hot : float
            Hot-side temperature [K].
        T_cold : float
            Cold-side temperature [K].
        area : float
            Blanket surface area [m²].

        Returns
        -------
        float
            Q = q · A  [W]
        """
        return self.heat_flux(T_hot, T_cold) * area

    def thermal_resistance(self, T_hot: float, T_cold: float) -> float:
        """Apparent radiative thermal resistance (K W⁻¹ m²).

        R = ΔT / q

        Returns ``inf`` when heat flux is zero (i.e., T_hot == T_cold).
        """
        q = self.heat_flux(T_hot, T_cold)
        if q == 0.0:
            return float("inf")
        return abs(T_hot - T_cold) / abs(q)

    def equivalent_n_layers(
        self,
        reference_material: str = "double_aluminized_mylar",
        reference_emissivity: Optional[float] = None,
    ) -> float:
        """Express this (possibly mixed) blanket as an equivalent number of
        uniform interior layers of *reference_material* with the same ε*.

        Rearranged from the full formula:

            n_eq = (1/ε* − 1/ε₁ − 1/ε₂ + 1) / (2/ε_ref − 1)
        """
        canonical = _resolve_material(reference_material)
        eps_ref = reference_emissivity or MATERIAL_LIBRARY[canonical]["emissivity"]
        if eps_ref is None:
            raise ValueError("Provide reference_emissivity for custom materials.")
        return (
            1.0 / self.effective_emissivity
            - 1.0 / self.outer_surface_emissivity
            - 1.0 / self.inner_surface_emissivity
            + 1.0
        ) / (2.0 / eps_ref - 1.0)

    # ==================================================================
    # Layer manipulation
    # ==================================================================

    def add_layer(self, layer: Layer, position: Optional[int] = None) -> None:
        """Insert *layer* at *position* (default: append to end).

        Parameters
        ----------
        layer : Layer
        position : int, optional
            0-indexed position in the stack. Negative indices are supported.
        """
        if position is None:
            self._layers.append(layer)
        else:
            self._layers.insert(position, layer)

    def remove_layer(self, index: int) -> Layer:
        """Remove and return the layer at *index*."""
        return self._layers.pop(index)

    def replace_layer(self, index: int, layer: Layer) -> Layer:
        """Replace the layer at *index* with *layer*, returning the old layer."""
        old = self._layers[index]
        self._layers[index] = layer
        return old

    # ==================================================================
    # Serialisation
    # ==================================================================

    def to_dict(self) -> dict:
        """Serialise the blanket to a plain dictionary."""
        return {
            "name": self.name,
            "notes": self.notes,
            "layers": [layer.to_dict() for layer in self._layers],
        }

    def to_json(
        self,
        path: Optional[Union[str, Path]] = None,
        indent: int = 2,
    ) -> str:
        """Serialise the blanket to a JSON string, optionally writing to *path*.

        Parameters
        ----------
        path : str or Path, optional
            If given, the JSON is written to this file.
        indent : int
            JSON indentation level.

        Returns
        -------
        str
            The serialised JSON string.
        """
        json_str = json.dumps(self.to_dict(), indent=indent)
        if path is not None:
            Path(path).write_text(json_str, encoding="utf-8")
        return json_str

    # ==================================================================
    # Display helpers
    # ==================================================================

    def summary(self) -> str:
        """Return a formatted multi-line summary string."""
        lines = [
            "=" * 60,
            f"  MLI Blanket: {self.name}",
            "=" * 60,
            f"  Total layers (N)   : {self.n_layers}",
            f"  Interior shields   : {self.n_interior_layers}  (n = N − 2)",
            f"  Spacer layers      : {self.n_spacer_layers}",
            f"  Outer surface ε₁   : {self.outer_surface_emissivity:.4f}  [{self._layers[0].name}]",
            f"  Inner surface ε₂   : {self.inner_surface_emissivity:.4f}  [{self._layers[-1].name}]",
            f"  Effective ε*       : {self.effective_emissivity:.6f}",
            f"  Mean interior ε    : {self.mean_layer_emissivity:.6f}",
        ]
        if self.notes:
            lines.append(f"  Notes              : {self.notes}")
        lines.append("")
        lines.append("  Material composition:")
        for mat, count in self.material_counts.items():
            desc = MATERIAL_LIBRARY[mat]["description"]
            lines.append(f"    {count:3d} × {desc}")
        lines.append("")
        lines.append("  Layer stack (outer → inner):")
        for i, layer in enumerate(self._layers):
            tag = " [spacer] " if layer.is_spacer else "          "
            boundary = " ← ε₁" if i == 0 else (" ← ε₂" if i == self.n_layers - 1 else "")
            lines.append(
                f"    [{i:02d}]{tag}ε={layer.emissivity:.4f}  {layer.name}{boundary}"
            )
        lines.append("=" * 60)
        return "\n".join(lines)

    def __repr__(self) -> str:
        return (
            f"MLIBlanket(name={self.name!r}, n_layers={self.n_layers}, "
            f"n_interior={self.n_interior_layers}, ε*={self.effective_emissivity:.5f})"
        )

    def __len__(self) -> int:
        return self.n_layers


# ---------------------------------------------------------------------------
# Quick usage examples (run as script)
# ---------------------------------------------------------------------------

if __name__ == "__main__":

    print("\n── Example 1: 25-layer uniform DAM blanket ──────────────────")
    # 25 total: layers[0] = ε₁, layers[24] = ε₂, 23 interior shields
    b1 = MLIBlanket.from_layer_count(25, material="double_aluminized_mylar")
    print(b1.summary())
    print(f"  n_interior (should be 23): {b1.n_interior_layers}")
    print(f"  Heat flux (300 K → 100 K): {b1.heat_flux(300, 100):.4f} W/m²")

    print("\n── Example 2: back-calculate total layers from ε*=0.005 ─────")
    b2 = MLIBlanket.from_effective_emissivity(
        effective_emissivity=0.005,
        layer_material="double_aluminized_mylar",
    )
    print(b2.summary())
    print(f"  Achieved ε*: {b2.effective_emissivity:.6f}")

    print("\n── Example 3: mixed-material blanket ────────────────────────")
    b3 = MLIBlanket.from_layer_types(
        {
            "aluminized_kapton": {"count": 5, "thickness_um": 25},
            "double_aluminized_mylar": {"count": 10, "thickness_um": 6},
            "dacron_net": {"count": 3, "is_spacer": True},
        },
        name="Hybrid Kapton/DAM Blanket",
    )
    print(b3.summary())

    print("\n── Example 4: per-layer emissivity list (7 total layers) ────")
    # [0]=ε₁, [6]=ε₂, [1..5]=5 interior shields
    b4 = MLIBlanket.from_emissivity_list(
        [0.85, 0.03, 0.035, 0.03, 0.04, 0.03, 0.85],
        name="Custom ε-list Blanket",
    )
    print(b4.summary())

    print("\n── Example 5: JSON round-trip ───────────────────────────────")
    json_str = b3.to_json()
    b5 = MLIBlanket.from_json(json.loads(json_str))
    print(f"  Reconstructed ε*: {b5.effective_emissivity:.6f}")
    print(f"  Match: {abs(b3.effective_emissivity - b5.effective_emissivity) < 1e-10}")

    print("\n── Formula verification (25-layer DAM, all same ε=0.03) ─────")
    # layers[0] = layers[24] = ε, 23 interior shields
    # 1/ε* = 1/ε + 1/ε - 1 + 23*(2/ε - 1) = (2/ε - 1)*(1 + 23) = 24*(2/ε - 1)
    N, eps = 25, 0.03
    n_int = N - 2
    expected = 1.0 / (2.0 / eps - 1.0 + 2.0 / eps - 1.0 + n_int * (2.0 / eps - 1.0))
    got = MLIBlanket.from_layer_count(N).effective_emissivity
    print(f"  Manual formula (n={n_int} interior): {expected:.8f}")
    print(f"  Class result                       : {got:.8f}")
    print(f"  Match: {abs(expected - got) < 1e-12}")