"""Screening-level cislunar optical-downlink sizing study.

This script converts payload data rate to physical data rate before sizing the
transmitter. That correction is important: forward-error-correction (FEC)
and protocol overhead consume physical-link capacity and therefore increase
the required received and transmitted optical power.

Sources for numeric inputs (also cited in the accompanying report):
  * Nayak (ed.), 2025, Table 7.1: 0.30 m L1 transmitter, 1.5 m Earth
    receiver, and the 356,510 km L1-to-Earth screening range.
  * Condon et al., 2014, slide 3: 384,400 km Earth-Moon distance and
    64,166 km Moon-to-L2 distance; used to derive the 448,566 km L2
    point-center range.
  * Giggenbach, Knopp, and Fuchs, 2023, Eq. (20) and Table V: 1550 nm,
    250 photons/physical bit, -1 dB transmitter loss, -3 dB pointing loss,
    -1 dB atmospheric loss, and -4.1 dB receiver internal/splitting loss.
  * Downey, 2024: 2/3 is an implemented CCSDS HPE code-rate option.
  * Burleigh et al., 1998: a below-5% protocol-overhead objective; 5% is
    used here as a conservative study allocation.

The code is intentionally a screening model, not a flight link budget. It
uses ideal aperture gain and a clear-sky loss allocation; it does not model
weather, scintillation, telescope efficiency, adaptive modulation/coding,
ground-network availability, or spacecraft electrical/thermal power.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


# Physical constants
SPEED_OF_LIGHT_MPS = 2.998e8
PLANCK_CONSTANT_J_S = 6.626e-34

# Common terminal and payload assumptions
WAVELENGTH_M = 1550e-9
TRANSMITTER_DIAMETER_M = 0.30
RECEIVER_DIAMETER_M = 1.50
TOTAL_PAYLOAD_RATE_BPS = 50.0e9      # desired goal in the solicitation
MINIMUM_PAYLOAD_RATE_BPS = 5.0e9     # minimum requirement in the solicitation

# Receiver/model inputs adapted from Giggenbach et al. (2023).
PHOTONS_PER_PHYSICAL_BIT = 250.0
TRANSMITTER_INTERNAL_LOSS_DB = -1.0
POINTING_LOSS_DB = -3.0
ATMOSPHERIC_LOSS_DB = -1.0
RECEIVER_INTERNAL_LOSS_DB = -4.1

# FEC and overhead selections. The 6 dB margin is a conservative screening
# allocation, rounded from the approximately 5.9 dB example margin reported by
# Giggenbach et al. (2023); it is not a qualified flight margin.
FEC_CODE_RATE = 2.0 / 3.0
PROTOCOL_OVERHEAD_FRACTION = 0.05
LINK_MARGIN_DB = 6.0

# Range inputs. The lunar-orbit value is a rough range proxy only: Nayak's
# Table 7.1 gives 393,860 km for a lunar-surface-to-Earth case, not ELFO.
L1_TO_EARTH_M = 356_510.0e3
LUNAR_ORBIT_PROXY_TO_EARTH_M = 393_860.0e3
L2_TO_EARTH_M = (384_400.0 + 64_166.0) * 1e3


@dataclass(frozen=True)
class ArchitectureCase:
    """A simultaneous direct-to-Earth sharing case used for a sensitivity."""

    name: str
    paths: dict[str, int]
    note: str = ""


PATH_RANGES_M = {
    "L1 halo": L1_TO_EARTH_M,
    "Lunar-orbit range proxy": LUNAR_ORBIT_PROXY_TO_EARTH_M,
    "L2 point-center": L2_TO_EARTH_M,
}

# These are simultaneous usable direct-to-Earth paths, not installed-terminal
# counts. G3's direct L2 leg remains conditional on a GMAT geometry result.
CASES = (
    ArchitectureCase("G1: one L1 gateway", {"L1 halo": 1}),
    ArchitectureCase(
        "G2: one L1 + two lunar-orbit gateways",
        {"L1 halo": 1, "Lunar-orbit range proxy": 2},
    ),
    ArchitectureCase(
        "G3 baseline: one L1 + four lunar-orbit gateways",
        {"L1 halo": 1, "Lunar-orbit range proxy": 4},
        "Five direct paths; no direct L2-to-Earth path assumed.",
    ),
    ArchitectureCase(
        "G3 conditional: add one direct L2 gateway",
        {
            "L1 halo": 1,
            "Lunar-orbit range proxy": 4,
            "L2 point-center": 1,
        },
        "Six direct paths only if GMAT validates the L2 Earth path.",
    ),
    ArchitectureCase("G4: two L1 gateways", {"L1 halo": 2}),
    ArchitectureCase(
        "G5: two L1 + two lunar-orbit gateways",
        {"L1 halo": 2, "Lunar-orbit range proxy": 2},
    ),
)


def watts_to_dbw(power_w: float | np.ndarray) -> float | np.ndarray:
    """Convert optical power in watts to dBW."""
    return 10.0 * np.log10(power_w)


def dbw_to_watts(power_dbw: float | np.ndarray) -> float | np.ndarray:
    """Convert optical power in dBW to watts."""
    return 10.0 ** (power_dbw / 10.0)


def energy_per_photon_j(wavelength_m: float = WAVELENGTH_M) -> float:
    """Return photon energy, hc/lambda, in joules."""
    return PLANCK_CONSTANT_J_S * SPEED_OF_LIGHT_MPS / wavelength_m


def aperture_gain_db(diameter_m: float, wavelength_m: float = WAVELENGTH_M) -> float:
    """Return ideal circular-aperture gain in dB (unity efficiency)."""
    return 20.0 * np.log10(np.pi * diameter_m / wavelength_m)


def free_space_path_gain_db(distance_m: float, wavelength_m: float = WAVELENGTH_M) -> float:
    """Return the negative free-space path gain in dB."""
    return 20.0 * np.log10(wavelength_m / (4.0 * np.pi * distance_m))


def payload_to_physical_rate_bps(payload_rate_bps: float | np.ndarray) -> float | np.ndarray:
    """Apply FEC and protocol overhead to a useful payload data rate."""
    return payload_rate_bps / (FEC_CODE_RATE * (1.0 - PROTOCOL_OVERHEAD_FRACTION))


def physical_to_payload_rate_bps(physical_rate_bps: float | np.ndarray) -> float | np.ndarray:
    """Return the useful payload rate represented by a physical data rate."""
    return physical_rate_bps * FEC_CODE_RATE * (1.0 - PROTOCOL_OVERHEAD_FRACTION)


def required_receiver_power_w(payload_rate_bps: float | np.ndarray) -> float | np.ndarray:
    """Size detector-plane optical power from the FEC/overhead-adjusted rate."""
    physical_rate_bps = payload_to_physical_rate_bps(payload_rate_bps)
    return physical_rate_bps * PHOTONS_PER_PHYSICAL_BIT * energy_per_photon_j()


def total_link_gain_db(distance_m: float) -> float:
    """Return gain/loss from transmitter optical output to detector input."""
    return (
        aperture_gain_db(TRANSMITTER_DIAMETER_M)
        + aperture_gain_db(RECEIVER_DIAMETER_M)
        + free_space_path_gain_db(distance_m)
        + TRANSMITTER_INTERNAL_LOSS_DB
        + POINTING_LOSS_DB
        + ATMOSPHERIC_LOSS_DB
        + RECEIVER_INTERNAL_LOSS_DB
    )


def required_transmitter_power_w(distance_m: float, payload_rate_bps: float | np.ndarray) -> float | np.ndarray:
    """Return spacecraft optical output needed to meet the receiver requirement."""
    receiver_power_dbw = watts_to_dbw(required_receiver_power_w(payload_rate_bps))
    transmitter_power_dbw = receiver_power_dbw - total_link_gain_db(distance_m) + LINK_MARGIN_DB
    return dbw_to_watts(transmitter_power_dbw)


def payload_capacity_bps(distance_m: float, transmitter_power_w: float) -> float:
    """Invert the sizing model to find payload capacity for an optical power class."""
    receiver_power_dbw = watts_to_dbw(transmitter_power_w) + total_link_gain_db(distance_m) - LINK_MARGIN_DB
    receiver_power_w = dbw_to_watts(receiver_power_dbw)
    physical_rate_bps = receiver_power_w / (PHOTONS_PER_PHYSICAL_BIT * energy_per_photon_j())
    return float(physical_to_payload_rate_bps(physical_rate_bps))


def summarize_case(case: ArchitectureCase) -> dict[str, float | str]:
    """Calculate equal-share requirements for one architecture sensitivity."""
    active_paths = sum(case.paths.values())
    payload_per_path_bps = TOTAL_PAYLOAD_RATE_BPS / active_paths
    physical_per_path_bps = payload_to_physical_rate_bps(payload_per_path_bps)

    row: dict[str, float | str] = {
        "case": case.name,
        "active_paths": active_paths,
        "payload_per_path_gbps": payload_per_path_bps / 1e9,
        "physical_per_path_gbps": physical_per_path_bps / 1e9,
        "max_per_terminal_optical_w": 0.0,
        "total_simultaneous_optical_w": 0.0,
        "note": case.note,
    }
    for path_name, count in case.paths.items():
        per_terminal_w = float(required_transmitter_power_w(PATH_RANGES_M[path_name], payload_per_path_bps))
        row[f"{path_name}_per_terminal_w"] = per_terminal_w
        row["max_per_terminal_optical_w"] = max(float(row["max_per_terminal_optical_w"]), per_terminal_w)
        row["total_simultaneous_optical_w"] = float(row["total_simultaneous_optical_w"]) + count * per_terminal_w
    return row


def print_case_summary(rows: list[dict[str, float | str]]) -> None:
    """Print readable results for traceability and use in the report."""
    print("50 Gbps payload sizing case (solicitation desired goal)")
    print(
        f"Physical aggregate rate: {payload_to_physical_rate_bps(TOTAL_PAYLOAD_RATE_BPS) / 1e9:.3f} Gbps "
        f"(FEC={FEC_CODE_RATE:.3f}, overhead={PROTOCOL_OVERHEAD_FRACTION:.0%})"
    )
    for row in rows:
        print(f"\n{row['case']}")
        print(f"  Active direct paths: {row['active_paths']}")
        print(f"  Payload / path: {row['payload_per_path_gbps']:.3f} Gbps")
        print(f"  Physical / path: {row['physical_per_path_gbps']:.3f} Gbps")
        for path_name in PATH_RANGES_M:
            key = f"{path_name}_per_terminal_w"
            if key in row:
                print(f"  {path_name}: {row[key]:.2f} W optical output per terminal")
        print(f"  Maximum terminal output: {row['max_per_terminal_optical_w']:.2f} W")
        print(f"  Sum of simultaneous optical outputs: {row['total_simultaneous_optical_w']:.2f} W")
        if row["note"]:
            print(f"  Note: {row['note']}")


def print_minimum_requirement_case() -> None:
    """Print the single-gateway (G1) output needed for the 5 Gbps minimum."""
    print("\n5 Gbps minimum-requirement check (one gateway carries everything)")
    for path_name, distance_m in PATH_RANGES_M.items():
        power_w = float(required_transmitter_power_w(distance_m, MINIMUM_PAYLOAD_RATE_BPS))
        print(f"  {path_name}: {power_w:.2f} W optical output")


def write_case_csv(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    """Write case results so table values can be checked outside the PDF."""
    output_path = output_dir / "link_budget_case_results.csv"
    fieldnames = sorted({key for row in rows for key in row})
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def plot_power_vs_rate(output_dir: Path) -> Path:
    """Plot required output over the 1--50 Gbps payload-rate range."""
    rates_gbps = np.linspace(1.0, 50.0, 400)
    fig, ax = plt.subplots(figsize=(7.0, 4.4))
    colors = {"L1 halo": "#1565c0", "Lunar-orbit range proxy": "#ef6c00", "L2 point-center": "#6a1b9a"}
    for path_name, distance_m in PATH_RANGES_M.items():
        powers_w = required_transmitter_power_w(distance_m, rates_gbps * 1e9)
        ax.plot(rates_gbps, powers_w, label=path_name, color=colors[path_name], linewidth=2.0)

    for power_w, color in ((75, "#607d8b"), (125, "#455a64"), (225, "#263238")):
        ax.axhline(power_w, color=color, linestyle="--", linewidth=1.0)
        ax.text(50.2, power_w, f"{power_w} W", va="center", fontsize=8, color=color)

    ax.set_xlim(0, 55)
    ax.set_ylim(0, 240)
    ax.set_xlabel("Payload rate carried by one direct optical path (Gbps)")
    ax.set_ylabel("Required spacecraft optical output (W)")
    ax.set_title("Screening optical-output requirement with FEC and overhead")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="upper left", frameon=True)
    fig.tight_layout()

    output_path = output_dir / "optical_power_vs_rate.png"
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return output_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures"),
        help="Directory for the CSV results and PNG figure (default: figures).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rows = [summarize_case(case) for case in CASES]
    print_case_summary(rows)
    print_minimum_requirement_case()
    csv_path = write_case_csv(rows, args.output_dir)
    figure_path = plot_power_vs_rate(args.output_dir)
    print(f"\nWrote {csv_path}")
    print(f"Wrote {figure_path}")


if __name__ == "__main__":
    main()
