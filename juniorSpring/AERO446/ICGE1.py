import numpy as np
from rich.console import Console
from rich.table import Table

console = Console()

SAWE_MGA_HIGH = {
    "estimated_layout": {
        "electronics": {
            (0, 5): 35,
            (5, 15): 25,
            (15, float("inf")): 20,
        },
        "structure": 25,
        "thermal": 50,
        "propulsion": 25,
        "batteries": 25,
        "wire harness": 100,
        "solar array": 35,
        "eclss crew": 30,
        "brackets clips fasteners": 35,
        "mechanisms": 25,
        "instrumentation": 75,
    },
    "layout": {
        "electronics": {
            (0, 5): 30,
            (5, 15): 20,
            (15, float("inf")): 15,
        },
        "structure": 20,
        "thermal": 30,
        "propulsion": 20,
        "batteries": 20,
        "wire harness": 45,
        "solar array": 20,
        "eclss crew": 20,
        "brackets clips fasteners": 25,
        "mechanisms": 20,
        "instrumentation": 30,
    },
    "preliminary_design": {
        "electronics": {
            (0, 5): 20,
            (5, 15): 15,
            (15, float("inf")): 12,
        },
        "structure": 15,
        "thermal": 15,
        "propulsion": 15,
        "batteries": 15,
        "wire harness": 25,
        "solar array": 15,
        "eclss crew": 15,
        "brackets clips fasteners": 15,
        "mechanisms": 15,
        "instrumentation": 25,
    },
    "released_design": {
        "electronics": {
            (0, 5): 10,
            (5, 15): 10,
            (15, float("inf")): 10,
        },
        "structure": 6,
        "thermal": 8,
        "propulsion": 7,
        "batteries": 7,
        "wire harness": 10,
        "solar array": 5,
        "eclss crew": 8,
        "brackets clips fasteners": 8,
        "mechanisms": 4,
        "instrumentation": 5,
    },
    "existing_hardware": {
        "electronics": {
            (0, 5): 5,
            (5, 15): 3,
            (15, float("inf")): 3,
        },
        "structure": 3,
        "thermal": 3,
        "propulsion": 3,
        "batteries": 3,
        "wire harness": 5,
        "solar array": 3,
        "eclss crew": 4,
        "brackets clips fasteners": 5,
        "mechanisms": 3,
        "instrumentation": 3,
    },
}

cubesatProperties = np.array([
    ["structure", 3.7, "layout"],
    ["thermal", 0.7, "estimated_layout"],
    ["propulsion", 1.0, "layout"],
    ["electronics", 2.5, "nte"],
    ["batteries", 0.75, "estimated_layout"],
    ["wire harness", 0.15, "estimated_layout"],
    ["mechanisms", 0.10, "existing_hardware"],
], dtype=object)


def get_mga_percent(subsystem, mass, maturity):
    if maturity in {"nte", "specification_value", "spec", "cfe"}:
        return 0.0

    entry = SAWE_MGA_HIGH[maturity][subsystem]

    if subsystem == "electronics":
        for (low_mass, high_mass), percent in entry.items():
            if low_mass <= mass < high_mass:
                return percent

    return entry


table = Table(title="CubeSat Mass Growth Allowance")
table.add_column("Subsystem", justify="left")
table.add_column("Base Mass [kg]", justify="right")
table.add_column("Maturity", justify="left")
table.add_column("MGA %", justify="right")
table.add_column("MGA Mass [kg]", justify="right")
table.add_column("Total Mass [kg]", justify="right")

total_base_mass = 0.0
total_mga_mass = 0.0
total_mass_with_mga = 0.0

for subsystem, mass, maturity in cubesatProperties:
    mass = float(mass)

    mga_percent = get_mga_percent(subsystem, mass, maturity)
    mga_mass = mass * mga_percent / 100.0
    total_mass = mass + mga_mass

    total_base_mass += mass
    total_mga_mass += mga_mass
    total_mass_with_mga += total_mass

    table.add_row(
        subsystem,
        f"{mass:.3f}",
        maturity,
        f"{mga_percent:.1f}",
        f"{mga_mass:.3f}",
        f"{total_mass:.3f}",
    )

table.add_section()
table.add_row(
    "[bold]TOTAL[/bold]",
    f"[bold]{total_base_mass:.3f}[/bold]",
    "-",
    "-",
    f"[bold]{total_mga_mass:.4f}[/bold]",
    f"[bold]{total_mass_with_mga:.4f}[/bold]",
)

console.print(table)
console.print()
console.print(f"[bold]Nominal estimated mass:[/bold] {total_base_mass:.3f} kg")
console.print(f"[bold]Total MGA added:[/bold] {total_mga_mass:.4f} kg")
console.print(f"[bold]Worst-case mass with MGA:[/bold] {total_mass_with_mga:.3f} kg")