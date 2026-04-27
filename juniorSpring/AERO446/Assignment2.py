# =========================================
# AERO 446 - Assignment 2
#
# Written by Henry Flushman
# =========================================

import numpy as np
from rich.console import Console
from rich.table import Table

console = Console()

# ── Helper functions ──────────────────────
def radiationFlux(emissivity, surfaceArea, temperature):
    """Total radiated power (W) from a surface."""
    sigma = 5.670374419e-8  # W/m^2 K^4
    return emissivity * sigma * surfaceArea * temperature**4

def incidentFlux(emissivity, temperature, radiusSource, distance):
    """Incident radiative flux (W/m^2) at distance from a spherical source."""
    sigma = 5.670374419e-8  # W/m^2 K^4
    return emissivity * sigma * temperature**4 * (radiusSource / distance)**2

def equilibriumTemperature(alphaIn, epsilonBlanket, epsilonRadiator, areaIn, areaRadiator, fluxIn):
    """
    Spacecraft equilibrium temperature from absorbed solar power
    = emitted thermal radiation.
    """
    sigma = 5.670374419e-8  # W/m^2 K^4
    return ((alphaIn * fluxIn * areaIn) /
            (sigma * (epsilonBlanket * areaIn + 2 * epsilonRadiator * areaRadiator)))**(1/4)

def requiredRadiatorArea(alphaIn, epsilonBlanket, epsilonRadiator, areaIn, fluxIn, temperature):
    """Required radiator area (m^2) to maintain a given temperature."""
    sigma = 5.670374419e-8  # W/m^2 K^4
    return ((alphaIn * fluxIn * areaIn) / (sigma * temperature**4) - epsilonBlanket * areaIn) / (2 * epsilonRadiator)

def section(title):
    console.rule(title)
# ─────────────────────────────────────────

# ─── Problem 4 ───────────────────────────
section("Problem 4 — Total Radiation Flux from Venus")

emissivity = 1
surfaceArea = 4 * np.pi * (6050e3)**2  # m^2
temperature = 327  # K
flux = radiationFlux(emissivity, surfaceArea, temperature)

table = Table()
table.add_column("Planet", justify="left")
table.add_column("Flux", justify="right")
table.add_column("Units", justify="right")
table.add_row("Venus", f"{flux:.6e}", "W")
console.print(table)
# ─────────────────────────────────────────

# ─── Problem 6 ───────────────────────────
section("Problem 6 — Incident Solar Flux on Earth from the Sun")

emissivity = 1
radiusSun = 695700e3      # m
distanceEarthSun = 149.6e9  # m
temperatureSun = 5800     # K
sunFlux = incidentFlux(emissivity, temperatureSun, radiusSun, distanceEarthSun)

table = Table()
table.add_column("Quantity", justify="left")
table.add_column("Value", justify="right")
table.add_column("Units", justify="right")
table.add_row("Incident Solar Flux", f"{sunFlux:.6e}", "W/m^2")
console.print(table)
# ─────────────────────────────────────────

# ─── Problem 7 ───────────────────────────
section("Problem 7 — Equilibrium Temperature of Geocentric Spacecraft")

areaBlanket = 10  # m^2
areaRadiator = 1  # m^2

# Example cases
alphaBlanket = np.array([1.0, 0.2])
epsilonBlanket = np.array([1.0, 0.2])
epsilonRadiator = np.array([1.0, 0.5])

equilTemp = [0] * 2
for i in range(2):
    equilTemp[i] = equilibriumTemperature(
        alphaBlanket[i],
        epsilonBlanket[i],
        epsilonRadiator[i],
        areaBlanket,
        areaRadiator,
        sunFlux
    )

reqAreaTemp = equilTemp[1] + 10
reqArea = requiredRadiatorArea(
    alphaBlanket[1],
    epsilonBlanket[1],
    epsilonRadiator[1],
    areaBlanket,
    sunFlux,
    reqAreaTemp
)

table = Table()
table.add_column("Case", justify="left")
table.add_column("Value", justify="right")
table.add_column("Units", justify="right")
table.add_row("Ideal Radiator Temperature", f"{equilTemp[0]:.2f}", "K")
table.add_row("Non-Ideal Radiator Temperature", f"{equilTemp[1]:.2f}", "K")
table.add_row("Required Radiator Area", f"{reqArea:.2f}", "m^2")
console.print(table)
# ─────────────────────────────────────────