# =========================================
# AERO 356 - Assignment 1
#
# Written by Henry Flushman
# =========================================

import numpy as np
import math
from MLI.MLI import MLIBlanket
from heatflux.heatflux import SpacecraftThermal
from rich.console import Console
from rich.table import Table

console = Console()

# ── Helper functions ──────────────────────
def section(title):
    console.rule(title)
# ─────────────────────────────────────────

# ─── Problem 2 ───────────────────────────
section("Problem 2 — Number of Interior Layers in MLI")

emissivity = 0.2
effectiveEmissivity = 0.005
numLayers = MLIBlanket.from_effective_emissivity(effective_emissivity=effectiveEmissivity, layer_emissivity=emissivity).n_layers

table = Table()
table.add_column("Effective Emittance", justify="left")
table.add_column("Emittance", justify="left")
table.add_column("Number of Interior Layers", justify="right")
table.add_row(f"{effectiveEmissivity:.3f}", f"{emissivity:.3f}", f"{numLayers-2}")
console.print(table)
# ─────────────────────────────────────────

# ─── Problem 3 ───────────────────────────
section("Problem 3 — Spacecraft Heat Flux")

sigma = 5.67e-8
diameterSC = 1  # m
r = diameterSC / 2
A_proj = np.pi * r**2
A_surf = 4 * np.pi * r**2
Se = 1366  # W/m^2
bondAlbedo = 0.3  # Earth
tempEarth = 255  # K

sc = SpacecraftThermal(
    name="Spherical Spacecraft",
    area=A_surf,
    emissivity=0.1,
)

# Solar
sc.addSolar(0.1, Se, A_proj, 
            notes="Absorbed solar flux for half of the s/c")

# Albedo
RE = 6378e3
h = 400e3

rho = (RE + h) / RE

G = (2/3) * ((2*rho + rho**-2) - (2 + rho**-2) * math.sqrt(rho**2 - 1))

sc.addAlbedo(0.1, Se, bondAlbedo, A_proj, G, 0.0,
             notes="Earth albedo, using max case (theta_s=0)")

# Earth IR
F = 0.5 * (1 - math.sqrt(1 - (RE / (RE + h))**2))

sc.addEarthIR(0.1, 0.9, tempEarth, F,
              notes="Earth IR using view factor")

# Internal Heat
# sc.addInternal(5, notes="Given heat gain problem")

# Temperature
T = sc.temperatureForNetGain(5.0)

# Output
table = Table(title="Problem 3 — Spacecraft Temperature")

table.add_column("Quantity", justify="left")
table.add_column("Value", justify="right")
table.add_column("Units", justify="left")

table.add_row("Solar Flux", f"{sc.inputs['solar']['Q']:.2f}", "W")
table.add_row("Albedo Flux", f"{sc.inputs['albedo']['Q']:.2f}", "W")
table.add_row("Earth IR Flux", f"{sc.inputs['earth_ir']['Q']:.2f}", "W")
# table.add_row("Internal Heat", f"{sc.inputs['internal']['Q']:.2f}", "W")
table.add_row("Total Input", f"{sc.totalInput():.2f}", "W")

table.add_row("Equilibrium Temp", f"{T:.2f}", "K")
table.add_row("Final Answer (rounded)", f"{round(T)}", "K")

console.print(table)
# ─────────────────────────────────────────

# ─── Problem 4 ───────────────────────────
section("Problem 4 — Flat Plate Temperature")

absorptivity = 0.5
emissivity = 0.1
thetaS = 5 # degrees
# Ignore heat transfer with the Earth

tempEquil = ((absorptivity * Se * np.cos(np.deg2rad(thetaS)))/(emissivity * sigma))**(1/4)

table = Table(title="Problem 4 — Flat Plate Temperature")
table.add_column("Quantity", justify="left")
table.add_column("Value", justify="right")
table.add_column("Units", justify="left")

table.add_row("Solar Absorption", f"{absorptivity}", "N/A")
table.add_row("Emittance", f"{emissivity}", "N/A")
table.add_row("Lambda", f"{thetaS}", "Degrees")
table.add_row("Equilibrium Temperature", f"{round(tempEquil)}", "K")

console.print(table)
# ─────────────────────────────────────────