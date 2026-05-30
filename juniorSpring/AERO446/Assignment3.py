# =========================================
# AERO 446 - Assignment 3
#
# Written by Henry Flushman
# =========================================

import numpy as np
from functions import truePowerFromRef, temperatureCoefficent

# ── Helper functions ──────────────────────
def section(title):
    width = 52
    print(f"\n{'═' * width}")
    print(f"  {title}")
    print(f"{'═' * width}")

def row(label, value, unit=""):
    print(f"  {label:<28} {value:>10.4f}  {unit}")
# ─────────────────────────────────────────

# ─────────────────────────────────────────
# Problem 1 — Solar Power Sizing
# ─────────────────────────────────────────
section("Problem 1 — Solar Power Sizing")

# Subsystem power requirements (W)
subsysPower = {
    "Payload":  138,
    "Structure": 3,
    "Thermal":   10,
    "Power":     5,
    "TTC":       36,
    "Computer":  5,
    "ADCS":      30,
}

# LEO orbit
orbitalPeriod = 90 * 60   # s

# Temperature range
tempDay     =  67   # C
tempEclipse = -65   # C

# Spacecraft spends one-third of orbit eclipsed
timeEclipse = orbitalPeriod / 3               # s
timeDay     = orbitalPeriod - timeEclipse     # s

# Imaging payload operation
payloadTimeDay     = 2 * 60   # s
payloadTimeEclipse = 3 * 60   # s

# ADCS operation
adcsTimeDay     = 30   # s
adcsTimeEclipse = 30   # s

# Power efficiency
effDay     = 0.8
effEclipse = 0.9

# Solar cell values
effSolarCell     = 0.302       # efficiency
tempRef          = 28          # C  (standard reference temperature)
Voc              = 0.99        # V  (open-circuit voltage)
l, w             = 6, 5        # cm (cell dimensions)
refUnitSolarCell = 135.3       # mW/cm²
dVdT             = -6.7        # mV/C (temperature coefficient)

# ─────────────────────────────────────────
# a) Temperature Coefficient
# ─────────────────────────────────────────

solarCellTC = temperatureCoefficent(Voc=Voc, dVdT=dVdT/1000)

row("Temperature Coefficient", solarCellTC, "%/C")

# ─────────────────────────────────────────
# b) Average Power per Orbit
# ─────────────────────────────────────────

basePower = (subsysPower["Structure"]
             + subsysPower["Thermal"]
             + subsysPower["Power"]
             + subsysPower["TTC"]
             + subsysPower["Computer"])

energyDay     = basePower * timeDay     + subsysPower["Payload"] * payloadTimeDay     + subsysPower["ADCS"] * adcsTimeDay
energyEclipse = basePower * timeEclipse + subsysPower["Payload"] * payloadTimeEclipse + subsysPower["ADCS"] * adcsTimeEclipse

avgPowerDay     = energyDay     / timeDay
avgPowerEclipse = energyEclipse / timeEclipse

avgPowerOrbitReq = ((avgPowerEclipse * timeEclipse / effEclipse)
                  + (avgPowerDay     * timeDay     / effDay)) / timeDay

row("Avg Power (day)",     avgPowerDay,      "W")
row("Avg Power (eclipse)", avgPowerEclipse,  "W")
row("Avg Power (orbit)",   avgPowerOrbitReq, "W")

# ─────────────────────────────────────────
# c) Power Generated per Cell
# ─────────────────────────────────────────

cellArea     = l * w                                          # cm²
refCellPower = refUnitSolarCell * cellArea / 1000             # W  (at 28 C)
dayCellPower = truePowerFromRef(refPower=refCellPower,
                                tempRef=tempRef,
                                tempCoef=solarCellTC/100,
                                tempActual=tempDay)           # W
orbitPowerGen = dayCellPower * timeDay / 3600                 # W·hr per orbit

row("Cell Area",            cellArea,      "cm²")
row("Ref Cell Power",       refCellPower,  "W")
row("Daylight Cell Power",  dayCellPower,  "W")
row("Power Gen per Orbit",  orbitPowerGen, "W·hr")

# ─────────────────────────────────────────
# d) Number of Solar Cells Required
# ─────────────────────────────────────────

nSolarCell = np.ceil(avgPowerOrbitReq / orbitPowerGen)

row("Solar Cells Required", nSolarCell, "cells")

# ─────────────────────────────────────────
# e) Area per Solar Panel (2 panels)
# ─────────────────────────────────────────

totalArea = nSolarCell * cellArea
panelArea = totalArea / 2

row("Total Cell Area",  totalArea, "cm²")
row("Area per Panel",   panelArea, "cm²")

print(f"\n{'═' * 52}\n")