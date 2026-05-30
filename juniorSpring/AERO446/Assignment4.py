# ╔═════════════════════════════════════════════════════════════╗
# ║                                                             ║
# ║        .o.       oooooooooooo ooooooooo.     .oooooo.       ║
# ║       .888.      `888'     `8 `888   `Y88.  d8P'  `Y8b      ║
# ║      .8"888.      888          888   .d88' 888      888     ║
# ║     .8' `888.     888oooo8     888ooo88P'  888      888     ║
# ║    .88ooo8888.    888    "     888`88b.    888      888     ║
# ║   .8'     `888.   888       o  888  `88b.  `88b    d88'     ║
# ║  o88o     o8888o o888ooooood8 o888o  o888o  `Y8bood8P'      ║
# ║                                                             ║
# ║                ── CALIFORNIA POLYTECHNIC ──                 ║
# ║                                                             ║
# ╠═════════════════════════════════════════════════════════════╣
# ║   Author      :  Henry Flushman                             ║
# ║   Course      :  AERO446 - Spacecraft Electrical and        ║
# ║                                     Electric Systems        ║
# ║   Assignment  :  Homework 4                                 ║
# ║   Date        :  May 6, 2026                                ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
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
# Problem 1 — Battery Capacity
# ─────────────────────────────────────────
section("Problem 1 — Battery Capacity")

# === Problem givens ===
T = 6*3600      # s
powerReq = 500  # W
DOD = 0.5
betaAngle = 0.0 # rad

def radiusFromPeriod(period, mu=398600):
    """
    Calculates the radius of a circular orbit, given a specified period length in seconds
    """
    return ((mu*period**2)/(4*np.pi**2))**(1/3)

def timeEclipse(R_orbit, betaAngle, R_earth=6378.18):
    """
    Calculates the total time the spacecraft spends eclipsed each orbit
    """
    return 0.5 + (1/np.pi)*np.arcsin((1-(R_earth/R_orbit)**2)**(1/2)/np.cos(betaAngle))

def fractionEclipse(R_orbit, betaAngle, R_earth=6378.18):
    """
    Calculates the fractional time in eclipse
    """
    return (1/np.pi)*np.arcsin(((R_earth/R_orbit)**2-np.sin(betaAngle)**2)**(1/2)/np.cos(betaAngle))

def batteryCapacity(P_out, t_discharge, DOD, eff=1):
    """
    Calculates the required battery capacity given the time in eclipse, depth of discharge,
    efficiency of the battery, and total power requirements
    """
    return (P_out*t_discharge)/(DOD*eff)

# Solve for orbit radius
R_orbit = radiusFromPeriod(T)

# Solve for time in eclipse
T_eclipse = timeEclipse(R_orbit, betaAngle)

# Solve for fractional time in eclipse
frac_eclipse = fractionEclipse(R_orbit, betaAngle)
T_eclipse_fractional = frac_eclipse*T

# Solve for battery capacity
E_battery = batteryCapacity(powerReq, T_eclipse_fractional/3600, DOD)

row("Fractional Eclipse Time", T_eclipse_fractional/60, "minutes")
row("Battery Capacity", E_battery, "A-hr")


# ─────────────────────────────────────────
# Problem 2 — Solar Panel Surface Area
# ─────────────────────────────────────────
section("Problem 2 — Solar Panel Surface Area")

# === Problem Givens ===
powerGen_solarcell = 245    # W/m2

powerGen = (powerReq*T_eclipse_fractional)/(T-T_eclipse_fractional)

areaSolarPanel = (powerGen+powerReq)/powerGen_solarcell

row("Total Power Required", powerGen+powerReq, "W")
row("Solar Panel Area Required", areaSolarPanel, "m^2")


# ─────────────────────────────────────────
# Problem 3 — Added Panel Area for Heater
# ─────────────────────────────────────────
section("Problem 3 — Added Panel Area for Heater")

powerGen_heater = ((powerReq+100)*T_eclipse_fractional)/(T-T_eclipse_fractional)

areaSolarPanel_heater = (powerGen_heater+powerReq)/powerGen_solarcell

row("Total Power Required with Heater", powerGen_heater+powerReq, "W")
row("Additional Panel Area", areaSolarPanel_heater-areaSolarPanel, "m^2")