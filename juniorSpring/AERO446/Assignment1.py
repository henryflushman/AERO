# =========================================
# AERO 446 - Assignment 1
#
# Written by Henry Flushman
# =========================================

import numpy as np

# ── Helper functions ──────────────────────
def FOR(r, Re=6378):
    """Field of Regard (degrees)"""
    return np.rad2deg(2 * np.arcsin(Re / r))

def T(r, mu=398600):
    """Orbital period (seconds)"""
    return 2 * np.pi * np.sqrt(r**3 / mu)

def section(title):
    width = 52
    print(f"\n{'═' * width}")
    print(f"  {title}")
    print(f"{'═' * width}")

def row(label, value, unit=""):
    print(f"  {label:<28} {value:>10.4f}  {unit}")
# ─────────────────────────────────────────

# ─────────────────────────────────────────
# Problem 1 — LEO Constellation Comparison
# ─────────────────────────────────────────
section("Problem 1 — LEO Constellations")

satellites = {
    "Kuiper Phase 2":  450,
    "OneWeb":          1200,
    "O3b mPower":      8000,
}

col = [28, 12, 14, 14]
header = f"  {'Constellation':<{col[0]}} {'Alt (km)':>{col[1]}} {'FOR (deg)':>{col[2]}} {'Period (min)':>{col[3]}}"
divider = "  " + "─" * (sum(col) + 3)

print(f"\n{divider}")
print(header)
print(divider)
for name, alt in satellites.items():
    r = alt + 6378
    print(f"  {name:<{col[0]}} {alt:>{col[1]}.0f} {FOR(r):>{col[2]}.4f} {T(r)/60:>{col[3]}.4f}")
print(divider)

# ─────────────────────────────────────────
# Problem 2 — GEO Orbit (Earth)
# ─────────────────────────────────────────
section("Problem 2 — GEO Orbit (Earth)")

sidereal_day = 86164.1   # s
mu_earth     = 398600    # km³/s²

r_GEO = ((mu_earth * sidereal_day**2) / (4 * np.pi**2))**(1/3)

row("GEO Altitude",         r_GEO - 6378, "km")
row("Field of Regard (FOR)", FOR(r_GEO),   "deg")

# ─────────────────────────────────────────
# Problem 3 — Synchronous Orbit (Mars)
# ─────────────────────────────────────────
section("Problem 3 — Synchronous Orbit (Mars)")

sidereal_day_mars = 88642.7   # s
mu_mars           = 4.28284e4  # km³/s²
Re_mars           = 3390       # km

r_MSO = ((mu_mars * sidereal_day_mars**2) / (4 * np.pi**2))**(1/3)

row("MSO Altitude",          r_MSO - Re_mars,       "km")
row("Field of Regard (FOR)", FOR(r_MSO, Re=Re_mars), "deg")

# ─────────────────────────────────────────
# Problem 4 — WorldView-3 FOV
# ─────────────────────────────────────────
section("Problem 4 — WorldView-3 FOV")

alt_WV3 = 617   # km
sw_WV3  = 13    # km  (swath width)

fov = np.rad2deg(2 * np.arctan(sw_WV3 / (2 * alt_WV3)))

row("Altitude",      alt_WV3, "km")
row("Swath Width",   sw_WV3,  "km")
row("FOV",           fov,     "deg")

print(f"\n{'═' * 52}\n")