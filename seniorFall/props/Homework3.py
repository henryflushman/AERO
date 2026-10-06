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
# ║   Course      :  AERO3304 - Propulsions                     ║
# ║                                                             ║
# ║   Assignment  :  Homework 3                                 ║
# ║   Date        :  September 23, 2026                         ║
# ╚═════════════════════════════════════════════════════════════╝


import numpy as np
import matplotlib.pyplot as plt


# === Config ================================

SPECIFIC_IMPULSE = {
    "cold gas": 60,
    "chemical": 320,
    "nuclear": 950,
    "electric": 3000,
}

INERT_MASS_FRACTION = np.linspace(1e-12, 1.0, 1000)

# ===========================================

# === Helper Functions ======================

def nonFeasibleDeltaV(propulsion_type, inert_mass_fraction, g0=9.80665):
    Isp = SPECIFIC_IMPULSE.get(propulsion_type)
    if Isp is None:
        return None
    return Isp * g0 * np.log(1 / inert_mass_fraction)

# ===========================================


# === Analysis ===============================

deltaV = {prop_type: [nonFeasibleDeltaV(prop_type, imf) for imf in INERT_MASS_FRACTION] for prop_type in SPECIFIC_IMPULSE}

plt.figure()
for prop_type, delta_v_values in deltaV.items():
    plt.semilogy(INERT_MASS_FRACTION, delta_v_values, label=prop_type)
plt.axvline(x=0.05, color='r', linestyle='--', label='Structural Limit')
plt.axvline(x=0.4, color='g', linestyle='--', label='Part b')
plt.xlabel("Inert Mass Fraction")
plt.ylabel("Delta V (m/s)")
plt.title("Delta V vs Inert Mass Fraction for Different Propulsion Types")
plt.legend()
plt.grid(True)
plt.show()

print({prop_type: nonFeasibleDeltaV(prop_type, 0.4) for prop_type in SPECIFIC_IMPULSE})

# ============================================