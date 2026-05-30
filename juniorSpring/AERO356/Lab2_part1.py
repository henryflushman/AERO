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
# ║   Course      :  AERO356 - Space Environments II            ║
# ║                                                             ║
# ║   Assignment  :  Lab 2 - Part 1                             ║
# ║   Date        :  May 19, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports =================
import numpy as np
import matplotlib.pyplot as plt
# =============================


# -----------------------------
# Group 1 Data Set: Short-distance aluminum data
# Gap distance = 0.5 in
# Voltage recorded in kV, converted to volts
# -----------------------------
pressure_group1 = np.array([
    1.80, 1.50, 1.30, 1.10, 1.00, 0.92, 0.87, 0.70, 0.65, 0.52, 0.41, 0.31
])

voltage_group1_kv = np.array([
    0.557, 0.525, 0.513, 0.515, 0.540, 0.511, 0.522, 0.538, 0.545, 0.574, 0.623, 0.757
])

voltage_group1 = voltage_group1_kv * 1000

gap_group1 = 0.5
torr_dist_group1 = pressure_group1 * gap_group1


# -----------------------------
# Group 2 Data Set: Short-distance aluminum data
# Gap distance = 1.35 in
# Voltage recorded in kV, converted to volts
# -----------------------------
pressure_group2 = np.array([
    .120, .160, .200, .240, .300, .330, .350, .470, .510, .570, .620, .680, .740, .780, .880, 1.000
])

voltage_group2 = np.array([
    901, 882, 850, 730, 730, 709, 661, 565, 580, 630, 660, 648, 660, 710, 730, 750
])

gap_group2 = 1.35
torr_dist_group2 = pressure_group2 * gap_group2


# -----------------------------
# Group 3 Data Set: Short-distance copper
# Gap distance = 0.5 in
# Voltage recorded in volts
# -----------------------------
pressure_group3 = np.array([
    1.00, 0.88, 0.80, 0.74, 0.64, 0.57, 0.28, 0.24, 0.19, 0.19
])

voltage_group3 = np.array([
    521, 526, 515, 514, 524, 520, 580, 680, 765, 835
])

gap_group3 = 0.5
torr_dist_group3 = pressure_group3 * gap_group3


# -----------------------------
# Group 4 Data Set: Long-distance copper
# Gap distance = 1.35 in
# Voltage recorded in volts
# -----------------------------
pressure_group4 = np.array([
    1.00, 1.10, 2.00, 3.00, 0.50, 0.75, 0.65, 1.30, 0.34, 0.24, 0.15
])

voltage_group4 = np.array([
    764, 760, 860, 930, 670, 675, 680, 760, 660, 705, 830
])

gap_group4 = 1.35
torr_dist_group4 = pressure_group4 * gap_group4


# -----------------------------
# Sort each data set by torr-distance for cleaner connected curves
# -----------------------------
def sort_by_x(x, y):
    idx = np.argsort(x)
    return x[idx], y[idx]

torr_dist_group4, voltage_group4 = sort_by_x(torr_dist_group4, voltage_group4)
torr_dist_group3, voltage_group3 = sort_by_x(torr_dist_group3, voltage_group3)
torr_dist_group2, voltage_group2 = sort_by_x(torr_dist_group2, voltage_group2)
torr_dist_group1, voltage_group1 = sort_by_x(torr_dist_group1, voltage_group1)


# -----------------------------
# Plot Paschen curves
# -----------------------------
plt.figure(figsize=(9, 6))

plt.plot(
    torr_dist_group4,
    voltage_group4,
    marker='o',
    linestyle='--',
    label='Group 4, Copper, d = 1.35 in',
    color='black'
)

plt.plot(
    torr_dist_group3,
    voltage_group3,
    marker='s',
    linestyle='--',
    label='Group 3, Copper, d = 0.5 in',
    color='black'
)

plt.plot(
    torr_dist_group2,
    voltage_group2,
    marker='D',
    linestyle='-',
    label='Group 2, Aluminum, d = 1.35 in',
    color='black'
)

plt.plot(
    torr_dist_group1,
    voltage_group1,
    marker='^',
    linestyle='-',
    label='Group 1, Aluminum, d = 0.5 in',
    color='black'
)

plt.xlabel('Torr-Distance [torr-in]', fontsize=14)
plt.ylabel('Breakdown Voltage [V]', fontsize=14)
plt.title('Paschen Curves', fontsize=20)
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()