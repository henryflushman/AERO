import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Problem 1 - AERO 446
# Spacecraft EPS / Solar Array / Battery / Cost Analysis
# ============================================================

# -----------------------------
# Constants
# -----------------------------
Re = 6378.0                 # Earth radius, km
mu = 398600.4418            # Earth gravitational parameter, km^3/s^2
P0 = 367.0                  # Solar cell output at normal incidence, W/m^2

beta = 0.0                  # Beta angle, radians
altitude = np.arange(300, 2000 + 20, 20, dtype=float)   # km
r = Re + altitude           # orbital radius, km

# -----------------------------
# Power requirements
# -----------------------------
P_payload = 138.0
P_structure = 30.0
P_thermal = 80.0
P_power = 25.0
P_comms = 40.0
P_obc = 25.0

# Payload operates 2 min around noon.
# Data downlink operates 5 min during eclipse/night.
t_payload = 2.0 * 60.0      # seconds
t_downlink = 5.0 * 60.0     # seconds

# All other systems active all the time
P_base = P_structure + P_thermal + P_power + P_obc

# ============================================================
# Orbit timing
# ============================================================

# Orbital period
T = 2.0 * np.pi * np.sqrt(r**3 / mu)      # seconds

# Eclipse fraction for circular orbit:
# f_e = (1/pi) asin( sqrt((Re/r)^2 - sin^2(beta)) / cos(beta) )
sin_beta = np.sin(beta)
cos_beta = np.cos(beta)

beta_crit = np.arcsin(Re / r)
has_eclipse = np.abs(beta) < beta_crit

f_eclipse = np.zeros_like(r)
arg = np.zeros_like(r)

arg[has_eclipse] = np.sqrt((Re / r[has_eclipse])**2 - sin_beta**2) / cos_beta
arg = np.clip(arg, 0.0, 1.0)

f_eclipse[has_eclipse] = (1.0 / np.pi) * np.arcsin(arg[has_eclipse])

t_eclipse = T * f_eclipse
t_daylight = T - t_eclipse

# Fixed one-sided anti-Earth panel only generates when its front side faces the Sun.
# For beta = 0 and radial panel orientation, this is half of the orbit.
t_generation = T / 2.0

# ============================================================
# Daylight and eclipse power
# ============================================================

# Average daylight power including 2-min payload operation
P_daylight_avg = (P_base * t_daylight + P_payload * t_payload) / t_daylight

# Average eclipse power including 5-min downlink operation
P_eclipse_avg = (P_base * t_eclipse + P_comms * t_downlink) / t_eclipse

# Energy per orbit
E_daylight_Wh = (P_base * t_daylight + P_payload * t_payload) / 3600.0
E_eclipse_Wh = (P_base * t_eclipse + P_comms * t_downlink) / 3600.0
E_total_Wh = E_daylight_Wh + E_eclipse_Wh

# ============================================================
# Solar panel area
# ============================================================

# The fixed one-sided anti-Earth panel generates power according to:
#
# P_solar(theta) = P0 * A * max(cos(theta - pi/2), 0)
#
# theta = 0 deg   -> Point A, bottom of orbit
# theta = 90 deg  -> noon, panel points most directly at Sun
# theta = 270 deg -> midnight / eclipse center
#
# Integral from 0 to 2pi of max(cos(theta - pi/2), 0) dtheta = 2.
#
# Since dt = T/(2pi) dtheta:
#
# E_solar_per_m2 = P0 * T/(2pi) * 2 / 3600
#                = P0 * T/(pi * 3600)

E_solar_per_m2_Wh = P0 * T / (np.pi * 3600.0)

# Solar array must cover daylight loads directly and must also recharge
# the battery energy used during eclipse. Because the battery is not 100%
# efficient, the array must generate E_eclipse_Wh / battery_efficiency.

E_solar_needed_Wh = E_daylight_Wh + E_eclipse_Wh

A_required = E_solar_needed_Wh / E_solar_per_m2_Wh

# ============================================================
# 10-year degraded solar area
# ============================================================

degradation_rate = 0.02
lifetime_years = 10.0

degradation_factor = (1.0 - degradation_rate) ** lifetime_years
A_required_10yr = A_required / degradation_factor

# ============================================================
# Battery capacity
# ============================================================

battery_efficiency = 0.80
max_DOD = 0.50

battery_capacity_Wh = E_eclipse_Wh / (battery_efficiency * max_DOD)

# ============================================================
# Mission cost
# ============================================================

C_baseline = 3_000_000.0
k_altitude = 50_000.0
k_solar = 5_000.0
a = 1.2
b = 1.5

# Use 10-year degraded solar panel area for mission cost.
# Change this to A_required if your instructor wants non-degraded area in the cost model.
A_for_cost = A_required_10yr

C_mission = C_baseline + k_altitude * altitude**a + k_solar * A_for_cost**b

idx_opt = np.argmin(C_mission)

h_opt = altitude[idx_opt]
C_min = C_mission[idx_opt]
A_opt = A_required[idx_opt]
A_opt_10yr = A_required_10yr[idx_opt]
battery_opt = battery_capacity_Wh[idx_opt]

# ============================================================
# Print results
# ============================================================

print("========== Optimal Mission Design ==========")
print(f"Optimal altitude: {h_opt:.0f} km")
print(f"Minimum mission cost: ${C_min:,.2f}")
print()
print("========== Values at Optimal Altitude ==========")
print(f"Orbital period: {T[idx_opt] / 60:.4f} min")
print(f"Eclipse time: {t_eclipse[idx_opt] / 60:.4f} min")
print(f"Daylight time: {t_daylight[idx_opt] / 60:.4f} min")
print(f"Power generation time: {t_generation[idx_opt] / 60:.4f} min")
print(f"Average daylight power required: {P_daylight_avg[idx_opt]:.4f} W")
print(f"Average eclipse power required: {P_eclipse_avg[idx_opt]:.4f} W")
print(f"Solar panel area, no degradation: {A_opt:.4f} m^2")
print(f"Solar panel area, 10-year degraded: {A_opt_10yr:.4f} m^2")
print(f"Battery capacity required: {battery_opt:.4f} W-hr")

# ============================================================
# Plot helper
# ============================================================

def make_plot(x, y, xlabel, ylabel, title):
    plt.figure()
    plt.plot(x, y, linewidth=1.8)
    plt.grid(True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.tight_layout()

# ============================================================
# 1. Eclipse Time vs Altitude
# ============================================================

make_plot(
    altitude,
    t_eclipse / 60.0,
    "Orbital Altitude (km)",
    "Eclipse Time (min)",
    "Eclipse Time vs Orbital Altitude"
)

# ============================================================
# 2. Daylight Time vs Altitude
# ============================================================

make_plot(
    altitude,
    t_daylight / 60.0,
    "Orbital Altitude (km)",
    "Daylight Time (min)",
    "Daylight Time vs Orbital Altitude"
)

# ============================================================
# 3. Power Generation Time vs Altitude
# ============================================================

make_plot(
    altitude,
    t_generation / 60.0,
    "Orbital Altitude (km)",
    "Power Generation Time (min)",
    "Power Generation Time vs Orbital Altitude"
)

# ============================================================
# 4. Power Required During Daylight vs Altitude
# ============================================================

make_plot(
    altitude,
    P_daylight_avg,
    "Orbital Altitude (km)",
    "Average Daylight Power Required (W)",
    "Average Daylight Power Required vs Orbital Altitude"
)

# ============================================================
# 5. Power Required During Eclipse vs Altitude
# ============================================================

make_plot(
    altitude,
    P_eclipse_avg,
    "Orbital Altitude (km)",
    "Average Eclipse Power Required (W)",
    "Average Eclipse Power Required vs Orbital Altitude"
)

# ============================================================
# 6. Solar Panel Area Required vs Altitude
# ============================================================

make_plot(
    altitude,
    A_required,
    "Orbital Altitude (km)",
    "Solar Panel Area Required (m$^2$)",
    "Solar Panel Area Required vs Orbital Altitude"
)

# ============================================================
# 7. 10-Year Degraded Solar Panel Area vs Altitude
# ============================================================

make_plot(
    altitude,
    A_required_10yr,
    "Orbital Altitude (km)",
    "10-Year Solar Panel Area Required (m$^2$)",
    "10-Year Degraded Solar Panel Area vs Orbital Altitude"
)

# ============================================================
# 8. Battery Capacity vs Altitude
# ============================================================

make_plot(
    altitude,
    battery_capacity_Wh,
    "Orbital Altitude (km)",
    "Battery Capacity Required (W-hr)",
    "Battery Capacity vs Orbital Altitude"
)

# ============================================================
# 9. Mission Cost vs Altitude
# ============================================================

make_plot(
    altitude,
    C_mission / 1e6,
    "Orbital Altitude (km)",
    "Mission Cost ($M)",
    "Mission Cost vs Orbital Altitude"
)

# Mark optimum on cost plot
plt.scatter(h_opt, C_min / 1e6, zorder=5)
plt.annotate(
    f"Min: {h_opt:.0f} km\n${C_min/1e6:.2f}M",
    xy=(h_opt, C_min / 1e6),
    xytext=(h_opt + 120, C_min / 1e6),
    arrowprops=dict(arrowstyle="->")
)

# ============================================================
# 11. Solar Power Generated vs Spacecraft Position
# ============================================================

theta_deg = np.linspace(0.0, 360.0, 1000)
theta = np.deg2rad(theta_deg)

# theta = 0 deg starts at Point A, bottom of orbit.
# theta = 90 deg is noon, where the anti-Earth panel points toward the Sun.
# theta = 270 deg is midnight/eclipse center.
#
# Use cosine measured relative to the noon direction:
# incidence = cos(theta - 90 deg)
incidence_factor = np.maximum(np.cos(theta - np.pi / 2.0), 0.0)

# Eclipse mask at optimal altitude
r_opt = Re + h_opt
alpha_eclipse = np.arcsin(Re / r_opt)   # eclipse half-angle, rad

def angle_difference(angle, center):
    return np.abs(np.arctan2(np.sin(angle - center), np.cos(angle - center)))

eclipse_center = 3.0 * np.pi / 2.0
in_eclipse = angle_difference(theta, eclipse_center) <= alpha_eclipse

# Power generated by the 10-year-sized panel
P_solar_BOL = P0 * A_opt_10yr * incidence_factor
P_solar_EOL = degradation_factor * P_solar_BOL

# Force power to zero during eclipse
P_solar_BOL[in_eclipse] = 0.0
P_solar_EOL[in_eclipse] = 0.0

plt.figure()
plt.plot(theta_deg, P_solar_BOL, linewidth=1.8, label="Beginning of Life")
plt.plot(theta_deg, P_solar_EOL, "--", linewidth=1.8, label="End of Life after 10 Years")
plt.grid(True)
plt.xlabel("Spacecraft Position from Point A (deg)")
plt.ylabel("Solar Power Generated (W)")
plt.title(f"Solar Power Generated Around Orbit at h = {h_opt:.0f} km")
plt.xlim(0, 360)
plt.legend()
plt.tight_layout()

plt.show()