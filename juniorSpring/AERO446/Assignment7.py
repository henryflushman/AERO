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
# ║   Assignment  :  Homework 7                                 ║
# ║   Date        :  May 15, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# === Constants ===
c = 3e8
T0 = 290.0


# === Helper functions ===

def section(title):
    width = 52
    print(f"\n{'═' * width}")
    print(f"  {title}")
    print(f"{'═' * width}")


def row(label, value, unit=""):
    print(f"  {label:<28} {value:>10.4f}  {unit}")
    

def wavelength(f_hz):
    return c / f_hz


def db_to_linear(x_db):
    return 10 ** (x_db / 10)


def linear_to_db(x):
    return 10 * np.log10(x)


def free_space_path_loss(range_m, f_hz):
    lam = wavelength(f_hz)
    return -20 * np.log10(4 * np.pi * range_m / lam)


def parabolic_gain_linear(area_m2, efficiency, f_hz):
    lam = wavelength(f_hz)
    return 4 * np.pi * efficiency * area_m2 / lam**2


def parabolic_area_from_gain_db(gain_db, efficiency, f_hz):
    lam = c / f_hz
    gain_linear = db_to_linear(gain_db)
    return gain_linear * lam**2 / (4 * np.pi * efficiency)


def diameter_from_area(area_m2):
    return np.sqrt(4 * area_m2 / np.pi)


def radius_from_area(area_m2):
    return np.sqrt(area_m2 / np.pi)


# === Problem 1 ===========================================


# === GIVEN ===
f_p1 = 10e9  # Hz
P_baseline = 100.0
A_baseline = 1.0

P_values = np.arange(10, 201, 10)

link_constant = P_baseline * A_baseline

A_required = link_constant / P_values
D_required = diameter_from_area(A_required)

# Table
df_p1 = pd.DataFrame({
    "Transmit Power (W)": P_values,
    "Required Area (m^2)": A_required,
    "Equivalent Diameter (m)": D_required
})

print("Problem 1 Data Table")
print("--------------------")
print(df_p1.to_string(index=False))
print()

row("Baseline Frequency", f_p1 / 1e9, "GHz")
row("Baseline Power", P_baseline, "W")
row("Baseline Area", A_baseline, "m^2")
row("Link Constant P*A", link_constant, "W*m^2")

# Plot: Antenna area vs. power
plt.figure(figsize=(8, 5))
plt.plot(P_values, A_required, marker="o")
plt.grid(True, linestyle="--", linewidth=0.6)
plt.xlabel("Transmit Power (W)")
plt.ylabel("Required Antenna Area (m²)")
plt.title("Problem 1: Required Antenna Area vs. Transmit Power")
plt.tight_layout()
plt.show()

# Plot: Equivalent dish diameter vs. power
plt.figure(figsize=(8, 5))
plt.plot(P_values, D_required, marker="o")
plt.grid(True, linestyle="--", linewidth=0.6)
plt.xlabel("Transmit Power (W)")
plt.ylabel("Equivalent Parabolic Dish Diameter (m)")
plt.title("Problem 1: Equivalent Dish Diameter vs. Transmit Power")
plt.tight_layout()
plt.show()


# === Problem 2 ===========================================

section("Problem 2: Earth Recieve Antenna Size")



# === GIVEN ===
fc = 4e9
B = 50e6

P_tx_w = 75.0
OBO_db = 3.0
L_line_db = 2.0

R_earth_km = 6378.18
R_geo_km = 42164.18
phi_deg = 60.0

G_rx_existing_db = 21.0
NF_rx_db = 3.0

CN_required_db = 10.0

eta_tx = 0.55
eta_rx = 0.55

L_pol_db = 0.3
L_radome_db = 1.0

L_o2_db_per_km = 0.0050
L_h2o_db_per_km = 0.0015

h_atm_km = 100.0
min_elevation_deg = 5.0



# === Tx Antenna ===

lam = wavelength(fc)

phi_rad = np.radians(phi_deg)

z_km = R_earth_km * np.cos(phi_rad)
y_km = R_earth_km * np.sin(phi_rad)
x_km = R_geo_km - z_km

slant_range_km = np.sqrt(x_km**2 + y_km**2)
slant_range_m = slant_range_km * 1e3

FOV_rad = 2 * np.arctan(y_km / x_km)
FOV_deg = np.rad2deg(FOV_rad)

# BW_deg = 70 * lam / D
D_tx_m = 70 * lam / FOV_deg
A_tx_m2 = np.pi * D_tx_m**2 / 4

G_tx_linear = parabolic_gain_linear(A_tx_m2, eta_tx, fc)
G_tx_db = linear_to_db(G_tx_linear)

P_tx_dbw = linear_to_db(P_tx_w)

transmit_db = P_tx_dbw + G_tx_db - OBO_db - L_line_db

# === Path loss ===

L_fs_db = free_space_path_loss(slant_range_m, fc)

atm_range_km = h_atm_km / np.cos(np.deg2rad(90 - min_elevation_deg))

L_o2_db = L_o2_db_per_km * atm_range_km
L_h2o_db = L_h2o_db_per_km * atm_range_km
L_atm_db = L_o2_db + L_h2o_db

L_path_db = L_fs_db - L_atm_db - L_pol_db - L_radome_db

print(L_path_db)

# === Rx noise temperature ===

F_rx_linear = db_to_linear(NF_rx_db)

T_rx_k = (F_rx_linear - 1) * T0
T_ambient_k = 290.0
T_sys_k = T_rx_k + T_ambient_k

T_sys_db = linear_to_db(T_sys_k)

G_over_T_current_db = G_rx_existing_db - T_sys_db


# === Req G/T ===
B_db = 10 * np.log10(B)

G_over_T_required_db = (
    CN_required_db
    - transmit_db
    - L_path_db
    - 228.6
    + B_db
)

additional_rx_gain_required_db = G_over_T_required_db - G_over_T_current_db

G_rx_antenna_required_db = additional_rx_gain_required_db

A_rx_required_m2 = parabolic_area_from_gain_db(
    G_rx_antenna_required_db,
    eta_rx,
    fc
)

r_rx_required_m = radius_from_area(A_rx_required_m2)
D_rx_required_m = diameter_from_area(A_rx_required_m2)


# -----------------------------
# Results
# -----------------------------

print("Problem 2 Input Values")
print("----------------------")
row("Carrier Frequency", fc / 1e9, "GHz")
row("Bandwidth", B / 1e6, "MHz")
row("Transmit Power", P_tx_w, "W")
row("Transmit Power", P_tx_dbw, "dBW")
row("Output Back Off", OBO_db, "dB")
row("Line Loss", L_line_db, "dB")
row("Receiver Gain", G_rx_existing_db, "dB")
row("Receiver Noise Figure", NF_rx_db, "dB")
row("Required C/N", CN_required_db, "dB")
print()

print("Orbit and Transmit Antenna")
print("--------------------------")
row("z", z_km, "km")
row("y", y_km, "km")
row("x", x_km, "km")
row("Slant Range", slant_range_km, "km")
row("Required FOV", FOV_deg, "deg")
row("TX Antenna Diameter", D_tx_m, "m")
row("TX Antenna Area", A_tx_m2, "m^2")
row("TX Antenna Gain", G_tx_db, "dB")
row("Transmit Term", transmit_db, "dB")
print()

print("Path Losses")
print("-----------")
row("Free Space Path Loss", L_fs_db, "dB")
row("Atmospheric Range", atm_range_km, "km")
row("Oxygen Loss", L_o2_db, "dB")
row("Water Vapor Loss", L_h2o_db, "dB")
row("Atmospheric Loss", L_atm_db, "dB")
row("Polarization Loss", L_pol_db, "dB")
row("Radome Loss", L_radome_db, "dB")
row("Total Path Term", L_path_db, "dB")
print()

print("Receiver Noise")
print("--------------")
row("Receiver Noise Factor", F_rx_linear, "")
row("Receiver Noise Temp", T_rx_k, "K")
row("Ambient/Antenna Temp", T_ambient_k, "K")
row("System Noise Temp", T_sys_k, "K")
row("System Noise Temp", T_sys_db, "dBK")
row("Current G/T", G_over_T_current_db, "dB/K")
print()

print("Required Receive Antenna")
print("------------------------")
row("Required G/T", G_over_T_required_db, "dB/K")
row("Additional RX Gain Needed", additional_rx_gain_required_db, "dB")
row("RX Antenna Gain Required", G_rx_antenna_required_db, "dB")
row("Required RX Antenna Area", A_rx_required_m2, "m^2")
row("Required RX Antenna Radius", r_rx_required_m, "m")
row("Required RX Antenna Diameter", D_rx_required_m, "m")
print()

print("Final Answer")
print("------------")
print(f"Required Earth receive antenna area     = {A_rx_required_m2:.3f} m^2")
print(f"Required Earth receive antenna radius   = {r_rx_required_m:.3f} m")
print(f"Required Earth receive antenna diameter = {D_rx_required_m:.3f} m")
print()