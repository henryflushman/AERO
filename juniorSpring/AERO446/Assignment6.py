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
# ║   Assignment  :  Homework 5                                 ║
# ║   Date        :  May 15, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# === Constants ================================================
c = 3e8


# ── Helper functions ──────────────────────
def section(title):
    width = 52
    print(f"\n{'═' * width}")
    print(f"  {title}")
    print(f"{'═' * width}")

def row(label, value, unit=""):
    print(f"  {label:<28} {value:>10.4f}  {unit}")

def wavelength(f_hz):
    return c / f_hz

def parabolic_gain_linear(area_m2, efficiency, f_hz):
    lam = wavelength(f_hz)
    a_eff = efficiency * area_m2
    return 4 * np.pi * a_eff / lam**2

def linear_to_db(x):
    return 10 * np.log10(x)

def free_space_path_loss_db(range_m, f_hz):
    lam = wavelength(f_hz)
    return -20 * np.log10(4 * np.pi * range_m / lam)
# ─────────────────────────────────────────


# === Problem 1 ===

# Given
D = 1.0         # antenna diameter
f_nom = 300e6   # nominal frequency
eta = 0.60      # efficiency

area = np.pi * D**2 / 4

# Analysis
G1_linear = parabolic_gain_linear(area, eta, f_nom)
G1_db = linear_to_db(G1_linear)

print("Problem 1")
print("---------")
print(f"Antenna area = {area:.4f} m^2")
print(f"Wavelength = {wavelength(f_nom):.4f} m")
print(f"Gain = {G1_linear:.4f} linear")
print(f"Gain = {G1_db:.2f} dB")
print()


# === Problem 2 ===

f_low = 0.95 * f_nom
f_high = 1.05 * f_nom

# Analysis
G_low_db = linear_to_db(parabolic_gain_linear(area, eta, f_low))
G_nom_db = linear_to_db(parabolic_gain_linear(area, eta, f_nom))
G_high_db = linear_to_db(parabolic_gain_linear(area, eta, f_high))

print("Problem 2")
print("---------")
print(f"Low frequency = {f_low/1e6:.1f} MHz")
print(f"Nominal frequency = {f_nom/1e6:.1f} MHz")
print(f"High frequency = {f_high/1e6:.1f} MHz")
print()
print(f"Gain at {f_low/1e6:.1f} MHz = {G_low_db:.3f} dB")
print(f"Gain at {f_nom/1e6:.1f} MHz = {G_nom_db:.3f} dB")
print(f"Gain at {f_high/1e6:.1f} MHz = {G_high_db:.3f} dB")
print()
print(f"Change below nominal = {G_low_db - G_nom_db:.3f} dB")
print(f"Change above nominal = {G_high_db - G_nom_db:.3f} dB")
print()
print("The gain is not exactly symmetric because G is proportional to f^2.")
print()


# === Problem 3 ===


# Freq initialize
freqs_hz = []

for decade_start in [1e8, 1e9, 1e10]:
    for multiplier in range(1, 10):
        freqs_hz.append(multiplier * decade_start)
        
freqs_hz.append(1e11)

freqs_hz = np.array(freqs_hz)
freqs_ghz = freqs_hz / 1e9

# Given
area_p3 = 1.0
eta_p3 = 1.0
range_m = 1000e3

gain_p3_linear = parabolic_gain_linear(area_p3, eta_p3, freqs_hz)
gain_p3_db = linear_to_db(gain_p3_linear)

fspl_db = free_space_path_loss_db(range_m, freqs_hz)

# Table
df = pd.DataFrame({
    "Frequency (Hz)": freqs_hz,
    "Frequency (GHz)": freqs_ghz,
    "Gain (dB)": gain_p3_db,
    "Free Space Path Loss (dB)": fspl_db
})

print("Problem 3 Data Table")
print("--------------------")
print(df.to_string(index=False))
print()


# -----------------------------
# Plot 1: Gain vs frequency
# -----------------------------
plt.figure(figsize=(8, 5))
plt.semilogx(freqs_ghz, gain_p3_db, marker="o")
plt.grid(True, which="both", linestyle="--", linewidth=0.6)
plt.xlabel("Frequency (GHz)")
plt.ylabel("Gain (dB)")
plt.title("Parabolic Reflector Gain vs. Frequency\nA = 1 m², Efficiency = 100%")
plt.tight_layout()
plt.show()


# -----------------------------
# Plot 2: Free-space path loss
# -----------------------------
plt.figure(figsize=(8, 5))
plt.semilogx(freqs_ghz, fspl_db, marker="o")
plt.grid(True, which="both", linestyle="--", linewidth=0.6)
plt.xlabel("Frequency (GHz)")
plt.ylabel("Free Space Path Loss (dB)")
plt.title("Free Space Path Loss vs. Frequency\nSlant Range = 1,000 km")
plt.tight_layout()
plt.show()


# -----------------------------
# Written conclusions
# -----------------------------
print("Problem 3 Conclusions")
print("---------------------")
print("Gain increases with frequency because G is proportional to 1/lambda^2.")
print("Since lambda = c/f, gain is proportional to f^2.")
print("Therefore, gain increases by about 20 dB per decade of frequency.")
print()
print("Free-space path loss also increases with frequency.")
print("For fixed range, FSPL is proportional to f^2 in linear scale,")
print("so it also increases by about 20 dB per decade.")