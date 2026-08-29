"""
Langmuir probe I-V analysis script.

This script:
1. Loads voltage/current data from Excel.
2. Optionally flips the current sign.
3. Smooths the current data using a Savitzky-Golay filter.
4. Estimates:
   - ion saturation current
   - floating potential
   - electron temperature
   - plasma potential
   - electron saturation current
   - plasma density
   - Debye length
   - plasma frequency
5. Creates useful plots for the lab report.

Adjust the voltage ranges near the top after looking at your plots.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------

DATA_FILE = Path("juniorSpring/AERO356/labData/LangmuirProfessional.xlsx")
SHEET_NAME = "Sheet1"

# If the curve is upside down compared to the class-note curve, set this to True.
FLIP_CURRENT_SIGN = False

# Use smoothed data for calculations.
# Set to False if you want calculations based on raw data.
USE_SMOOTHED_DATA_FOR_ANALYSIS = True

# Probe geometry from your notes.
PROBE_DIAMETER_CM = 0.0508
PROBE_EXPOSED_LENGTH_CM = 0.1905

# Xenon atomic mass in amu.
XENON_ATOMIC_MASS_AMU = 131.1

# Tolerance used to identify ion saturation region.
# Class note method: take min current, add small number, average all values smaller.
ION_SAT_TOLERANCE = 0.0001

# Savitzky-Golay smoothing settings.
# Window length must be odd.
SMOOTH_WINDOW = 31
SMOOTH_POLY_ORDER = 3

# Region II: log-linear electron retardation region.
# Adjust this after looking at the semilog plot.
RETARDING_REGION_BOUNDS = (10, 30)

# Window near I = 0 for floating potential fit.
FLOATING_FIT_BOUNDS = (-3, 3)

# Steep part of I-V curve used for plasma potential tangent fit.
PLASMA_TANGENT_BOUNDS = (25, 38)

# Electron saturation region.
ELECTRON_SAT_BOUNDS = (55, 75)


# ------------------------------------------------------------
# CONSTANTS
# ------------------------------------------------------------

e = 1.602176634e-19          # elementary charge, C
k_B = 1.380649e-23           # Boltzmann constant, J/K
eps0 = 8.8541878128e-12      # vacuum permittivity, F/m
m_e = 9.1093837015e-31       # electron mass, kg
amu = 1.66053906660e-27      # atomic mass unit, kg


# ------------------------------------------------------------
# HELPER FUNCTIONS
# ------------------------------------------------------------

def between(x, bounds):
    """Return Boolean mask for values inside inclusive bounds."""
    lo, hi = bounds
    return (x >= lo) & (x <= hi)


def linear_fit(x, y):
    """
    Fit y = m*x + b.
    Returns m, b.
    """
    m, b = np.polyfit(x, y, 1)
    return m, b


def line_intersection(m1, b1, m2, b2):
    """
    Solve for x where:
        m1*x + b1 = m2*x + b2
    """
    if np.isclose(m1, m2):
        raise ValueError("Lines are nearly parallel; plasma potential intersection is unreliable.")

    return (b2 - b1) / (m1 - m2)


def require_points(mask, name, min_points=2):
    """Make sure a selected fitting region has enough points."""
    count = np.sum(mask)

    if count < min_points:
        raise ValueError(
            f"Not enough points in {name}. "
            f"Found {count}, need at least {min_points}. "
            f"Adjust the voltage bounds for that region."
        )


def make_valid_savgol_window(requested_window, data_length, poly_order):
    """
    Savitzky-Golay window must be:
    - odd
    - greater than poly_order
    - less than or equal to number of data points
    """
    window = requested_window

    if window > data_length:
        window = data_length

    if window % 2 == 0:
        window -= 1

    minimum_window = poly_order + 2

    if minimum_window % 2 == 0:
        minimum_window += 1

    if window < minimum_window:
        window = minimum_window

    if window > data_length:
        raise ValueError("Not enough data points for the requested smoothing settings.")

    return window


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

df = pd.read_excel(DATA_FILE, sheet_name=SHEET_NAME)

# Clean column names.
df.columns = [str(c).strip().lower() for c in df.columns]

# Expect columns called voltage and current.
V = df["voltage"].to_numpy(dtype=float)
I_raw = df["current"].to_numpy(dtype=float)

# Remove bad rows.
valid = np.isfinite(V) & np.isfinite(I_raw)
V = V[valid]
I_raw = I_raw[valid]

# Sort by voltage.
order = np.argsort(V)
V = V[order]
I_raw = I_raw[order]

# Flip current sign if needed.
if FLIP_CURRENT_SIGN:
    I_raw = -1 * I_raw


# ------------------------------------------------------------
# SMOOTH CURRENT DATA
# ------------------------------------------------------------

smooth_window = make_valid_savgol_window(
    requested_window=SMOOTH_WINDOW,
    data_length=len(I_raw),
    poly_order=SMOOTH_POLY_ORDER
)

I_smooth = savgol_filter(
    I_raw,
    window_length=smooth_window,
    polyorder=SMOOTH_POLY_ORDER
)

if USE_SMOOTHED_DATA_FOR_ANALYSIS:
    I_analysis = I_smooth
else:
    I_analysis = I_raw

print("\n--- Data Settings ---")
print(f"Current sign flipped = {FLIP_CURRENT_SIGN}")
print(f"Using smoothed data for analysis = {USE_SMOOTHED_DATA_FOR_ANALYSIS}")
print(f"Smoothing window used = {smooth_window}")
print(f"Smoothing polynomial order = {SMOOTH_POLY_ORDER}")


# ------------------------------------------------------------
# PROBE AREA AND ION MASS
# ------------------------------------------------------------

# Convert probe dimensions from cm to m.
d_m = PROBE_DIAMETER_CM * 1e-2
L_m = PROBE_EXPOSED_LENGTH_CM * 1e-2
r_m = d_m / 2

# Cylindrical side area:
# A = 2*pi*r*L = pi*d*L
A_probe = 2 * np.pi * r_m * L_m

# Xenon ion mass.
m_i = XENON_ATOMIC_MASS_AMU * amu

print("\n--- Probe / Gas Properties ---")
print(f"Probe diameter = {PROBE_DIAMETER_CM:.6f} cm")
print(f"Probe length = {PROBE_EXPOSED_LENGTH_CM:.6f} cm")
print(f"Probe area = {A_probe:.6e} m^2")
print(f"Xenon ion mass = {m_i:.6e} kg")


# ------------------------------------------------------------
# FLOATING POTENTIAL
# ------------------------------------------------------------

# Floating potential occurs where I = 0.
floating_mask = between(V, FLOATING_FIT_BOUNDS)
require_points(floating_mask, "floating potential fit region")

m_float, b_float = linear_fit(V[floating_mask], I_analysis[floating_mask])
V_f = -b_float / m_float

print("\n--- Floating Potential ---")
print(f"Floating potential, V_f = {V_f:.6f} V")


# ------------------------------------------------------------
# ION SATURATION CURRENT
# ------------------------------------------------------------

# Class-note method:
# Find current min, add tolerance, average all current values below threshold.
I_min = np.min(I_analysis)
I_cutoff = I_min + ION_SAT_TOLERANCE

ion_sat_mask = I_analysis <= I_cutoff
require_points(ion_sat_mask, "ion saturation region")

I_ion_sat_measured = np.mean(I_analysis[ion_sat_mask])
I_is = abs(I_ion_sat_measured)

print("\n--- Ion Saturation Current ---")
print(f"Minimum current = {I_min:.6e} A")
print(f"Ion saturation cutoff = {I_cutoff:.6e} A")
print(f"Average ion saturation current = {I_ion_sat_measured:.6e} A")
print(f"I_is magnitude = {I_is:.6e} A")
print(f"Number of points used = {np.sum(ion_sat_mask)}")


# ------------------------------------------------------------
# ELECTRON CURRENT
# ------------------------------------------------------------

# The measured current includes ion current and electron current.
# The ion current baseline is approximately the ion saturation current.
# Subtract it to isolate electron current.
#
# If I_ion_sat_measured is negative:
# I_e = I - (-I_is) = I + I_is
I_e = I_analysis - I_ion_sat_measured

# Keep positive values only because log(current) requires positive current.
positive_mask = I_e > 0


# ------------------------------------------------------------
# ELECTRON TEMPERATURE
# ------------------------------------------------------------

# In the electron retardation region:
#
# I_e ~ exp(V / T_e)
#
# ln(I_e) = (1 / T_e) V + b
#
# Therefore:
#
# T_e in eV = 1 / slope

retarding_mask = between(V, RETARDING_REGION_BOUNDS) & positive_mask
require_points(retarding_mask, "electron retardation region")

ln_Ie = np.log(I_e[retarding_mask])

m_log, b_log = linear_fit(V[retarding_mask], ln_Ie)

T_e_eV = 1 / abs(m_log)
T_e_K = T_e_eV * e / k_B

print("\n--- Electron Temperature ---")
print(f"Slope of ln(I_e) vs V = {m_log:.6e} 1/V")
print(f"T_e = {T_e_eV:.6f} eV")
print(f"T_e = {T_e_K:.6e} K")


# ------------------------------------------------------------
# PLASMA POTENTIAL
# ------------------------------------------------------------

# Plasma potential is found by intersecting:
# 1. a tangent/linear fit through the steep electron retardation region
# 2. a linear fit through the electron saturation region

steep_mask = between(V, PLASMA_TANGENT_BOUNDS)
esat_mask = between(V, ELECTRON_SAT_BOUNDS)

require_points(steep_mask, "plasma tangent region")
require_points(esat_mask, "electron saturation region")

m_steep, b_steep = linear_fit(V[steep_mask], I_analysis[steep_mask])
m_esat, b_esat = linear_fit(V[esat_mask], I_analysis[esat_mask])

V_p = line_intersection(m_steep, b_steep, m_esat, b_esat)

print("\n--- Plasma Potential ---")
print(f"V_p = {V_p:.6f} V")


# ------------------------------------------------------------
# ELECTRON SATURATION CURRENT
# ------------------------------------------------------------

I_es_mean = np.mean(I_analysis[esat_mask])

print("\n--- Electron Saturation Current ---")
print(f"Mean electron saturation current = {I_es_mean:.6e} A")


# ------------------------------------------------------------
# PLASMA DENSITY, DEBYE LENGTH, PLASMA FREQUENCY
# ------------------------------------------------------------

# Bohm ion speed:
#
# c_s = sqrt(k*T_e / m_i)
#
# Since T_e is in eV:
#
# k*T_e = e*T_e_eV
#
# So:
#
# c_s = sqrt(e*T_e_eV / m_i)

c_s = np.sqrt(e * T_e_eV / m_i)

# Bohm ion current:
#
# I_is = 0.6 * e * n_i * A_probe * c_s
#
# Solve:
#
# n_i = I_is / (0.6 * e * A_probe * c_s)

n_i = I_is / (0.6 * e * A_probe * c_s)
n_e = n_i

# Debye length:
#
# lambda_D = sqrt(eps0 * k*T_e / (n_e * e^2))
#
# Since k*T_e = e*T_e_eV:
#
# lambda_D = sqrt(eps0 * e*T_e_eV / (n_e * e^2))

lambda_D = np.sqrt(eps0 * e * T_e_eV / (n_e * e**2))

# Electron plasma angular frequency:
#
# omega_pe = sqrt(n_e * e^2 / (eps0 * m_e))

omega_pe = np.sqrt(n_e * e**2 / (eps0 * m_e))
f_pe = omega_pe / (2 * np.pi)

print("\n--- Plasma Density and Derived Parameters ---")
print(f"Bohm speed = {c_s:.6e} m/s")
print(f"n_i = n_e = {n_i:.6e} m^-3")
print(f"n_i = n_e = {n_i / 1e6:.6e} cm^-3")
print(f"Debye length = {lambda_D:.6e} m")
print(f"Electron plasma angular frequency = {omega_pe:.6e} rad/s")
print(f"Electron plasma frequency = {f_pe:.6e} Hz")


# ------------------------------------------------------------
# PLOTS
# ------------------------------------------------------------

V_fit = np.linspace(np.min(V), np.max(V), 500)


# Plot 1: raw and smoothed I-V curve.
plt.figure()
plt.plot(V, I_raw, ".", markersize=3, alpha=0.45, label="Raw I-V data")
plt.plot(V, I_smooth, "-", linewidth=2, label="Smoothed I-V data")
plt.xlabel("Voltage, V")
plt.ylabel("Current, A")
plt.title("Raw and Smoothed Langmuir Probe I-V Curve")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("iv_curve_raw_and_smooth.png", dpi=200)


# Plot 2: I-V curve with analysis lines.
plt.figure()
plt.plot(V, I_raw, ".", markersize=3, alpha=0.35, label="Raw I-V data")
plt.plot(V, I_analysis, "-", linewidth=2, label="Analysis data")
plt.axhline(I_ion_sat_measured, linestyle="--", label="Ion saturation avg")
plt.axvline(V_f, linestyle="--", label=f"Floating potential = {V_f:.2f} V")
plt.axvline(V_p, linestyle="--", label=f"Plasma potential = {V_p:.2f} V")
plt.plot(V_fit, m_steep * V_fit + b_steep, label="Steep-region fit")
plt.plot(V_fit, m_esat * V_fit + b_esat, label="Electron saturation fit")
plt.xlabel("Voltage, V")
plt.ylabel("Current, A")
plt.title("Langmuir Probe I-V Analysis")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("iv_curve_analysis.png", dpi=200)


# Plot 3: semilog electron current for electron temperature.
plt.figure()
plt.semilogy(
    V[positive_mask],
    I_e[positive_mask],
    ".",
    markersize=3,
    label="Electron current"
)

V_log_fit = np.linspace(
    RETARDING_REGION_BOUNDS[0],
    RETARDING_REGION_BOUNDS[1],
    200
)

Ie_log_fit = np.exp(m_log * V_log_fit + b_log)

plt.semilogy(
    V_log_fit,
    Ie_log_fit,
    linewidth=2,
    label=f"Retardation fit, Te = {T_e_eV:.2f} eV"
)

plt.axvspan(
    RETARDING_REGION_BOUNDS[0],
    RETARDING_REGION_BOUNDS[1],
    alpha=0.15,
    label="Fit region"
)

plt.xlabel("Voltage, V")
plt.ylabel("Electron current, A")
plt.title("Electron Retardation Region Semilog Plot")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig("electron_temperature_fit.png", dpi=200)


plt.show()