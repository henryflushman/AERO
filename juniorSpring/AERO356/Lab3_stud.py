"""
Langmuir probe I-V analysis script for Lab 3 Data.xlsx.

This script:
1. Loads voltage/current data from your Excel file.
2. Uses the exact spreadsheet columns:
      Voltage (V)
      Current (A)
      neg current
3. Uses "neg current" for the Langmuir probe analysis.
4. Smooths the analysis current using a Savitzky-Golay filter.
5. Calculates:
      - floating potential
      - ion saturation current
      - electron temperature
      - plasma potential
      - electron saturation current
      - plasma density
      - Debye length
      - plasma frequency
6. Creates plots:
      - Current vs Voltage
      - Negative Current vs Voltage
      - Raw and smoothed Langmuir I-V curve
      - I-V curve with analysis lines
      - Semilog electron-current plot
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter


# ------------------------------------------------------------
# USER SETTINGS
# ------------------------------------------------------------

# This assumes Lab3.py is in:
# juniorSpring/AERO356/
#
# and your Excel file is in:
# juniorSpring/AERO356/labData/
DATA_FILE = Path(__file__).resolve().parent / "labData" / "Lab 3 Data.xlsx"
SHEET_NAME = "Sheet1"


ANALYSIS_CURRENT_COLUMN = "neg current"

# Probe geometry from your notes.
PROBE_DIAMETER_CM = 0.0508
PROBE_EXPOSED_LENGTH_CM = 0.1905

# Xenon atomic mass in amu.
XENON_ATOMIC_MASS_AMU = 131.1

# Smoothing settings.
USE_SMOOTHED_DATA_FOR_ANALYSIS = True

# Increase this for more smoothing: 11, 13, 15, etc.
# Your student data only has about 25 points, so do not go too huge.
SMOOTH_WINDOW = 9
SMOOTH_POLY_ORDER = 3

# Region choices for your student data.
# Adjust these after looking at your plots.
FLOATING_FIT_BOUNDS = (-5, 6)
ION_SAT_V_BOUNDS = (-100, -80)
RETARDING_REGION_BOUNDS = (0, 20)
PLASMA_TANGENT_BOUNDS = (10, 50)
ELECTRON_SAT_BOUNDS = (60, 100)


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

def normalize_name(name):
    """Normalize column names for easier matching."""
    return str(name).strip().lower().replace("_", " ")


def get_column(df, target_name):
    """
    Find a column by name, ignoring capitalization and extra spaces.
    """
    target = normalize_name(target_name)

    for col in df.columns:
        if normalize_name(col) == target:
            return col

    raise KeyError(
        f"Could not find column '{target_name}'.\n"
        f"Available columns are:\n{list(df.columns)}"
    )


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
        raise ValueError("Lines are nearly parallel; plasma potential is unreliable.")

    return (b2 - b1) / (m1 - m2)


def require_points(mask, name, min_points=2):
    """Make sure a selected fitting region has enough points."""
    count = int(np.sum(mask))

    if count < min_points:
        raise ValueError(
            f"Not enough points in {name}. "
            f"Found {count}, need at least {min_points}. "
            f"Adjust that voltage range."
        )


def make_valid_savgol_window(requested_window, n_points, poly_order):
    """
    Savitzky-Golay window must be:
    - odd
    - greater than poly_order
    - less than or equal to number of data points
    """
    window = min(requested_window, n_points)

    if window % 2 == 0:
        window -= 1

    min_window = poly_order + 2

    if min_window % 2 == 0:
        min_window += 1

    if window < min_window:
        window = min_window

    if window > n_points:
        raise ValueError("Not enough data points for smoothing settings.")

    return window


# ------------------------------------------------------------
# LOAD DATA
# ------------------------------------------------------------

df = pd.read_excel(DATA_FILE, sheet_name=SHEET_NAME)

# Clean column names, but keep original text for lookup.
df.columns = [str(c).strip() for c in df.columns]

print("\n--- Spreadsheet Columns Found ---")
for col in df.columns:
    print(col)

voltage_col = get_column(df, "Voltage (V)")
current_col = get_column(df, "Current (A)")
neg_current_col = get_column(df, "neg current")
analysis_col = get_column(df, ANALYSIS_CURRENT_COLUMN)

V = pd.to_numeric(df[voltage_col], errors="coerce").to_numpy()
I_current = pd.to_numeric(df[current_col], errors="coerce").to_numpy()
I_neg_current = pd.to_numeric(df[neg_current_col], errors="coerce").to_numpy()
I_raw_analysis = pd.to_numeric(df[analysis_col], errors="coerce").to_numpy()

# Remove rows with invalid numbers.
valid = (
    np.isfinite(V)
    & np.isfinite(I_current)
    & np.isfinite(I_neg_current)
    & np.isfinite(I_raw_analysis)
)

V = V[valid]
I_current = I_current[valid]
I_neg_current = I_neg_current[valid]
I_raw_analysis = I_raw_analysis[valid]

# Sort by voltage.
order = np.argsort(V)
V = V[order]
I_current = I_current[order]
I_neg_current = I_neg_current[order]
I_raw_analysis = I_raw_analysis[order]


# ------------------------------------------------------------
# SMOOTH DATA
# ------------------------------------------------------------

smooth_window = make_valid_savgol_window(
    requested_window=SMOOTH_WINDOW,
    n_points=len(I_raw_analysis),
    poly_order=SMOOTH_POLY_ORDER
)

I_smooth = savgol_filter(
    I_raw_analysis,
    window_length=smooth_window,
    polyorder=SMOOTH_POLY_ORDER
)

if USE_SMOOTHED_DATA_FOR_ANALYSIS:
    I_analysis = I_smooth
else:
    I_analysis = I_raw_analysis


# ------------------------------------------------------------
# PROBE AREA AND ION MASS
# ------------------------------------------------------------

d_m = PROBE_DIAMETER_CM * 1e-2
L_m = PROBE_EXPOSED_LENGTH_CM * 1e-2
r_m = d_m / 2

# Cylindrical side area:
# A = 2*pi*r*L = pi*d*L
A_probe = 2 * np.pi * r_m * L_m

m_i = XENON_ATOMIC_MASS_AMU * amu

print("\n--- Data Settings ---")
print(f"Data file = {DATA_FILE}")
print(f"Analysis current column = {analysis_col}")
print(f"Using smoothed data for analysis = {USE_SMOOTHED_DATA_FOR_ANALYSIS}")
print(f"Smoothing window used = {smooth_window}")
print(f"Smoothing polynomial order = {SMOOTH_POLY_ORDER}")

print("\n--- Probe / Gas Properties ---")
print(f"Probe diameter = {PROBE_DIAMETER_CM:.6f} cm")
print(f"Probe exposed length = {PROBE_EXPOSED_LENGTH_CM:.6f} cm")
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
print(f"Floating potential V_f = {V_f:.6f} V")


# ------------------------------------------------------------
# ION SATURATION CURRENT
# ------------------------------------------------------------

# Use the low-voltage ion saturation region.
ion_sat_mask = between(V, ION_SAT_V_BOUNDS)
require_points(ion_sat_mask, "ion saturation region")

I_ion_sat_measured = np.mean(I_analysis[ion_sat_mask])
I_is = abs(I_ion_sat_measured)

print("\n--- Ion Saturation Current ---")
print(f"Ion saturation voltage range = {ION_SAT_V_BOUNDS}")
print(f"Average ion saturation current = {I_ion_sat_measured:.6e} A")
print(f"I_is magnitude = {I_is:.6e} A")
print(f"Number of points used = {np.sum(ion_sat_mask)}")


# ------------------------------------------------------------
# ELECTRON CURRENT
# ------------------------------------------------------------

# The measured probe current contains ion current and electron current.
# Approximate ion current as the ion-saturation baseline, then subtract.
I_e = I_analysis - I_ion_sat_measured

# Log requires positive current.
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
print(f"Retarding fit voltage range = {RETARDING_REGION_BOUNDS}")
print(f"Slope of ln(I_e) vs V = {m_log:.6e} 1/V")
print(f"T_e = {T_e_eV:.6f} eV")
print(f"T_e = {T_e_K:.6e} K")


# ------------------------------------------------------------
# PLASMA POTENTIAL
# ------------------------------------------------------------

# Plasma potential is estimated by intersecting:
# 1. a linear fit through the steep part of the curve
# 2. a linear fit through the high-voltage electron-current region
steep_mask = between(V, PLASMA_TANGENT_BOUNDS)
esat_mask = between(V, ELECTRON_SAT_BOUNDS)

require_points(steep_mask, "plasma tangent region")
require_points(esat_mask, "electron saturation region")

m_steep, b_steep = linear_fit(V[steep_mask], I_analysis[steep_mask])
m_esat, b_esat = linear_fit(V[esat_mask], I_analysis[esat_mask])

V_p = line_intersection(m_steep, b_steep, m_esat, b_esat)

print("\n--- Plasma Potential ---")
print(f"Plasma tangent voltage range = {PLASMA_TANGENT_BOUNDS}")
print(f"Electron saturation voltage range = {ELECTRON_SAT_BOUNDS}")
print(f"Plasma potential V_p = {V_p:.6f} V")


# ------------------------------------------------------------
# ELECTRON SATURATION CURRENT
# ------------------------------------------------------------

I_es_mean = np.mean(I_analysis[esat_mask])

print("\n--- Electron Saturation Current ---")
print(f"Mean electron saturation current = {I_es_mean:.6e} A")


# ------------------------------------------------------------
# PLASMA DENSITY, DEBYE LENGTH, PLASMA FREQUENCY
# ------------------------------------------------------------

# Bohm speed:
# c_s = sqrt(k*T_e/m_i)
#
# Since T_e is in eV:
# k*T_e = e*T_e_eV
c_s = np.sqrt(e * T_e_eV / m_i)

# Bohm ion current:
# I_is = 0.6*e*n_i*A_probe*c_s
#
# Solve for n_i:
n_i = I_is / (0.6 * e * A_probe * c_s)
n_e = n_i

# Debye length:
lambda_D = np.sqrt(eps0 * e * T_e_eV / (n_e * e**2))

# Electron plasma frequency:
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


# Plot 1: Current vs Voltage, matching spreadsheet column C.
plt.figure()
plt.plot(V, I_current, "o-", markersize=4, label="Current (A)")
plt.xlabel("Voltage [V]")
plt.ylabel("Current [A]")
plt.title("Current vs Voltage")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("current_vs_voltage.png", dpi=200)


# Plot 2: Negative Current vs Voltage, matching spreadsheet column D.
plt.figure()
plt.plot(V, I_neg_current, "o-", markersize=4, label="neg current")
plt.xlabel("Voltage [V]")
plt.ylabel("Negative Current [A]")
plt.title("Negative Current vs Voltage")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("negative_current_vs_voltage.png", dpi=200)


# Plot 3: raw and smoothed analysis current.
plt.figure()
plt.plot(V, I_raw_analysis, "o", markersize=4, alpha=0.55, label="Raw analysis current")
plt.plot(V, I_smooth, "-", linewidth=2, label="Smoothed analysis current")
plt.xlabel("Voltage [V]")
plt.ylabel("Current [A]")
plt.title("Raw and Smoothed Langmuir Probe I-V Curve")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("iv_curve_raw_and_smooth.png", dpi=200)


# Plot 4: I-V curve with analysis lines.
plt.figure()
plt.plot(V, I_raw_analysis, "o", markersize=4, alpha=0.45, label="Raw analysis current")
plt.plot(V, I_analysis, "-", linewidth=2, label="Analysis current")
plt.axhline(I_ion_sat_measured, linestyle="--", label="Ion saturation avg")
plt.axvline(V_f, linestyle="--", label=f"Floating potential = {V_f:.2f} V")
plt.axvline(V_p, linestyle="--", label=f"Plasma potential = {V_p:.2f} V")
plt.plot(V_fit, m_steep * V_fit + b_steep, label="Steep-region fit")
plt.plot(V_fit, m_esat * V_fit + b_esat, label="Electron saturation fit")
plt.xlabel("Voltage [V]")
plt.ylabel("Current [A]")
plt.title("Langmuir Probe I-V Analysis")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("iv_curve_analysis.png", dpi=200)


# Plot 5: semilog electron current for electron temperature.
plt.figure()
plt.semilogy(
    V[positive_mask],
    I_e[positive_mask],
    "o",
    markersize=4,
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

plt.xlabel("Voltage [V]")
plt.ylabel("Electron current [A]")
plt.title("Electron Retardation Region Semilog Plot")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.savefig("electron_temperature_fit.png", dpi=200)


plt.show()