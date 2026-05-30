import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Radiation Disk Analysis
# 1) Plot distance versus Counts Per Minute (CPM)
# 2) Plot shield thickness versus CPM for each material on one set of axes
# 3) Calculate dose in rad/hr at a 3 inch distance
#
# Given conversions:
#   1 CPM = 0.001 mR/hr
#   1 R/hr = 0.877 rad/hr
# ============================================================

# ============================================================
# PUT YOUR TWO CSV FILE PATHS HERE
# Example:
# distance_csv = r"C:\Users\YourName\Documents\labData\Lab1b_radiation.csv"
# shield_csv   = r"C:\Users\YourName\Documents\labData\Lab1c_radShield.csv"
# ============================================================
distance_csv = r"juniorSpring/AERO356/labData/Lab1c_varyDist.csv"
shield_csv   = r"juniorSpring/AERO356/labData/Lab1c_radShield.csv"


# -----------------------------
# Load CSV files
distance_data = pd.read_csv(distance_csv)
shield_data = pd.read_csv(shield_csv)

# Clean column names in case they include quotes/spaces
distance_data.columns = distance_data.columns.str.strip().str.replace('"', '', regex=False)
shield_data.columns = shield_data.columns.str.strip().str.replace('"', '', regex=False)

# Clean string values in shield material column
if "Shield Material" in shield_data.columns:
    shield_data["Shield Material"] = (
        shield_data["Shield Material"]
        .astype(str)
        .str.strip()
        .str.replace('"', '', regex=False)
    )

# Convert numeric columns
distance_data["Time [s]"] = pd.to_numeric(distance_data["Time [s]"], errors="coerce")
distance_data["Distance [in]"] = pd.to_numeric(distance_data["Distance [in]"], errors="coerce")
distance_data["Total Counts"] = pd.to_numeric(distance_data["Total Counts"], errors="coerce")

shield_data["Time [s]"] = pd.to_numeric(shield_data["Time [s]"], errors="coerce")
shield_data["Distance [in]"] = pd.to_numeric(shield_data["Distance [in]"], errors="coerce")
shield_data["Number of Shields"] = pd.to_numeric(shield_data["Number of Shields"], errors="coerce")
shield_data["Thickness of Shield"] = pd.to_numeric(shield_data["Thickness of Shield"], errors="coerce")
shield_data["Total Count"] = pd.to_numeric(shield_data["Total Count"], errors="coerce")

# Drop rows with missing required numeric data
distance_data = distance_data.dropna(subset=["Time [s]", "Distance [in]", "Total Counts"]).copy()
shield_data = shield_data.dropna(
    subset=["Time [s]", "Thickness of Shield", "Total Count", "Shield Material"]
).copy()

# -----------------------------
# Compute CPM
distance_data["CPM"] = distance_data["Total Counts"] / (distance_data["Time [s]"] / 60.0)
shield_data["CPM"] = shield_data["Total Count"] / (shield_data["Time [s]"] / 60.0)

# -----------------------------
# Plot 1: Distance vs CPM
plt.figure(figsize=(8, 5))

plt.plot(
    distance_data["Distance [in]"],
    distance_data["CPM"],
    linestyle="-",
    marker="o",
    markevery=1,       # distance dataset is small, so every point is okay
    linewidth=1.2,
    markersize=5,
    color="black",
)

plt.xlabel("Distance [in]")
plt.ylabel("Counts Per Minute (CPM)")
plt.title("Distance vs Counts Per Minute (CPM)")
plt.grid(True, linestyle=":", linewidth=0.5)
plt.tight_layout()
plt.show()


# -----------------------------
# Plot 2: Shield Thickness vs CPM for all materials on one set of axes
plt.figure(figsize=(10, 6))

line_styles = ["-", "-", "-", "-"]
markers = ["o", "s", "^", "v", "D", "x", "+"]

materials = shield_data["Shield Material"].dropna().unique()

for i, material in enumerate(materials):
    material_data = shield_data[shield_data["Shield Material"] == material].copy()
    material_data = material_data.sort_values("Thickness of Shield")

    plt.plot(
        material_data["Thickness of Shield"],
        material_data["CPM"],
        linestyle=line_styles[i % len(line_styles)],
        marker=markers[i % len(markers)],
        linewidth=1.2,
        markersize=5,
        color="black",
        label=material
    )

plt.xlabel("Shield Thickness")
plt.ylabel("Counts Per Minute (CPM)")
plt.title("Shield Thickness vs CPM for Each Shielding Material")
plt.legend()
plt.grid(True, linestyle=":", linewidth=0.5)
plt.tight_layout()
plt.show()

# -----------------------------
# Dose calculation at 3 inches
row_3in = distance_data[np.isclose(distance_data["Distance [in]"], 3.0)]

if not row_3in.empty:
    cpm_3in = row_3in["CPM"].iloc[0]

    # 1 CPM = 0.001 mR/hr
    dose_mR_per_hr = cpm_3in * 0.001

    # Convert mR/hr to R/hr
    dose_R_per_hr = dose_mR_per_hr / 1000.0

    # 1 R/hr = 0.877 rad/hr
    dose_rad_per_hr = dose_R_per_hr * 0.877

    print("----- Dose Calculation at 3 inches -----")
    print(f"CPM at 3 inches: {cpm_3in:.3f} CPM")
    print(f"Exposure rate: {dose_mR_per_hr:.6f} mR/hr")
    print(f"Exposure rate: {dose_R_per_hr:.9f} R/hr")
    print(f"Dose rate: {dose_rad_per_hr:.9f} rad/hr")
else:
    print("No data point found at 3 inches.")