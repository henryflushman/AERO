""" LAB 2: SSWT Flow Vis """

# === Imports =============================
import numpy as np
import matplotlib.pyplot as plt
import os
# =========================================


# === Load Data .txt File =================
# Variable constructor helper function
def VarsFromData(dataName):
    data = np.genfromtxt(
        dataName,
        delimiter='\t',
        skip_header=1
    )
    t = data[:, 0]
    pitotPress = data[:, 1]
    plenumPress = data[:, 2]
    tankPress = data[:, 3]
    staticPress = data[:, 4]
    
    return t, pitotPress, plenumPress, tankPress, staticPress

# Load and save data to assigned variables
t1, pitotPress1, plenumPress1, tankPress1, staticPress1 = VarsFromData('Lab2Data_run1.txt')
t2, pitotPress2, plenumPress2, tankPress2, staticPress2 = VarsFromData('Lab2Data_run2.txt')
t3, pitotPress3, plenumPress3, tankPress3, staticPress3 = VarsFromData('Lab2Data_run3.txt')
# =========================================


# === Analysis ============================
ax1 = plt.subplot(2, 3, 1)
ax1.plot(t1, staticPress1)
ax1.axvline(11.5, linestyle='--', color='black')
ax1.axvline(14.5, linestyle='--', color='black')
ax1.axvspan(11.5, 14.5, color='green', alpha=0.3)
ax1.set_ylabel('Static Pressure [psig]')
ax1.set_title('Run 1')
ax1.grid(True)

ax2 = plt.subplot(2, 3, 2)
ax2.plot(t2, staticPress2)
ax2.axvline(9, linestyle='--', color='black')
ax2.axvline(12, linestyle='--', color='black')
ax2.axvspan(9, 12, color='green', alpha=0.3)
ax2.set_title('Run 2')
ax2.grid(True)

ax3 = plt.subplot(2, 3, 3)
ax3.plot(t3, staticPress3)
ax3.axvline(8.5, linestyle='--', color='black')
ax3.axvline(11, linestyle='--', color='black')
ax3.axvspan(8.5, 11, color='green', alpha=0.3)
ax3.set_title('Run 3')
ax3.grid(True)

# ----- Plenum pressure row -----
ax4 = plt.subplot(2, 3, 4)
ax4.plot(t1, plenumPress1)
ax4.axvline(11.5, linestyle='--', color='black')
ax4.axvline(14.5, linestyle='--', color='black')
ax4.axvspan(11.5, 14.5, color='green', alpha=0.3)
ax4.set_ylabel('Plenum Pressure [psig]')
ax4.set_xlabel('Time [s]')
ax4.grid(True)

ax5 = plt.subplot(2, 3, 5)
ax5.plot(t2, plenumPress2)
ax5.axvline(9, linestyle='--', color='black')
ax5.axvline(12, linestyle='--', color='black')
ax5.axvspan(9, 12, color='green', alpha=0.3)
ax5.set_xlabel('Time [s]')
ax5.grid(True)

ax6 = plt.subplot(2, 3, 6)
ax6.plot(t3, plenumPress3)
ax6.axvline(8.5, linestyle='--', color='black')
ax6.axvline(11, linestyle='--', color='black')
ax6.axvspan(8.5, 11, color='green', alpha=0.3)
ax6.set_xlabel('Time [s]')
ax6.grid(True)

plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

steadyStateRange = np.array([
    [11.5, 14.5],
    [9, 12],
    [8.5, 11]
])

tRuns = [t1, t2, t3]
plenumRuns = [plenumPress1, plenumPress2, plenumPress3]
staticRuns = [staticPress1, staticPress2, staticPress3]

# Helper function
def MachFromPressRatio(P0_P, gamma=1.4):
    P0_P = np.asarray(P0_P, dtype=float)
    return np.sqrt((P0_P**((gamma-1)/gamma)-1.0) * (2.0/(gamma-1)))

# Steady State Samples
P0_all_psig = []
Ps_all_psig = []
ratioAll = []

for i, (tStart, tEnd) in enumerate(steadyStateRange):
    t = tRuns[i]
    P0_abs = plenumRuns[i] + 14.7
    Ps_abs = staticRuns[i] + 14.7

    mask = (t >= tStart) & (t <= tEnd)
    if not np.any(mask):
        raise ValueError(f"No points in steady window for run {i+1}: [{tStart}, {tEnd}]")

    R = P0_abs[mask]/Ps_abs[mask]
    R = R[np.isfinite(R)]
    if R.size == 0:
        raise ValueError(f"All ratios invalied in ss window")

    # Append all steady samples
    ratioAll.append(R)

ratioAll = np.concatenate(ratioAll)

R_mean = np.mean(ratioAll)
R_max = np.max(ratioAll)
R_min = np.min(ratioAll)

M_mean = MachFromPressRatio(R_mean)
M_max = MachFromPressRatio(R_max)
M_min = MachFromPressRatio(R_min)

mu_mean = np.degrees(np.arcsin(1.0 / M_mean))
mu_max  = np.degrees(np.arcsin(1.0 / M_max))
mu_min  = np.degrees(np.arcsin(1.0 / M_min))

print(f"Ratio stats over all SS windows: mean={R_mean:.4f}, min={R_min:.4f}, max={R_max:.4f}")
print(f"Mach stats: M_mean={M_mean:.4f}, M_min={M_min:.4f}, M_max={M_max:.4f}")
print(f"Mach angle stats (deg): mu_mean={mu_mean:.3f}, mu_min={mu_min:.3f}, mu_max={mu_max:.3f}")

# --- Theta-beta-Mach Plot -----------
A_Astar = 4.8
# Using table A.1
machNumber1_theory = 3.13
machAngle1_theory = np.arcsin(1/machNumber1_theory)
print("\nTheoretical mach angle (from A/A*):", np.rad2deg(machAngle1_theory))
# Using oblique shock graph @M=3.2
betaMax = 66        # degrees
thetaMax = 35.5     # degrees

# Theta-Beta-Mach Graph
def ThetaFromBetaMach(beta):
    Mn1_sq = (machNumber1_theory*np.sin(beta))**2
    numerator = 2*(1/np.tan(beta))*(Mn1_sq-1)
    denominator = (machNumber1_theory**2)*(1.4+np.cos(2*beta)) + 2
    theta = np.arctan(numerator/denominator)
    return theta

beta_min = np.arcsin(1/machNumber1_theory) + np.deg2rad(0.1)
beta_max = np.deg2rad(89.9)

beta_rad = np.linspace(beta_min, beta_max, 2000)
theta_rad = ThetaFromBetaMach(beta_rad)

valid = theta_rad > 0
beta_deg = np.rad2deg(beta_rad[valid])
theta_deg = np.rad2deg(theta_rad[valid])

# Experimental Beta Value
betaExperimental = np.array([
    90, 82.53, 69.04, 64.15, 58.86, 71.26, 48.76, 40.29, 35.09, 28.32
])
betaRad_experimental = np.deg2rad(betaExperimental)
thetaRad_experimental = ThetaFromBetaMach(betaRad_experimental)
thetaDeg_experimental = np.rad2deg(thetaRad_experimental)

plt.figure(figsize=(8, 5))
plt.plot(theta_deg, beta_deg)
plt.plot(thetaDeg_experimental, betaExperimental, 'o')
plt.xlabel("Deflection angle theta [deg]")
plt.ylabel("Shock angle beta [deg]")
plt.title(f"θ–β–M Relation (M = {machNumber1_theory:.2f}, γ = 1.4)")
plt.grid(True)
plt.tight_layout()
plt.legend(['Theoretical', 'Experimental'])
plt.show()

# =========================================