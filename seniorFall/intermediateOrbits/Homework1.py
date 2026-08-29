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
# ║   Course      :  AERO4452 - Intermediate Orbits             ║
# ║                                                             ║
# ║   Assignment  :  Homework 1                                 ║
# ║   Date        :  August 26, 2026                            ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ===
import matplotlib.pyplot as plt
import numpy as np

from IntermediateOrbitsTools import (
    ClassicalOrbitalElements,
    StateVector,
    eci_to_lvlh,
    propagate_two_body,
    R_EARTH_KM
)


names = ["h", "e", "i", "raan", "argp", "theta"]
A = ClassicalOrbitalElements(
    [51400, 0.0006387, 51.65, 15, 157, 15],
    names=names,
    angle_unit="deg",
)
B = ClassicalOrbitalElements(
    [51398, 0.0072696, 50, 15, 140, 15],
    names=names,
    angle_unit="deg",
)

# Propagate both spacecraft for 10 periods of A.
t_final = 10.0 * A.period_s
prop_A = propagate_two_body(A.to_state(), (0.0, t_final), method="universal_anomaly", time_step_s=30.0)
prop_B = propagate_two_body(B.to_state(), (0.0, t_final), method="universal_anomaly", time_step_s=30.0)

t = prop_A.times_s

# B relative to A, expressed in A's LVLH/RTN frame.
rho = np.array([
    eci_to_lvlh(r_B - r_A, r_A, v_A)
    for r_A, v_A, r_B in zip(
        prop_A.positions_km,
        prop_A.velocities_km_s,
        prop_B.positions_km,
    )
])

distance = np.linalg.norm(rho, axis=1)
k = np.argmin(distance)

print(f"A period: {A.period_min:.3f} min")
print(f"Closest approach: {distance[k]:.3f} km")
print(f"Time of closest approach: {t[k]:.3f} s ({t[k] / 3600.0:.3f} hr)")
print(f"Relative position at closest approach [R, T, N]: {rho[k]} km")

# One required graph: relative position history with closest-approach time marked.
fig = plt.figure()

ax = fig.add_subplot(
    111,
    projection="3d",
)

ax.plot(
    rho[:, 0],
    rho[:, 1],
    rho[:, 2],
    linewidth=0.5,
    label="Spacecraft B"
)

# sc A
ax.scatter(
    0,
    0,
    0,
    label="Spacecraft A"
)

# sc B
ax.scatter(
    rho[k, 0],
    rho[k, 1],
    rho[k, 2],
    label="Closest Approach"
)

ax.set_xlabel("Radial (R) [km]")
ax.set_ylabel("Along-track (T) [km]")
ax.set_zlabel("Orbit-Normal (N) [km]")

ax.set_title(
    "Spacecraft B Relative Position in Spacecraft A LVLH Frame"
)

ax.legend()

plt.tight_layout()
plt.show()

plt.figure()
plt.plot(t / A.period_s, distance[:], 'k')
plt.axvline(t[k] / A.period_s, linestyle="--", label="Closest Approach")
plt.xlabel("Time [Periods of Spacecraft A]")
plt.ylabel("Distance Between Satellites [km]")
plt.title("Spacecraft B Distance to Spacecraft A")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# Problem 2

TARGET_POS = np.array([0, 6678, 0])    # ECI
TARGET_VEL = np.array([0, 0, 7.7258])  # ECI
TARGET_ACC = np.array([0, -0.0089, 0]) # ECI

CHASER_POS = np.array([0, 0, 6628])    # ECI
CHASER_VEL = np.array([0, -7.7549, 0]) # ECI
CHASER_ACC = np.array([0, 0, -0.0091]) # ECI

# Build LVLH DCM
R_lvlh = TARGET_POS / np.linalg.norm(TARGET_POS)
N_lvlh = np.cross(TARGET_POS, TARGET_VEL)/(np.linalg.norm(np.cross(TARGET_POS, TARGET_VEL)))
T_lvlh = np.cross(N_lvlh, R_lvlh)

C_lvlh_eci = np.vstack([R_lvlh, T_lvlh, N_lvlh])


print(C_lvlh_eci)

dr_0 = CHASER_POS - TARGET_POS

rho = C_lvlh_eci @ dr_0


print(rho)