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
# ║   Assignment  :  Homework 2                                 ║
# ║   Date        :  September 6, 2026                          ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ===
import warnings

import matplotlib.pyplot as plt
import numpy as np

from IntermediateOrbitsTools import ClassicalOrbitalElements, R_EARTH_KM
from IntermediateOrbitsTools.OrbitalFrames import (
    eci_to_lvlh_dcm,
    relative_state_eci_to_lvlh,
    relative_state_lvlh_to_eci,
)
from IntermediateOrbitsTools.RelativeMotion import (
    cwhill_delr_delv,
    cwhill_propagate,
    exact_relative_motion,
    relative_separation_norm,
)
from IntermediateOrbitsTools.Propagators import propagate_two_body


# ============================================================
# Problem 3: Relative motion after 10 orbital periods
# ============================================================


# Step 1: define chief orbit using perigee altitude and eccentricity

perigee_altitude_km = 250.0
eccentricity = 0.1
perigee_radius_km = R_EARTH_KM + perigee_altitude_km
semi_major_axis_km = perigee_radius_km / (1.0 - eccentricity)

chief = ClassicalOrbitalElements(
    [semi_major_axis_km, eccentricity, 51.0, 0.0, 0.0, 0.0],
    names=["a", "e", "i", "raan", "argp", "nu"],
    angle_unit="deg",
)

# The CW equations are derived for a near-circular chief orbit. This scenario uses
# e = 0.1, so the result should be interpreted as an approximate comparison rather
# than a strict CW-valid solution.
if abs(eccentricity) > 0.05:
    warnings.warn(
        "CW is being used with a chief eccentricity of 0.1; this is an approximation and should be compared to the exact relative-motion solution.",
        RuntimeWarning,
    )

# chief state to ECI
chief_state = chief.to_state()

# orbital period and mean motion
mean_motion_rad_s = chief.mean_motion_rad_s
period_s = chief.period_s
final_time_s = 10.0 * period_s


# Step 2: define the relative state in the chief LVLH frame
delta_r0 = np.array([-1.0, -1.0, 0.0], dtype=float)
delta_v0 = np.array([0.0, 0.002, 0.0], dtype=float)


# Step 3: propagate the relative state using the CW equations
cw_r, cw_v = cwhill_propagate(
    delta_r0,
    delta_v0,
    final_time_s,
    mean_motion_rad_s,
    eccentricity=eccentricity,
)
phi_r, phi_v = cwhill_delr_delv(
    delta_r0,
    delta_v0,
    final_time_s,
    mean_motion_rad_s,
    eccentricity=eccentricity,
)

cw_separation = relative_separation_norm(cw_r)
phi_separation = relative_separation_norm(phi_r)

print("\nCW propagation after 10 periods:")
print(f"delta_r = {cw_r} km")
print(f"delta_v = {cw_v} km/s")
print(f"Relative separation magnitude = {cw_separation:.6f} km")

print("\nCW matrix-form propagation after 10 periods:")
print(f"delta_r = {phi_r} km")
print(f"delta_v = {phi_v} km/s")
print(f"Relative separation magnitude = {phi_separation:.6f} km")

# DCM from ECI to LVLH for the chief at t = 0.
C_eci_to_lvlh = eci_to_lvlh_dcm(chief_state.r_km, chief_state.v_km_s)

deputy_r_eci, deputy_v_eci = relative_state_lvlh_to_eci(
    chief_state.r_km,
    chief_state.v_km_s,
    delta_r0,
    delta_v0,
)

rho_exact, rho_dot_exact = exact_relative_motion(
    chief_state.r_km,
    chief_state.v_km_s,
    deputy_r_eci,
    deputy_v_eci,
)

exact_separation = relative_separation_norm(rho_exact)

print("\nExact relative state at the initial instant (snapshot):")
print(f"rho_exact = {rho_exact} km")
print(f"rho_dot_exact = {rho_dot_exact} km/s")
print(f"Exact relative separation magnitude = {exact_separation:.6f} km")


# Step 5: create a time history for plotting

num_points = 500
times_s = np.linspace(0.0, final_time_s, num_points)
relative_positions = np.zeros((num_points, 3))
relative_velocities = np.zeros((num_points, 3))
separations_km = np.zeros(num_points)
exact_relative_positions = np.zeros((num_points, 3))
exact_relative_velocities = np.zeros((num_points, 3))
exact_separations_km = np.zeros(num_points)

# Propagate both spacecraft in ECI so the exact result can be compared with CW.
chief_history = propagate_two_body(
    chief_state,
    (0.0, final_time_s),
    times_s=times_s,
)
deputy_history = propagate_two_body(
    [*deputy_r_eci, *deputy_v_eci],
    (0.0, final_time_s),
    times_s=times_s,
)

for i, t in enumerate(times_s):
    r_t, v_t = cwhill_propagate(
        delta_r0,
        delta_v0,
        t,
        mean_motion_rad_s,
        eccentricity=eccentricity,
    )
    relative_positions[i, :] = r_t
    relative_velocities[i, :] = v_t
    separations_km[i] = relative_separation_norm(r_t)

    rho_t, rho_dot_t = relative_state_eci_to_lvlh(
        chief_history.positions_km[i],
        chief_history.velocities_km_s[i],
        deputy_history.positions_km[i],
        deputy_history.velocities_km_s[i],
    )
    exact_relative_positions[i, :] = rho_t
    exact_relative_velocities[i, :] = rho_dot_t
    exact_separations_km[i] = relative_separation_norm(rho_t)

print("\nExact-vs-CW comparison after 10 periods:")
print(f"Exact relative position = {exact_relative_positions[-1]} km")
print(f"Exact relative velocity = {exact_relative_velocities[-1]} km/s")
print(f"Exact relative separation = {exact_separations_km[-1]:.6f} km")
print(f"Position-model difference = {np.linalg.norm(cw_r - exact_relative_positions[-1]):.6f} km")


# Step 6: graph the relative motion results

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

ax.plot(
    relative_positions[:, 0],
    relative_positions[:, 1],
    relative_positions[:, 2],
    linewidth=2,
    label="CW relative motion",
)
ax.plot(
    exact_relative_positions[:, 0],
    exact_relative_positions[:, 1],
    exact_relative_positions[:, 2],
    linewidth=1.5,
    linestyle="--",
    label="Exact two-body relative motion",
)

ax.plot(
    [0.0], [0.0], [0.0],
    marker="o",
    linestyle="None",
    color="black",
    markersize=6,
    label="Chief",
)
ax.plot(
    [relative_positions[0, 0]],
    [relative_positions[0, 1]],
    [relative_positions[0, 2]],
    marker="o",
    linestyle="None",
    color="tab:green",
    markersize=7,
    label="Initial deputy position",
)
ax.plot(
    [relative_positions[-1, 0]],
    [relative_positions[-1, 1]],
    [relative_positions[-1, 2]],
    marker="o",
    linestyle="None",
    color="tab:red",
    markersize=7,
    label="10 periods later",
)

ax.set_xlabel("Radial (R) [km]")
ax.set_ylabel("Along-track (T) [km]")
ax.set_zlabel("Cross-track (N) [km]")
ax.set_title("Relative Motion in Chief LVLH Frame")
ax.legend()
plt.tight_layout()
plt.show()


plt.figure(figsize=(10, 5))
plt.plot(times_s / period_s, separations_km, label="CW separation")
plt.plot(
    times_s / period_s,
    exact_separations_km,
    linestyle="--",
    label="Exact two-body separation",
)
plt.axvline(10.0, linestyle="--", color="k", label="10 periods")
plt.xlabel("Time [orbital periods]")
plt.ylabel("Relative separation [km]")
plt.title("Relative Separation vs Time")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()



# ============================================================
# Problem 4: del_r(t=15 minutes) given chief period and relative position and velocity
# ============================================================

T = 90.0 * 60.0
t = 15.0 * 60.0
n = 2.0 * np.pi / T

delr_0 = np.array([1.0, 0.0, 0.0])
delv_0 = np.array([0.0, 0.010, 0.0])

delr_3ans, delv_3ans = cwhill_propagate(
    delr_0,
    delv_0,
    t,
    n,
)

separation_3ans = relative_separation_norm(
    delr_3ans
)

print("\nTextbook Problem 7.7:")
print(f"Mean motion = {n:.6f} rad/s")
print(f"Relative position after 15 minutes = {delr_3ans}")
print(f"Relative velocity after 15 minutes = {delv_3ans}")
print(f"Relative separation after 15 minutes = {separation_3ans}")