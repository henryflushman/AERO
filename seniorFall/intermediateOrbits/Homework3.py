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
# ║   Assignment  :  Homework 3                                 ║
# ║   Date        :  September 12, 2026                         ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ===
import warnings

import matplotlib.pyplot as plt
import numpy as np

from IntermediateOrbitsTools import ClassicalOrbitalElements, R_EARTH_KM
from IntermediateOrbitsTools.RelativeMotion import (
    cwhill_transition_matrix,
    cwhill_delr_delv,
)



# ============================================================
# Problem 5: CW Relative Speed Propagation
# ============================================================

# Initial relative state at t = 0
delr_0 = np.array([0.0, 0.0, 0.0], dtype=float)     # m
delv_0 = np.array([1.0, -1.0, 1.0], dtype=float)/1000.0    # km/s

n_rad_s = 1.0
period_s = 2.0 * np.pi / n_rad_s
t_elapsed_s = period_s / 4.0

# Propagate relative state using the Clohessy-Wiltshire state transition
delr_quarter, delv_quarter = cwhill_delr_delv(
    delr_0,
    delv_0,
    t_elapsed_s,
    n_rad_s,
)

# Relative speed (norm of relative velocity vector)
relative_speed = np.linalg.norm(delv_quarter)

print("Problem 5:")
print(f"  Initial relative velocity:  {delv_0} km/s ({delv_0*1000.0} m/s)")
print(f"  Relative velocity at T/4:   {delv_quarter} km/s ({delv_quarter*1000.0} m/s)")
print(f"  Relative speed at T/4:      {relative_speed:.6f} km/s ({relative_speed*1000.0:.6f} m/s)\n")


# ============================================================
# Problem 6: Coplanar CW Relative Motion
# ============================================================

delr_0 = np.array([1.0, np.pi, 0.0], dtype=float)/1000.0    # km
delv_0 = np.array([np.pi/16.0, 7.0/4.0, 0.0], dtype=float)/1000.0    # km/s

n_rad_s = 1.0
period_s = 2.0 * np.pi / n_rad_s
t_elapsed_s = np.pi / (2.0 * n_rad_s)

delr, delv = cwhill_delr_delv(
    delr_0,
    delv_0,
    t_elapsed_s,
    n_rad_s,
)
print("Problem 6:")
print(f"  Initial relative position:  {delr_0} km")
print(f"  Initial relative velocity:  {delv_0} km/s")
print(f"  Relative position at t:    {delr} km")
print(f"  Relative velocity at t:    {delv} km/s\n")

# ============================================================
# Problem 7: Two-Impulse Rendezvous
# ============================================================

target = ClassicalOrbitalElements(
    a=(300 + R_EARTH_KM),  # km
    e=0.0,
    i=0.0,
    raan=0.0,
    argp=0.0,
    nu=0.0,
)

n_rad_s = target.mean_motion_rad_s

t_f = 30.0 * 60.0

delr_0 = np.array([-1.0, 0.0, 0.0], dtype=float)  # km

phi_rr_tf, phi_rv_tf, phi_vr_tf, phi_vv_tf = cwhill_transition_matrix(
    t_f,
    n_rad_s
)

delv_0_before = np.array([0.0, 0.0, 0.0], dtype=float)  # km/s
delv_0_after = np.linalg.inv(phi_rv_tf) @ (-phi_rr_tf) @ delr_0  # km/s
delv_f_before = phi_vr_tf @ delr_0 + phi_vv_tf @ delv_0_after  # km/s
delv_f_after = np.array([0.0, 0.0, 0.0], dtype=float)

delv_0 = np.linalg.norm(delv_0_before - delv_0_after)  # km/s
delv_f = np.linalg.norm(delv_f_before - delv_f_after)  # km/s

delv_total = delv_0 + delv_f  # km/s

print("Problem 7:")
print(f"  Initial relative velocity before impulse:  {delv_0_before} km/s")
print(f"  Initial relative velocity after impulse:   {delv_0_after} km/s ({delv_0_after*1000.0} m/s)")
print(f"  Final relative velocity before impulse:    {delv_f_before} km/s ({delv_f_before*1000.0} m/s)")
print(f"  Final relative velocity after impulse:     {delv_f_after} km/s")
print(f"  Delta-v at start:                          {delv_0:.6f} km/s ({delv_0*1000.0:.4f} m/s)")
print(f"  Delta-v at end:                            {delv_f:.6f} km/s ({delv_f*1000.0:.4f} m/s)")
print(f"  Total delta-v:                             {delv_total:.6f} km/s ({delv_total*1000.0:.4f} m/s)\n")