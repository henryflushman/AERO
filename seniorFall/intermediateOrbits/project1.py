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
# ║   Assignment  :  Project 1                                  ║
# ║   Date        :  September 17, 2026                         ║
# ╚═════════════════════════════════════════════════════════════╝


import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from IntermediateOrbitsTools import ClassicalOrbitalElements, R_EARTH_KM
from IntermediateOrbitsTools.RelativeMotion import (
    cwhill_transition_matrix,
    cwhill_delr_delv,
)


# === CONFIG =======================

TARGET_ECI_0  = np.array(
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0]
)   # Replace with the actual ECI coordinates of the target

A = ClassicalOrbitalElements(
    [51400, 0.0006387, 51.65, 15, 157, 15],
    names=names,
    angle_unit="deg",
)


CHASER_LVLH_0 = np.array(
    [0.0, 100.0, 0.0],
    [0.0, 0.0, 0.0]
) # 100.0 km away from target

RENDEZVOUS_1_PARAMS = np.array(
    [40.0, 20.0]
)

# ==================================


# Maneuver 1 will be a radially inward burn.

# a = 100km initially, therefore delv = (100.0/2.0)*mean_motion

delv_req_for_current = (CHASER_LVLH_0[0, 1]/2.0)*A.mean_motion_rad_s

delv_req_for_desired = (RENDEZVOUS_1_PARAMS[0]/2.0)*A.mean_motion_rad_s

delv_req_for_change = delv_req_for_desired - delv_req_for_current

# Now determine where trajectories intersect

# At next intersection
delv_req_for_same_traj = -delv_req_for_desired

# Now that trajectories have been aligned
delv_hop1 = # phasing maneuver large enough to reach 1km, try to do in multiple revs if possible

# Now sc are 1 km away from eachother

# Hop again to 300m same method as delv_hop1
delv_hop2 = # phasing maneuver large enough to reach 300m, try to do in multiple revs if possible

# Now sc are 300 m away from eachother

# hop approach to 20m with 
delv_hop3 = # phasic maneuver large enough to reach 20m, try to do in multiple revs if possible

# Begin vbar approach with cm/s relative speed
delv_vbar = # vbar approach to 0 km relative distance


solve_ivp(
    fun=lambda t, y: cwhill_delr_delv(t, y, A.mean_motion_rad_s),
    t_span=(0, 3600),  # 1 hour for example
    y0=CHASER_LVLH_0.flatten(),
    method="RK45",
)