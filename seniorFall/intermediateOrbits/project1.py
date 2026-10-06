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

CHASER_LVLH_0 = np.array(
    [0.0, 100.0, 0.0],
    [0.0, 0.0, 0.0]
) # 100.0 km away from target

# ==================================


# Maneuver 1 will be a radially inward burn.

# a = 100km initially, therefore delv = (100.0/2.0)*mean_motion

delv_1 = np.array([-1.0, 0.0, 0.0])