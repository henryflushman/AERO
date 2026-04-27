import Orbits_Functions as of
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta


import os
os.system('cls')

TLEs = [
    [
        "DELTA 2 DEB",
        "1 25138U 98002E   25321.77251655  .00000756  00000-0  47462-3 0  9993",
        "2 25138  26.8788 345.9874 0323619 168.6975 354.1864 13.72117421392909"
    ],
    [
        "DELTA 1 DEB",
        "1 10227U 77065BE  25321.71492541  .00006135  00000-0  20445-2 0  9995",
        "2 10227  28.9572 330.4974 0675997 286.1808  66.5427 13.52266693346021"
    ],
    [
        "ATLAS 5 CENTAUR DEB",
        "1 46656U 18079CJ  25321.49577821 -.00000102  00000-0  00000-0 0  9992",
        "2 46656  13.3737  20.3709 4531874 347.7674   4.2683  1.96953546 36789"
    ],
    [
        "TITAN 3C TRANSTAGE DEB",
        "1 38695U 68081V   25321.21743478 -.00000041  00000-0  00000+0 0  9998",
        "2 38695   1.8787  19.2505 0430881 271.1786  96.3341  0.99387212 61242"
    ]
]

"""
_, ECI1 = of.parse_tle(TLEs[0])
_, ECI2 = of.parse_tle(TLEs[1])

of.lamberts_porkchop(TLEs[0], TLEs[2], 6*3600, 12*3600, 300)
"""

def epoch_to_datetime(epoch_year, epoch_day):
    # Coerce to Python types
    y = int(float(epoch_year))
    d = float(epoch_day)

    # TLE epochs may use 2-digit years (00–56 => 2000–2056; 57–99 => 1957–1999)
    if y < 100:
        y = 2000 + y if y < 57 else 1900 + y

    day_int = int(d)                  # 1..365/366
    frac = d - day_int                # fractional day
    base = datetime(y, 1, 1) + timedelta(days=day_int - 1)
    return base + timedelta(seconds=frac * 86400.0)

fps = 30
secs = 12
tspan = [0.0, 24*3600.0]
frames = fps*secs
t_eval = np.linspace(tspan[0], tspan[1], frames)

solutions = []
labels = []
COEs_all = []
ECI_all = []
epochs = []

for tle in TLEs:
    _, _, time = of.parse_tle(tle)
    epoch_year, epoch_day = time
    epochs.append(epoch_to_datetime(epoch_year, epoch_day))

time_initial = max(epochs)

# Orbit Plotter
for tle, time_now in zip(TLEs, epochs):
    COEs, ECI, _ = of.parse_tle(tle)
    ECI_all.append(ECI)
    rx, ry, rz, vx, vy, vz = ECI
    r0 = np.array([rx, ry, rz])
    v0 = np.array([vx, vy, vz])
    
    ds = (time_initial - time_now).total_seconds()
    if ds != 0.0:
        sol_new = of.ODEprimer(r0, v0, [0.0, ds])
        rx, ry, rz, vx, vy, vz = sol_new.y[0:6, -1]
        r0 = np.array([rx, ry, rz])
        v0 = np.array([vx, vy, vz])
    
    y0 = np.hstack((r0, v0))
    COEs = of.ECI2COEs(y0)
    sol = of.ODEprimer(r0, v0, tspan, t_eval)
    solutions.append(sol)
    COEs_all.append(COEs)
    labels.append(tle[0].strip())

y0 = [sol.y for sol in solutions]

of.plot_orbit(y0, labels, show=True)
#of.animate_orbits(solutions, labels, save="Satellite_Trails.mp4",
#                  framerate=fps, max_trail=15)
