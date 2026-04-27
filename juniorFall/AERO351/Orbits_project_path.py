import Orbits_Functions as of
import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timedelta

import os

os.system('cls')

def stitch_segments(segments):
    out = segments[0]
    for seg in segments[1:]:
        if out.shape[1] > 0 and seg.shape[1] > 0 and np.allclose(out[:, -1], seg[:, 0]):
            out = np.concatenate((out, seg[:, 1:]), axis=1)
        else:
            out = np.concatenate((out, seg), axis=1)
    return out

def propagate_segments(y_segments, tspan, teval=None):
    """
    Propagate a list of orbit state segments forward by a given time span.

    Parameters
    ----------
    y_segments : list of np.ndarray
        Each element is a 6×N array of [r; v] states (output from ODEprimer.y[0:6]).
    tspan : list or tuple
        [t0, tf] time span (seconds) to propagate each orbit.
    teval : array-like or None
        Optional time evaluation points for ODEprimer.

    Returns
    -------
    new_segments : list of np.ndarray
        List of propagated 6×N arrays.
    new_COEs : list
        List of COEs (from of.ECI2COEs) at the final propagated state.
    """
    new_segments = []
    new_COEs = []

    for y in y_segments:
        # initial state for this segment
        r0 = y[0:3, 0]
        v0 = y[3:6, 0]

        # propagate
        sol = of.ODEprimer(r0, v0, tspan, teval=teval)

        # final state
        rx, ry, rz, vx, vy, vz = sol.y[0:6, -1]
        r_end = np.array([rx, ry, rz])
        v_end = np.array([vx, vy, vz])
        state_end = np.hstack((r_end, v_end))

        # store results
        new_segments.append(sol.y[0:6])
        new_COEs.append(of.ECI2COEs(state_end))

    return new_segments, new_COEs


""" MUST STAY IN SYNC WITH DEBRIS FOR 5 PERIODS """

TLEs_original = [
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


# --------------- PROPAGATE ORBITS TO COMMON EPOCH ---------------
def epoch_to_datetime(epoch_year, epoch_day):
    # Coerce to Python types
    y = int(float(epoch_year))
    d = float(epoch_day)

    if y < 100:
        y = 2000 + y if y < 57 else 1900 + y

    day_int = int(d)
    frac = d - day_int 
    base = datetime(y, 1, 1) + timedelta(days=day_int - 1)
    return base + timedelta(seconds=frac * 86400.0)

y0_t0_2_trans1 = []
labels = []
COEs_t0 = []
ECI_all = []
epochs = []
COEs_5periods = []

for tle in TLEs_original:
    _, _, time = of.parse_tle(tle)
    epoch_year, epoch_day = time
    epochs.append(epoch_to_datetime(epoch_year, epoch_day))

time_initial = max(epochs)
tspan_2t0 = [0.0, 24.1488333*5*3600.0] # 5 periods of first debris
t_eval_2t0 = None           # ADD IN TVAL AT THE END

# Orbit Plotter
for tle, time_now in zip(TLEs_original, epochs):
    COEs, _, _ = of.parse_tle(tle)
    if COEs[4] < 1e-3:
        COEs[4] = 0
    r0, v0 = of.COEs2ECI(COEs[0], COEs[1], COEs[2], COEs[3], COEs[4], COEs[5])
    
    ds = (time_initial - time_now).total_seconds()
    if ds != 0.0:
        sol_new = of.ODEprimer(r0, v0, [0.0, ds])
        
        rx, ry, rz, vx, vy, vz = sol_new.y[0:6, -1]
        r0 = np.array([rx, ry, rz])
        v0 = np.array([vx, vy, vz])
        
    y0 = np.hstack((r0, v0))
    COEs = of.ECI2COEs(y0)
    sol = of.ODEprimer(r0, v0, tspan_2t0, teval=t_eval_2t0)
    rx, ry, rz, vx, vy, vz = sol.y[0:6, -1]
    r0 = np.array([rx, ry, rz])
    v0 = np.array([vx, vy, vz])
    y0_5periods = np.hstack((r0, v0))
    COEs_5periods.append(of.ECI2COEs(y0_5periods))
    y0_t0_2_trans1.append(sol.y[0:6])
    COEs_t0.append(COEs)
    labels.append(tle[0].strip())
    
labels.append("SPACECRAFT")

# --------------- COAST UNTIL FIRST MANEUVER ---------------

sc_COEs = COEs_t0[3]
sc_TA = sc_COEs[5]
sc_ecc = np.linalg.norm(sc_COEs[1])
sc_h = np.linalg.norm(sc_COEs[0])

print(np.rad2deg(sc_TA))

dt_depart = ((sc_TA - np.deg2rad(89.67))/(2*np.pi))*24.1488333

print(dt_depart)

tspan_2depart1 = [0.0, dt_depart]
t_eval_2depart1 = None           # ADD IN TVAL AT THE END

y0_2depart1 = []
# Reset COEs
COEs_2depart1 = []

y0_2depart1, COEs_2depart1 = propagate_segments(
    y0_t0_2_trans1, tspan_2depart1, teval=t_eval_2depart1
)


y0_full = [stitch_segments([seg1, seg2]) 
           for seg1, seg2 in zip(y0_t0_2_trans1, y0_2depart1)]

p = 34217.027
a = 37858.048
e = (a-p)/(a+p)
print(e)