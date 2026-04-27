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

_, ECI2, _ = of.parse_tle(TLEs[2])
_, ECI1, _ = of.parse_tle(TLEs[1])
of.lamberts_porkchop(ECI2, ECI1, 2*3600, 6*3600, 150)