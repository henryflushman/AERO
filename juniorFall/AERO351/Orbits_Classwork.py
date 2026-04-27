""" In class example """

import numpy as np
import matplotlib.pyplot as plt
import Orbits_Functions as of
import os
os.system('cls')

ecc = 1.2
alt = 200
inc = np.deg2rad(50)
RAAN = np.deg2rad(75)
argp = np.deg2rad(80)
TA = np.deg2rad(80)
mu = 398600
h = np.sqrt(mu*(1+ecc)*(alt+6578))

r_eci, v_eci = of.COEs2ECI(alt, ecc, inc, RAAN, argp, TA, mu)

print(r_eci)
print(v_eci)

delv = of.phasechange_delV(13000, 11000, .5) # type: ignore
print(delv)