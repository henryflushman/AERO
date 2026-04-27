import os

from matplotlib.pylab import f
os.system('cls')

from matplotlib import pyplot as plt
from numpy import sin, cos, pi
import numpy as np
import datetime
import Orbits_Functions as of


""" Curtis 2.23 """
print("Curtis Exercise 2.23:\n")
# Constants
mu_earth = 398600 #km^3/s^2                   grav parameter
r_earth = 6378 #km                            radius of earth
r_p = r_earth + 500 #km                       radius of periapsis
v_p = 10 #km/s                                velocity at periapsis
nu = np.deg2rad(120)
h = r_p * v_p #km^2/s                         specific angular momentum
a = ((2/r_p) - (v_p**2/mu_earth))**-1 #km     semi-major axis
ecc = np.arccos((h**2/(mu_earth*r_p)) - 1) #eccentricity
flightpath = np.arctan((ecc*np.sin(nu)) / (1 + ecc*np.cos(nu))) #rad   flightpath angle
r = of.trueanom2pos(nu, ecc, h)
print(f"The altitude of the satellite: {r-r_earth} km")
print(f"Flightpath angle: {np.rad2deg(flightpath)} degrees\n\n\n")

""" Curtis 2.36 """
print("\nCurtis Exercise 2.36:\n")
r_p = r_earth + 250 #km
v_p = 11 #km/s
nu = np.deg2rad(100)
h = r_p * v_p #km^2/s
ecc = (v_p*h) / mu_earth - 1
a = (h**2 / mu_earth) * (1 / (1 - ecc**2))
r = of.trueanom2pos(nu, ecc, h)
exc_speed = np.sqrt(mu_earth/np.abs(a))
print(f"Eccentricity: {ecc}")
print(f"Excess speed: {exc_speed} km/s")
print(f"Radial position of the satellite at 100 degrees: {r} km")
v_r = (mu_earth / h) * (1+ecc*np.cos(nu))
v_az = (mu_earth / h) * (ecc*np.sin(nu))
v_mag = np.sqrt(v_r**2 + v_az**2)
print(f"Azimuthal velocity at 100 degrees: {v_az} km/s\nRadial velocity at 100 degrees: {v_r} km/s\nVelocity magnitude at 100 degrees: {v_mag} km/s\n\n\n")

""" Curtis 2.37 """
print("\nCurtis Exercise 2.37:\n")
r_obsv = 402000 #km
v_obsv = 2.23 #km/s
nu = np.deg2rad(150)
eps = (v_obsv**2)/2 - mu_earth/r_obsv
a = -mu_earth/(2*eps)
ecc = (-(r_obsv/a)*np.cos(nu) + np.sqrt(((r_obsv/a)*np.cos(nu))**2 - 4*((r_obsv/a) - 1)))/2
print(f"Eccentricity: {ecc}")
r_p = a*(1 - ecc)
print(f"Altitude of closest approach: {r_p-r_earth} km")
v_p = np.sqrt(mu_earth*((2/r_p) - (a**-1)))
print(f"Velocity at closest approach: {v_p} km/s\n\n\n")

""" Curtis 3.8 """
print("\nCurtis Exercise 3.8:\n")
r_p = 200 + r_earth
r_a = 600 + r_earth
ecc = (r_a - r_p) / (r_a + r_p)
a = (r_a + r_p) / 2
r = 400 + r_earth
E1 = np.arccos(ecc**-1 - r/(a*ecc))
E2 = E1 + np.pi
M1 = E1 - ecc * np.sin(E1)
M2 = E2 - ecc * np.sin(E2)
delM = M2 - M1
n = np.sqrt(mu_earth/a**3)
delT = delM / n
print(f"Time spent about 400 km altitude in orbit: {delT/(60)} minutes\n\n\n")

""" Curtis 3.10 """
print("\nCurtis Exercise 3.10:\n")
T = 14*3600
n = (2*np.pi) / T
t = 10*3600
r_p = 10000 #km
a = (mu_earth*(T/(2*np.pi))**2)**(1/3)
ecc = 1 - r_p/a
Me = n * t
if Me > np.pi:
    E = Me - ecc/2
else:
    E = Me + ecc/2
tol = 1
while tol > 10e-8:
    fE = Me - E + ecc * np.sin(E)
    dfE = ecc * np.cos(E) - 1
    E_new = E - fE / dfE
    tol = np.abs(E - E_new)
    E = E_new
r = a * (1 - ecc * np.cos(E))
v = np.sqrt(mu_earth*(2/r - 1/a))
v_r = (np.sqrt(mu_earth/a)*ecc*np.sin(E)) / (1 - ecc*np.cos(E))
print(f"Radial position: {r} km\nVelocity magnitude: {v} km/s\nRadial velocity: {v_r} km/s\n\n\n")

""" Curtis 3.20 """
print("\nCurtis Exercise 3.20:\n")
# Check answer with using ODE45 as well, do they agree?
# Please note, answers in book for 3.10 and 3.20 are slightly off (1-2 km)
r0 = np.array([20000, 105000, 19000])
v0 = np.array([.9, 3.4, 1.5])
dt = 2*3600
r, v, iterations = of.ECI_dt(r0, v0, dt) # type: ignore
print(f"Universal anomaly iterations: {iterations}")
print(f"Using the Universal Anomaly method,\n\nPosition: {r} km\nVelocity: {v} km/s\n\n")

sol = of.ODEprimer([20000,105000,19000],[.9,3.4,1.5], [0, 2*3600], mu_earth)
rODE = sol.y[0:3, -1]
vODE = sol.y[3:6, -1]
print(f"Using ODE45,\n\nPosition: {rODE} km\nVelocity: {vODE} km/s\n\n")

r_err_vec = r - rODE
v_err_vec = v - vODE
print(f"Differences between the two methods,\n\nPosition error: {r_err_vec} km\nVelocity error: {v_err_vec} km/s\n\n\n")

""" Curtis 4.5 """
print("\nCurtis Exercise 4.5:\n")
h_vec, ecc_vec, inc, RAAN, argp, TA = of.ECI2COEs([6500, -7500, -2500, 4, 3, -3])
print(f"H vector: {h_vec} km**2/s\nH magnitude: {np.linalg.norm(h_vec)} km**2/s\necc vector: {ecc_vec}\n")
print(f"Radius: {np.linalg.norm([6500, -7500, -2500])} km\n")
print(f"Velocity: {np.linalg.norm([4, 3, -3])} km/s\n")
print(f"Specific angular momentum: {round(np.linalg.norm(h_vec))} km**2/s\nEccentricity: {round(np.linalg.norm(ecc_vec), 4)}\nInclination: {round(np.rad2deg(inc), 2)} degrees\nRAAN: {round(np.rad2deg(RAAN), 1)} degrees\nArgument of perigee: {round(np.rad2deg(argp), 2)} degrees\nTrue anomaly: {round(np.rad2deg(TA), 1)} degrees\n\n\n")


""" Curtis 4.7 """
print("\nCurtis Exercise 4.7:\n")
r_vec = np.array([-6600, -1300, -5200])
e_vec = np.array([-.4, -.5, -.6])
r = np.linalg.norm(r_vec)
ecc = np.linalg.norm(e_vec)
TA = 2*np.pi - np.arccos((np.dot(e_vec, r_vec)/(ecc*r)))
h = np.sqrt(r*mu_earth*(1+ecc*np.cos(TA)))
v_r = (mu_earth/h)*ecc*np.sin(TA)
v_az = (mu_earth/h)*(1+ecc*np.cos(TA))
v = np.sqrt(v_r**2 + v_az**2)
v_vec = (1/(r*v_r))*((v**2 - (mu_earth/r))*r_vec - mu_earth*e_vec)
h_vec = np.cross(r_vec, v_vec)
inc = np.arccos(h_vec[2]/np.linalg.norm(h_vec))
print(f"Inclination: {round(np.rad2deg(inc), 2)} degrees\n")