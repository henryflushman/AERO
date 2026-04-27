import os
os.system('cls')

from Orbits_Functions import julian_date
from Orbits_Functions import sidereal_time
from Orbits_Functions import ODEprimer
from Orbits_Functions import plot_orbit
from matplotlib import pyplot as plt
from numpy import sin, cos, pi
import numpy as np
import datetime

print("Julian Date of 2025-09-22 08:12:00 (with USNO comparison):\n")
print(julian_date(2025, 9, 22, 8, 12, 0, True))

print("\nLocal Sidereal Time examples:\n")
angle_a = 144 + (58/60)
print(f"{sidereal_time(angle_a, 2007,12,21,10,0,0)} (Time from Melbourne)")
print(f"{sidereal_time(-120.653, 2025,7,4,12,30,22)}  (Time from SLO)")

print(f"\nODE Primer:\n")
sol = ODEprimer()
print(sol)
print(f"\nOrbit Plotter:\n")
plot_orbit(sol)
plt.show()

rf_vec = sol.y[0:3, -1]
rf = np.linalg.norm(rf_vec)
print(f"Position magnitude after 5 hours:\n{rf}\n")

vf_vec = sol.y[3:6, -1]
vf = np.linalg.norm(vf_vec)
print(f"Velocity magnitude after 5 hours:\n{vf}\n")


print("\nCurtis Exercise 2.9:\n")
print("The general velocity of a satellite in any type of orbit can be split into two components:\n A radial component and an azimuthal component.\n The radial component is the velocity directed along a line from the center of the Earth to the satellite.\n The azimuthal component is perpendicular to this line and lies in the plane of the orbit.")
print("The radial component of velocity is equivalent to r dot, which is also equal to:\n")
print("v_r = (mu/h) * e * sin(theta)\n")
print("The azimuthal component of velocity is equivalent to h/r, which is also equal to:\n")
print("v_perp = (mu/h) * (1 + e * cos(theta))\n")
print("where mu is the gravitational parameter of the Earth,\n h is the specific angular momentum of the satellite,\n e is the eccentricity of the orbit, and theta is the true anomaly of the satellite.\n")
print("Now using the values of v_r and v_perp, we can find the magnitude of the velocity.\n")
print("v = sqrt(()^2 + v_perp^2)\n")

print("\nCurtis Exercise 2.16:\n")
r = 3390 + 200 #km
mu_mars = 42828 #km^3/s^2
h = np.sqrt(mu_mars * r) #km^2/s
v = h / r #km/s
T = (2*pi*r**(3/2)) / np.sqrt(mu_mars) #s
print(T/(60*60)) #hours
print(f"Orbital period: {datetime.timedelta(seconds=T)} hours")
print(f"Orbital velocity: {v} km/s")

print("\nCurtis Exercise 2.20:\n")

rp = 10000
ra = 100000
mu = 398600 # km^3/s^2
a = (ra + rp) / 2
ecc = (ra - rp) / (ra + rp)
T = (2 * pi * a**(3/2)) / np.sqrt(398600)

print(f"Eccentricity: {ecc}")
print(f"Semi-major axis: {a} km")
print(T/(60*60))
print(f"Orbital period: {datetime.timedelta(seconds=T)} hours")

eta = -.5*(mu/a)
print(f"Specific orbital energy: {eta} km^2/s^2")

h = np.sqrt(mu * rp * (1 + ecc))
tru_anom = np.arccos((h**2 / ((mu * (rp+6378))) - 1) / ecc)
print(f"True anomaly at periapsis: {np.degrees(tru_anom)} degrees")

v_r = (mu/h) * ecc * sin(tru_anom)
v_az = (mu/h) * (1 + ecc * cos(tru_anom))
print(f"Radial velocity at the given true anomaly: {v_r} km/s")
print(f"Azimuthal velocity at the given true anomaly: {v_az} km/s")

v_p = (mu/h) * (1 + ecc * cos(0))
print(f"Velocity at periapsis: {v_p} km/s")
v_a = (mu/h) * (1 + ecc * cos(pi))
print(f"Velocity at apoapsis: {v_a} km/s")