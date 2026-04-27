""" Fluids Homework 1 """
import os
os.system('cls')
import numpy as np
import Fluids_Functions as ff
import matplotlib.pyplot as plt
from fpdf import FPDF
from scipy.integrate import quad
from scipy.integrate import solve_ivp

""" Fundamentals of Fluids 1 """
# Constants
visc_0 = 1.716e-5   # Pa.s
T_0 = 273   # Kelvin
S = 110.4  # Kelvin
temp = np.linspace(-50, 500, 1000) + 273.15
# Calculate dynamic viscosity using Sutherland's law
ff.sutherlands_law(temp)
# Plotting
plt.plot(temp - 273.15, ff.sutherlands_law(temp))
plt.xlabel("Temperature (C)")
plt.ylabel("Dynamic Viscosity (Pa.s)")
plt.title("Dynamic Viscosity of Air vs Temperature")
plt.grid(True)
plt.xscale('log')
plt.show()


""" Fundamentals of Fluids 2 """
# Open vessel (1 bar): ice → melt → liquid to 100°C → boil → superheat vapor
bounds  = [0, 100]
cp      = [2.10, 4.18, 2.01]   # kJ/kg-K
latent  = [334, 2260]          # kJ/kg
cv      = [2.10, 4.18, 1.41]   # use liquid cp above 0°C in this simple model
mass = 1.0  # kg
Qp500 = ff.heat_required(0, 500, bounds,  cp,  latent,  mass)  # ≈ 3816 kJ
Qv500 = ff.heat_required(0, 500, bounds, cv, latent, mass)  # ≈ 2424 kJ
T_graph = np.linspace(0, 500, 100)
Qp_graph = [ff.heat_required(0, T, bounds, cp, latent, mass) for T in T_graph]
Qv_graph = [ff.heat_required(0, T, bounds, cv, latent, mass) for T in T_graph]
fig, axs = plt.subplots(1, 2, figsize=(12, 6), sharey=True)
axs[0].plot(Qp_graph, T_graph)
axs[0].set_xlabel('Constant Pressure Energy [kJ]')
axs[0].set_ylabel('Temperature [°C]')
axs[0].grid(True)
axs[1].plot(Qv_graph, T_graph)
axs[1].set_xlabel('Constant Volume Energy [kJ]')
axs[1].grid(True)
plt.show()







""" Fluid Statics """
h = np.linspace(0, 100000, 100)
T = np.zeros_like(h, dtype=float)
P = np.zeros_like(h, dtype=float)
rho = np.zeros_like(h, dtype=float)
visc = np.zeros_like(h, dtype=float)

for i in range(len(h)):
    T[i], P[i], rho[i], visc[i], a = ff.standardatmosphere(h[i])

# Create subplots
fig, axs = plt.subplots(1, 4, figsize=(12, 6), sharey=True)

axs[0].plot(T - 273.15, h)
axs[0].set_xlabel('Temperature [°C]')
axs[0].set_ylabel('Altitude [m]')
axs[0].grid(True)

axs[1].plot(P / 1000, h)  # Convert Pa → kPa
axs[1].set_xlabel('Pressure [kPa]')
axs[1].grid(True)

axs[2].plot(rho, h)
axs[2].set_xlabel('Density [kg/m³]')
axs[2].grid(True)

axs[3].plot(visc, h)
axs[3].set_xlabel('Viscosity [Pa.s]')
axs[3].grid(True)

plt.suptitle('Standard Atmosphere Profiles')
plt.tight_layout()
plt.show()





""" Fluid Statics 2 """
def P_func(z):
    T, P, rho, visc, a = ff.standardatmosphere(float(z))
    return P
L = 100000  # meters
r = 1       # meters
I0, _ = quad(lambda zz: P_func(zz), 0, L)
I1, _ = quad(lambda zz: zz * P_func(zz), 0, L)
# Resultant force and center of pressure
F = 2*np.pi*r * I0
z_cp = I1 / I0
print(z_cp)



""" Fluid Statics 3 """
# Attach to line 1 and line 3
print(ff.height2flowvel(1))
print(ff.height2flowvel(1, 50))




""" Fluid Statics 4 """
# Constants
R = 100000
rho_air = 1.225
rho_h2 = .0899
g = 9.81

# Analysis
print(.75*np.pi*100000**2)




""" Flight Mechanics """
# Constants
g = 9.81
beta = 4800
V0 = 11200
h0 = 100000
def rho_func(h):
    T, P, rho, visc, a = ff.standardatmosphere(float(h))
    return rho
# ODE
def reentry_ode(t, y):
    V, h = y
    dVdt = -(g * rho_func(h) * V**2) / (2 * beta)
    dhdt = -V
    return [dVdt, dhdt]
# Stop when reach sea level

t_span = (0, 4000)
y0 = [V0, h0]

sol = solve_ivp(reentry_ode, t_span, y0, rtol=1e-8, atol=1e-9, max_step=1.0)

t = sol.t
V = sol.y[0]
h = sol.y[1]

plt.figure()
plt.plot(V, h)
plt.xlabel('Speed V [m/s]')
plt.ylabel('Altitude h [m]')
plt.title('Apollo Capsule Entry: Speed–Altitude Map')
plt.show()