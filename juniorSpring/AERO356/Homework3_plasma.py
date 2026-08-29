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
# ║   Course      :  AERO356 - Space Environments II            ║
# ║                                                             ║
# ║   Assignment  :  Homework 3 - Plasma                        ║
# ║   Date        :  June 6, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ===
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


# === Constants ===
kB = 1.38e-23   # J/K
qe = 1.6e-19    # C
eps0 = 8.8e-12  # m-3kg-1s4A2
mi = 1.67e-27   # kg
me = 9.11e-31   # kg
Re = 6378e3     # m
mu = 3.986e14   # m3/s2

# ASSUMED SINGLY IONIZED PLASMA
qi = qe


# === Helper Functions ===

def particleSpeed(electron_temperature):
    return np.sqrt(8 * kB * electron_temperature / (np.pi * me))


def orbitalSpeedCircular(altitude):
    return np.sqrt(mu / (Re + altitude))


def electronCurrent(
    electron_density, 
    electron_velocity, 
    surface_area, 
    voltage, 
    electron_temperature
):
    voltage = np.asarray(voltage)

    coefficient = (
        (1 / 4)
        * qe
        * electron_density
        * electron_velocity
        * surface_area
    )

    current = np.where(
        voltage < 0,
        coefficient * np.exp(qe * voltage / (kB * electron_temperature)),
        coefficient * (1 + qe * voltage / (kB * electron_temperature))
    )

    return current



def ionCurrent(
    ion_density,
    ion_velocity,
    ram_facing_area,
    voltage,
    ion_temperature,
    spacecraft_velocity
):
    voltage = np.asarray(voltage)

    # LEO ram-current approximation
    if ion_velocity < spacecraft_velocity:
        return qi * ion_density * spacecraft_velocity * ram_facing_area * np.ones_like(voltage)

    coefficient = (
        (1 / 4)
        * qi
        * ion_density
        * ion_velocity
        * ram_facing_area
    )

    current = np.where(
        voltage < 0,
        coefficient * (1 - qi * voltage / (kB * ion_temperature)),
        coefficient * np.exp(-qi * voltage / (kB * ion_temperature))
    )

    return current


def debyeLength(
    electron_temperature,
    electron_density
):
    return np.sqrt(eps0 * kB * electron_temperature / (electron_density * qe**2))


def plasmaParameter(
    debye_length,
    electron_density
):
    return (4/3) * np.pi * debye_length**3 * electron_density


# === Problem 1 ===


# Given
r_sc = 0.5
h = 300e3
Te = 1500
Ti = Te
ne = 5e11
ni = ne


# spacecraft speed
v_sc = orbitalSpeedCircular(h)

# particle speeds
v_e_th = particleSpeed(Te)
v_i_th = 5615   # assumed lower than LEO speeds

# collection areas
Ae = 4 * np.pi * r_sc**2
Ai = np.pi * r_sc**2

# voltage range
V = np.linspace(-1, 0, 500)

# current
Ie = electronCurrent(
    electron_density=ne,
    electron_velocity=v_e_th,
    surface_area=Ae,
    voltage=V,
    electron_temperature=Te
)
Ii = ionCurrent(
    ion_density=ni,
    ion_velocity=v_i_th,
    ram_facing_area=Ai,
    voltage=V,
    ion_temperature=Ti,
    spacecraft_velocity=v_sc
)

# equilibrium voltage
Ie0 = (1/4) * qe * ne * v_e_th * Ae
V_eq = (kB * Te / qe) * np.log(Ii[0] / ((1/4) * qe * ne * v_e_th * Ae))


current_difference = Ie - Ii

crossing_indices = np.where(np.diff(np.sign(current_difference)) != 0)[0]

if len(crossing_indices) == 0:
    print("No current intersection found in the plotted voltage range.")
    graph_agrees = False
else:
    i = crossing_indices[0]
    
    V1 = V[i]
    V2 = V[i + 1]
    d1 = current_difference[i]
    d2 = current_difference[i + 1]
    
    V_graph = V1 - d1 * (V2 - V1) / (d2 - d1)
    
    tolerance = 0.01
    graph_agrees = abs(V_graph - V_eq) <= tolerance
    
    
    # part a)
    plt.figure(figsize=(8, 5))
    plt.plot(V, Ie * 1e3, linewidth=2, label="Electron current magnitude")
    plt.plot(V, Ii * 1e3, linewidth=2, label="Ion current magnitude")
    plt.axvline(V_eq, linestyle="--", label=rf"$V_{{eq}}$ = {V_eq:.3f} V")
    plt.scatter(V_graph, Ii[0] * 1e3, zorder=5, label=rf"Graph intersection = {V_graph:.3f} V")

    plt.grid(True)
    plt.xlabel("Spacecraft Voltage, V [V]")
    plt.ylabel("Current Magnitude [mA]")
    plt.title("Electron and Ion Currents for Spherical Spacecraft in LEO Eclipse (Part a)")
    plt.legend()
    plt.tight_layout()

    print("=== QUESTION ONE ==================")

    plt.show()
    
    # part b)
    print("Part b)")
    print(f"    Floating Potential, V_fp = {V_eq:.4f} V")
    
    # part c)
    print("Part c)")
    print(f"    Graph floating potential = {V_graph:.4f} V")
    print(f"    Error = {abs(V_graph - V_eq):.6f} V")
    print(f"    Answer: {graph_agrees}")
    
# part d)
debye = debyeLength(
    electron_temperature=Te,
    electron_density=ne,
)
print("Part d)")
print(f"    Debye Length = {debye:.4f} m")

# part e)
plasma_parameter = plasmaParameter(
    debye_length=debye,
    electron_density=ne
)
print("Part e)")
print(f"    Plasma Parameter = {plasma_parameter:.0f}")

# part f)
alt_geo = 35786e3   # m


# Question has no mention of temperature changing.
# I will assume that temperature changes to common GEO values
Te_geo = 1e7

v_sc_geo = orbitalSpeedCircular(alt_geo)

Ii = ionCurrent(
    ion_density=ni,
    ion_velocity=v_i_th,
    ram_facing_area=Ai,
    voltage=V,
    ion_temperature=Te_geo,
    spacecraft_velocity=v_sc_geo
)

# at float potential I = 0
# so Ie = Ii
# for Vf < 0:
# - Ie = (1/4)*qe*ne*ve*Ae*exp(qe*Vf/(kB*Te))
# - Ii = (1/4)*qi*ni*vi*Ai*(1-qi*Vf/(kB*Ti))
# so:
# - ve*exp(qe*Vf/(kB*Te)) = vi*(1-qe*Vf/(kB*Te))
# thermal speed:
# - v = sqrt(8*kB*T/(pi*m))
# so:
# - vi/ve = sqrt(me/mi)
# and:
# - exp(qe*Vf/(kB*Te)) = sqrt(me/mi)*(1-qe*Vf/(kB*Te))
# let:
# - x = qe*Vf/(kB*Te)
# then:
# - e**x = sqrt(me/mi)(1-x)
# need to run a bisection method on:
# - f(x) = e**x*sqrt(me/mi)*(1-x)

def floatingPotentialGEO(x):
    return np.exp(x) - np.sqrt(me/mi)*(1-x)

x_low = -3
x_high = 3

for _ in range(100):
    x_mid = 0.5 * (x_low + x_high)
    
    if floatingPotentialGEO(x_low) * floatingPotentialGEO(x_mid) <= 0:
        x_high = x_mid
    else:
        x_low = x_mid
    
x = 0.5 * (x_low + x_high)

# because:
# - x = qe*Vf/(kB*Te)
# then:
# - Vf = x*kB*Te/qe

Vf_geo = x * kB * Te_geo / qe

print("Part f)")
print(f"    Floating Potential in GEO = {Vf_geo:.0f}")
print("===================================")
print()


# === QUESTION 3 ============================================


print("=== QUESTION THREE ================")

csv_path = Path("juniorSpring/AERO356/Ionosphere Data Final Exam.csv")
iono_df = pd.read_csv(csv_path)

iono_df.columns = iono_df.columns.str.replace("\ufeff", "", regex=False).str.strip()

iono_df = iono_df.rename(columns={
    iono_df.columns[0]: "altitude_km",
    iono_df.columns[1]: "electron_density_m3"
})

iono_df = iono_df.sort_values("altitude_km").reset_index(drop=True)

iono_df["altitude_m"] = iono_df["altitude_km"] * 1000.0

spacecraft_altitudes_km = [200, 400, 600, 1000, 2000]

frequency_MHz = np.logspace(1, 4, 500)
frequency_Hz = frequency_MHz * 1e6

def calculateTEC(target_altitude_km, df):
    target_altitude_km = float(target_altitude_km)
    
    min_alt = df["altitude_km"].min()
    max_alt = df["altitude_km"].max()
    
    subset = df[df["altitude_km"] <= target_altitude_km].copy()
    
    tec = np.trapezoid(
        subset["electron_density_m3"],
        subset["altitude_m"]
    )
    
    return tec

plt.figure(figsize=(8,5))

all_excess_ranges = []

for h_sc in spacecraft_altitudes_km:
    tec = calculateTEC(h_sc, iono_df)
    
    excess_range_m = 40.31 * tec / frequency_Hz**2
    
    plt.loglog(
        frequency_MHz,
        excess_range_m,
        linewidth=2,
        label=f"{h_sc} km"
    )
    
    all_excess_ranges.append(excess_range_m)
    
plt.grid(True, which="both")
plt.xlabel("Signal Frequency [MHz]")
plt.ylabel("Excess Range [m]")
plt.title("Ionospheric Excess Range vs Signal Frequency")
plt.legend(title="Spacecraft Altitude")
plt.tight_layout()

plt.show()

# part 3c)

all_excess_ranges = np.array(all_excess_ranges)

max_excess_range = np.max(all_excess_ranges)

print("Part c)")
print(f"    Maximum excess range = {max_excess_range:.0f} m")

# part 3d)

max_density = iono_df["electron_density_m3"].max()

max_density_altitude = iono_df.loc[
    iono_df["electron_density_m3"].idxmax(),
    "altitude_km"
]

fp_max_Hz = 8.98 * np.sqrt(max_density)
fp_max_MHz = fp_max_Hz / 1e6

print("Part d)")
print(f"    Max plasma frequency = {fp_max_MHz:.1f} MHz")


# part 3f)

