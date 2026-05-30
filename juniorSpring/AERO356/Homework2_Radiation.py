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
# ║   Assignment  :  Homework 2 - Radiation                     ║
# ║   Date        :  May 15, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
import numpy as np
import matplotlib.pyplot as plt
from rich.console import Console
from rich.table import Table

console = Console()

# ── Helper functions ──────────────────────
def section(title):
    console.rule(title)
# ─────────────────────────────────────────


# ─── Problem 1 ───────────────────────────
section("Problem 1 — Mean Penetration Depth of Particles in Silicon")

# Givens
silicon_density = 2.3    # g/cm^3
energy = 2.0             # MeV

# Total stopping power
silicon_stopping_power = {
    "electron": 1.567,
    "proton": 111.8,
    "alpha": 1024,
}   # MeV cm^2/g

def meanPenetrationDepth(stoppingPower, density, energy):
    """
    Calculates mean penetration depth in cm.
    """
    return energy / (stoppingPower * density)
    
table = Table()
table.add_column("Particle", justify="left")
table.add_column("Stopping Power (MeV cm^2/g)", justify="right")    
table.add_column("Mean Penetration Depth (cm)", justify="right")

for particle, stoppingPower in silicon_stopping_power.items():
    depth = meanPenetrationDepth(stoppingPower, silicon_density, energy)
    table.add_row(particle, f"{stoppingPower:.3f}", f"{depth:.3f}")

console.print(table)
    

aluminum_density = 2.7    # g/cm^3
aluminum_stopping_power = {
    "electron": 1.518,
    "proton": 109.5,
    "alpha": 985.9,
}   # MeV cm^2/g

graphite_density = 1.7    # g/cm^3
graphite_stopping_power = {
    "electron": 1.619,
    "proton": 140.9,
    "alpha": 1401,
}   # MeV cm^2/g

depth_aluminum = {}
depth_graphite = {}
depth_silicon = {}

for particle in silicon_stopping_power.keys():
    depth_silicon[particle] = meanPenetrationDepth(
        silicon_stopping_power[particle],
        silicon_density,
        energy
    )
    depth_aluminum[particle] = meanPenetrationDepth(
        aluminum_stopping_power[particle],
        aluminum_density,
        energy
    )
    depth_graphite[particle] = meanPenetrationDepth(
        graphite_stopping_power[particle],
        graphite_density,
        energy
    )

particles = list(silicon_stopping_power.keys())

silicon_vals = [depth_silicon[p] for p in particles]
aluminum_vals = [depth_aluminum[p] for p in particles]
graphite_vals = [depth_graphite[p] for p in particles]

x = np.arange(len(particles))
width = 0.25

plt.figure()
plt.bar(x - width, silicon_vals, width, label="Silicon")
plt.bar(x, aluminum_vals, width, label="Aluminum")
plt.bar(x + width, graphite_vals, width, label="Graphite")

plt.xlabel("Particle")
plt.ylabel("Mean Penetration Depth (cm)")
plt.title("Mean Penetration Depth by Particle and Material")
plt.xticks(x, particles)
plt.yscale("log")
plt.legend()
plt.show()

# ─────────────────────────────────────────


# ─── Problem 2 ───────────────────────────
section("Problem 2 - Cumulative Effective Dose")

radiation_weight_factors = {
    "photon": 1,
    "electron": 1,
    "proton": 2,
}

background_dose = {
    "photon": 6,
    "electron": 5,
    "proton": 5,
}   # mGy/month

cme_doses = {
    "May 2032": {
        "photon": 18,
        "electron": 15,
        "proton": 8,
    },
    "Dec 2032": {
        "photon": 21,
        "electron": 19,
        "proton": 9,
    },
}   # mGy/month

def effectiveDose(dose_dict):
    """
    Returns whole-body effective dose in mSv.

    Since exposure is over the whole body, the tissue weighting factors
    sum to 1. For mGy absorbed dose, multiplying by the radiation
    weighting factor gives mSv equivalent/effective dose.
    """
    return sum(
        radiation_weight_factors[rad_type] * dose
        for rad_type, dose in dose_dict.items()
    )

def scaleDose(dose_dict, multiplier):
    """
    Scales all radiation dose components by a multiplier.
    """
    return {
        rad_type: dose * multiplier
        for rad_type, dose in dose_dict.items()
    }

months = []
for year in [2032, 2033, 2034]:
    for month in ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]:
        months.append(f"{month} {year}")

# Mission is Jan 2032 through Dec 2034, so 36 months
months = months[:36]

monthly_results = []

# October 2033 through March 2034.
# October is 1.5x, March is 4.0x because a 300% increase means original + 300% = 4x.
degradation_months = [
    "Oct 2033",
    "Nov 2033",
    "Dec 2033",
    "Jan 2034",
    "Feb 2034",
    "Mar 2034",
]

degradation_multipliers = np.linspace(1.5, 4.0, len(degradation_months))
degradation_lookup = dict(zip(degradation_months, degradation_multipliers))

for i, month in enumerate(months, start=1):
    if month in cme_doses:
        dose_dict = cme_doses[month]
        event_type = "CME"

    elif month in degradation_lookup:
        multiplier = degradation_lookup[month]
        dose_dict = scaleDose(background_dose, multiplier)
        event_type = f"Shielding Degradation, {multiplier:.2f}x"

    else:
        dose_dict = background_dose
        event_type = "Normal"
    
    monthly_effective_dose = effectiveDose(dose_dict)
    
    monthly_results.append({
        "month_number": i,
        "month": month,
        "event": event_type,
        "dose_mSv": monthly_effective_dose,
    })
    
cumulative_mSv = 0
    
for result in monthly_results:
    cumulative_mSv += result["dose_mSv"]
    result["cumulative_mSv"] = cumulative_mSv
    
table = Table()
table.add_column("#", justify="right")
table.add_column("Month", justify="left")
table.add_column("Event", justify="left") 
table.add_column("Monthly Dose (mSv)", justify="right")
table.add_column("Cumulative Dose (mSv)", justify="right")
table.add_column("Cumulative Dose (Sv)", justify="right")

for result in monthly_results:
    table.add_row(
        str(result["month_number"]),
        result["month"],
        result["event"],
        f"{result['dose_mSv']:.2f}",
        f"{result['cumulative_mSv']:.2f}",
        f"{result['cumulative_mSv'] / 1000:.3f}"
    )
    
console.print(table)

total_dose_mSv = monthly_results[-1]["cumulative_mSv"]
total_dose_Sv = total_dose_mSv / 1000

summary_table = Table()
summary_table.add_column("Quantity", justify="left")
summary_table.add_column("Value", justify="right")

summary_table.add_row("Total Effective Dose", f"{total_dose_mSv:.2f} mSv")
summary_table.add_row("Total Effective Dose", f"{total_dose_Sv:.3f} Sv")

console.print(summary_table)

# 45-year-old female astronaut career limit from the class table
dose_limit_mSv = 900

month_exceeded = 0

for result in monthly_results:
    if result["cumulative_mSv"] > dose_limit_mSv:
        month_exceeded = result["month_number"]
        break

limit_table = Table()
limit_table.add_column("Dose Limit Check", justify="left")
limit_table.add_column("Result", justify="right")

limit_table.add_row("Dose Limit", f"{dose_limit_mSv / 1000:.2f} Sv")
limit_table.add_row("Month Exceeded", str(month_exceeded))

console.print(limit_table)

# ─────────────────────────────────────────


# ─── Problem 3 ───────────────────────────

section("Problem 3 - Larmor Radius of an Electron")

# Givens
electron_kinetic_energy_keV = 100   # keV
magnetic_flux_density = 0.15e-4     # T
electron_charge = 1.6e-19           # C
electron_mass = 9.11e-31            # kg

def kineticEnergyKeVToJoules(kinetic_energy_keV, charge):
    """
    Converts kinetic energy from keV to joules.
    """
    return kinetic_energy_keV * 1000 * charge

def velocityFromKineticEnergy(kinetic_energy_J, mass):
    """
    Converts kinetic energy to velocity
    """
    return np.sqrt(2 * kinetic_energy_J / mass)
    
def larmorRadius(mass, perpendicular_velocity, charge, magnetic_flux_density):
    """
    Calculates the larmur radius in meters:
    r_L = m * v_perp / (|q| * B)
    """
    return (mass * perpendicular_velocity) / (abs(charge) * magnetic_flux_density)

electron_kinetic_energy_J = kineticEnergyKeVToJoules(
    electron_kinetic_energy_keV,
    electron_charge
)

electron_velocity = velocityFromKineticEnergy(
    electron_kinetic_energy_J,
    electron_mass
)
    
electron_larmor_radius = larmorRadius(
    electron_mass,
    electron_velocity,
    electron_charge,
    magnetic_flux_density
)

table = Table()

table.add_column("Quantity", justify="left")
table.add_column("Value", justify="right")

table.add_row("Kinetic Energy", f"{electron_kinetic_energy_keV:.2f} keV")
table.add_row("Kinetic Energy", f"{electron_kinetic_energy_J:.3e} J")
table.add_row("Magnetic Flux Density", f"{magnetic_flux_density:.3e} T")
table.add_row("Electron Velocity", f"{electron_velocity:.3e} m/s")
table.add_row("Larmor Radius", f"{electron_larmor_radius:.2f} m")

console.print(table)

# ─────────────────────────────────────────


# ─── Problem 4 ───────────────────────────

section("Problem 4 - Linear Attenuation Coefficient")

I = 1984/30     # cpm
I0 = 2439/10    # cpm
x = 3           # cm

def linearAttenuationCoefficient(cpm_at_detector, cpm_source, distance):
    """
    Solves for the linear attenuation coefficient
    """
    return -(1/distance)*np.log(cpm_at_detector/cpm_source)
    
linear_attenuation_coefficient = linearAttenuationCoefficient(I, I0, x)

table = Table()
table.add_column("Quantity", justify="left")
table.add_column("Value", justify="right")

table.add_row("Count Rate at Detector", f"{I:.2f} cpm")
table.add_row("Count Rate at Source", f"{I0:.2f} cpm")
table.add_row("Distance", f"{x:.2f} cm")
table.add_row("Linear Attenuation Coefficient", f"{linear_attenuation_coefficient:.3f} cm^-1")

console.print(table)

# ─────────────────────────────────────────