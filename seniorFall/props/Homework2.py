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
# ║   Course      :  AERO3304 - Propulsions                     ║
# ║                                                             ║
# ║   Assignment  :  Homework 1                                 ║
# ║   Date        :  September 16, 2026                         ║
# ╚═════════════════════════════════════════════════════════════╝


from matplotlib.pyplot import figure
import numpy as np
from scipy.optimize import brentq
from ambiance import Atmosphere


# === CONFIG =============================


# Problem 6
EXHAUST_VELOCITY_P6 = 7000.0  # ft/s
PROPELLANT_MASS_FLOW_RATE_P6 = 280.0  # lbm/s
HEAT_OF_COMBUSTION_P6 = 2400.0  # Btu/lbm

BTU_PER_LBM_TO_FT2_PER_S2_P6 = 778.169 * 32.174

VELOCITY_RATIO_P6 = np.linspace(0, 1, 1000)[1:-1] # trims first and last value to fit :      0.0 < Rv < 1.0

RATED_VELOCITY_P6 = 5000.0 / 7000.0

TOTAL_TIME_P6 = 65  # seconds

TIME_SPAN_P6 = np.linspace(0, TOTAL_TIME_P6, 1000)


# Problem 7
GAMMA_P7 = 1.20
PRESSURE_RATIO_P7 = [2.0, 10.0, 100.0, 1000.0, np.inf]

MACH_SPAN_P7 = np.linspace(1.0, 5.5, 1000)


# Problem 8
ALTITUDE_P8 = 10_000  # m


# ========================================


# === Helper Functions ===================


# Problem 6 functions
def internalEfficiency(velocity_ratio, velocity_exhaust):
    return ((velocity_ratio**2 + 1.0) * velocity_exhaust**2) / (2.0 * HEAT_OF_COMBUSTION_P6 * BTU_PER_LBM_TO_FT2_PER_S2_P6 + (velocity_exhaust * velocity_ratio)**2)

def propulsiveEfficiency(velocity_ratio):
    return (2.0 * velocity_ratio) / (1.0 + velocity_ratio**2)

def totalEfficiency(internal_efficiency, propulsive_efficiency):
    return internal_efficiency * propulsive_efficiency


# Problem 7 functions
def areaRatio(M, gamma):
    return (1/M) * (
        (2/(gamma+1)) *
        (1 + (gamma-1)/2 * M**2)
    ) ** ((gamma+1)/(2*(gamma-1)))
    
def PeOverPc(M, gamma):
    return (
        1.0 + (gamma - 1.0) / 2.0 * M**2
    ) ** (-gamma / (gamma - 1.0))

def thrustCoefficientMomentum(pe_pc, gamma):
    return np.sqrt(
        (2*gamma**2/(gamma-1))
        * (2/(gamma+1))**((gamma+1)/(gamma-1))
        * (
            1
            - pe_pc**((gamma-1)/gamma)
        )
    )

# ========================================




# Problem 6 Analysis
internal_efficiencies = internalEfficiency(VELOCITY_RATIO_P6, EXHAUST_VELOCITY_P6)
propulsive_efficiencies = propulsiveEfficiency(VELOCITY_RATIO_P6)
total_efficiencies = totalEfficiency(internal_efficiencies, propulsive_efficiencies)

figure()
import matplotlib.pyplot as plt
plt.plot(VELOCITY_RATIO_P6, internal_efficiencies, label='Internal Efficiency')
plt.plot(VELOCITY_RATIO_P6, propulsive_efficiencies, label='Propulsive Efficiency')
plt.plot(VELOCITY_RATIO_P6, total_efficiencies, label='Total Efficiency')
plt.axvline(RATED_VELOCITY_P6, color='k', linestyle='--', label='Rated Velocity')
plt.xlabel('Velocity Ratio')
plt.ylabel('Efficiency')
plt.title('Rocket Efficiencies vs Velocity Ratio')
plt.grid()
plt.legend()
plt.show()


# Problem 7 Analysis
plt.figure(figsize=(11, 7))

area_ratio = areaRatio(MACH_SPAN_P7, GAMMA_P7)
pe_pc = PeOverPc(MACH_SPAN_P7, GAMMA_P7)

for R in PRESSURE_RATIO_P7:
    if np.isinf(R):
        p0_pc = 0.0
        label = r"$p_c/p_0=\infty$"
    else:
        p0_pc = 1.0 / R
        label = rf"$p_c/p_0={R}$"
    
    thrust_coefficient_momentum = thrustCoefficientMomentum(pe_pc, GAMMA_P7)
    
    
    thrust_coefficient = thrust_coefficient_momentum + area_ratio * (pe_pc - p0_pc)
    
    if np.isinf(R):
        plt.plot(area_ratio, thrust_coefficient, label=label)
        continue
    
    idx_perfect_expansion = np.argmin(
        np.abs(pe_pc - p0_pc)
    )
    
    plt.plot(
        area_ratio[idx_perfect_expansion],
        thrust_coefficient[idx_perfect_expansion],
        "o"
    )
    
    # Summerfield criterion
    #    pe < 0.4*p0
    separation_limit = 0.4 * p0_pc
    
    # attached and detached masks
    attached = pe_pc >= separation_limit
    separated = pe_pc < separation_limit
    
    plt.semilogx(
        area_ratio[attached],
        thrust_coefficient[attached],
        label=label
    )
    plt.semilogx(
        area_ratio[separated],
        thrust_coefficient[separated],
        linestyle="--"
    )
    
thrust_coefficient_max = np.sqrt(
    (2.0 * GAMMA_P7**2 / (GAMMA_P7 - 1.0))
    * (2.0 / (GAMMA_P7 + 1.0))
    ** ((GAMMA_P7 + 1.0) / (GAMMA_P7 - 1.0))
)

plt.axhline(
    thrust_coefficient_max,
    linestyle=":",
    label=rf"$C_{{F,\max}}={thrust_coefficient_max:.3f}$"
)

plt.xlabel(r"Expansion Ratio, $\epsilon=A_e/A^*$")
plt.ylabel(r"Thrust Coefficient, $C_F$")
plt.title("Thrust Coefficient vs Expansion Ratio")
plt.xlim(1, 180)
plt.ylim(0, 2.32)

plt.grid()
plt.legend()
plt.show()


# Problem 8 Analysis
atm = Atmosphere(ALTITUDE_P8)

pressure_p8 = atm.pressure

print(pressure_p8)