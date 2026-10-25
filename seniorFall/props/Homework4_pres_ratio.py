from math import sqrt
from scipy.optimize import brentq

gamma = 1.4
epsilon = 20.0  # Exit area / throat area

def pressure_ratio_error(q):
    # q = exit pressure / chamber pressure
    numerator = ((gamma - 1) / 2) * (
        2 / (gamma + 1)
    )**((gamma + 1) / (gamma - 1))

    denominator = q**(2 / gamma) * (
        1 - q**((gamma - 1) / gamma)
    )

    return sqrt(numerator / denominator) - epsilon

# This interval selects the supersonic solution.
q_sonic = (2 / (gamma + 1))**(gamma / (gamma - 1))
pressure_ratio = brentq(pressure_ratio_error, 1e-12, q_sonic)

print(pressure_ratio)  # 0.002619044