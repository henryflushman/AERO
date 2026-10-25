"""Problem 4: single cylindrical-port hybrid fuel regression.

Run: python hybrid_rocket_problem4.py
Requires: numpy, scipy, matplotlib

Assume uniform radial regression, fixed grain length, and no end-face burning.
The given coefficient is interpreted with oxidizer mass flux in kg/(m^2 s):
    dr/dt = a0 * G_ox**n, with regression rate in m/s.
Thus a0 has units (m/s) / (kg/(m^2 s))**n.
The plot stops at fuel burnout, immediately before any shutdown transient.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp, trapezoid


def main():
    mdot_ox = 4.0       # Constant oxidizer mass flow rate [kg/s]
    rho_f = 920.0      # HTPB fuel density [kg/m^3]
    a0 = 0.0018        # SI regression coefficient, as specified above
    n = 0.3            # Dimensionless regression exponent
    r0 = 0.04 / 2      # Initial port radius [m]
    r_outer = 0.20 / 2 # Grain outer radius [m]
    length = 0.50      # Grain length [m]

    def regression_rate(radius):
        G_ox = mdot_ox / (np.pi * radius**2)  # [kg/(m^2 s)]
        return a0 * G_ox**n                 # [m/s]

    def rhs(time, state):
        return [regression_rate(state[0])]

    def burnout(time, state):
        return state[0] - r_outer

    burnout.terminal = True
    burnout.direction = 1

    # Exact integral provides a time bound and an independent numerical check:
    # r(t)**(2*n+1) = r0**(2*n+1) + (2*n+1)*K*t.
    K = a0 * (mdot_ox / np.pi)**n
    exponent = 2*n + 1
    tb_exact = (r_outer**exponent - r0**exponent) / (exponent * K)

    solution = solve_ivp(
        rhs, (0.0, 1.1 * tb_exact), [r0], events=burnout,
        dense_output=True, rtol=1e-10, atol=1e-12,
    )
    if not solution.success or solution.t_events[0].size != 1:
        raise RuntimeError('Integration failed to reach fuel burnout.')

    tb = solution.t_events[0][0]
    time = np.linspace(0.0, tb, 2001)
    radius = solution.sol(time)[0]
    rdot = regression_rate(radius)
    burning_area = 2 * np.pi * radius * length
    mdot_fuel = rho_f * burning_area * rdot
    mdot_total = mdot_ox + mdot_fuel

    fuel_mass = rho_f * np.pi * (r_outer**2 - r0**2) * length
    integrated_fuel = trapezoid(mdot_fuel, time)
    if not np.isclose(tb, tb_exact, rtol=1e-8):
        raise RuntimeError('Numerical burnout time disagrees with exact result.')
    if not np.isclose(integrated_fuel, fuel_mass, rtol=1e-6):
        raise RuntimeError('Integrated fuel flow disagrees with grain mass.')

    print(f'Burnout time: {tb:.6f} s')
    print(f'Initial/final port diameter: {2*radius[0]:.6f} / {2*radius[-1]:.6f} m')
    print(f'Initial/final regression rate: {1e3*rdot[0]:.6f} / {1e3*rdot[-1]:.6f} mm/s')
    print(f'Initial/final fuel flow: {mdot_fuel[0]:.6f} / {mdot_fuel[-1]:.6f} kg/s')
    print(f'Initial/final total flow: {mdot_total[0]:.6f} / {mdot_total[-1]:.6f} kg/s')
    print(f'Fuel mass from geometry: {fuel_mass:.6f} kg')
    print(f'Fuel mass from integrated flow: {integrated_fuel:.6f} kg')

    plt.rcParams.update({'font.size': 11, 'axes.spines.top': False,
                         'axes.spines.right': False})
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                             layout='constrained')
    axes[0].plot(time, mdot_total, color='#1764a0', lw=2.6,
                 label='Total = oxidizer + fuel')
    axes[0].plot(time, mdot_fuel, color='#d27616', lw=2.2, label='Fuel')

    axes[0].axhline(mdot_ox, color='#6b7280', ls='--', lw=1.5, label='Oxidizer')
    axes[0].set_ylabel('Mass flow rate (kg/s)')
    axes[0].set_ylim(0, 1.12 * mdot_total[-1])
    axes[0].legend(loc='center right', frameon=False)
    axes[0].set_title('Single cylindrical-port hybrid rocket', loc='left',
                       fontsize=15, weight='bold')
    axes[1].plot(time, 1000 * rdot, color='#8b4294', lw=2.6)
    axes[1].set_ylabel('Regression rate (mm/s)')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylim(0, 1.12 * 1000 * rdot[0])
    for ax in axes:
        ax.grid(alpha=0.20)
        ax.set_xlim(0, tb)
        ax.axvline(tb, color='#6b7280', ls=':', lw=1.3)
    axes[1].text(0.98, 0.92, f'Fuel burnout: {tb:.3f} s',
                 transform=axes[1].transAxes, ha='right', va='top')
    plot_path = Path(__file__).resolve().parent / 'hybrid_rocket_problem4.png'
    fig.savefig(plot_path, dpi=220)
    plt.close(fig)
    print(f'Plot saved: {plot_path.name}')


if __name__ == '__main__':
    main()
