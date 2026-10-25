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
# ║   Course      :  AERO4452 - Intermediate Orbits             ║
# ║                                                             ║
# ║   Assignment  :  Homework 4                                 ║
# ║   Date        :  October 18, 2026                           ║
# ╚═════════════════════════════════════════════════════════════╝

import numpy as np
import matplotlib.pyplot as plt

# from toolkit
import time

from IntermediateOrbitsTools import (
    ClassicalOrbitalElements,
    ExponentialAtmosphere,
    OrbitingSpacecraft,
    R_EARTH_KM,
    MU_EARTH_KM3_S2,
    altitude_event,
    altitude_crossing_event,
    enckes_method,
    cowells_method,
    rv_to_coe,
)
from IntermediateOrbitsTools.Propagators import propagate_two_body


# Givens
PERIGEE_ALTITUDE_0 = 215  # km
APOGEE_ALTITUDE_0  = 939  # km

RAAN_0 = 340  # deg
INCLINATION_0 = 65.2  # deg
OMEGA_0 = 58  # deg
THETA_0 = 332  # deg

MASS_KG = 100.0
DIAMETER_M = 1.0
RAM_AREA_M2 = np.pi * (DIAMETER_M / 2.0) ** 2
DRAG_COEFFICIENT = 2.2  # ASSUMPTION: not given in the assignment

# ASSUMPTION: single-layer exponential model referenced near perigee altitude
# (replace with the course/Curtis table values if different)
ATMOSPHERE = ExponentialAtmosphere(
    reference_density_kg_m3=2.789e-10,
    reference_altitude_km=200.0,
    scale_height_km=37.105,
)

DURATION_S = 120 * 24 * 3600.0
TIME_STEP_S = 30.0
REENTRY_ALTITUDE_KM = 100.0
ENCKE_RECTIFICATION_TOLERANCE = 1e-6  # |dr|/|r| before the reference is reset
# ===


# Initial orbit
rp_km = R_EARTH_KM + PERIGEE_ALTITUDE_0
ra_km = R_EARTH_KM + APOGEE_ALTITUDE_0
a_km = 0.5 * (rp_km + ra_km)
e = (ra_km - rp_km) / (ra_km + rp_km)

coes_0 = ClassicalOrbitalElements.from_degrees(
    a_km=a_km,
    e=e,
    i_deg=INCLINATION_0,
    raan_deg=RAAN_0,
    argp_deg=OMEGA_0,
    nu_deg=THETA_0,
)
state_0 = coes_0.to_state()

spacecraft = OrbitingSpacecraft(
    orbital_elements=coes_0,
    mass_kg=MASS_KG,
    ram_area_m2=RAM_AREA_M2,
    coefficient_of_drag=DRAG_COEFFICIENT,
    name="HW4 satellite",
)
drag = spacecraft.drag_acceleration_model(ATMOSPHERE)
t_span = (0.0, DURATION_S)


# Two-body (no perturbations)
tic = time.perf_counter()
two_body = propagate_two_body(
    state_0, t_span, times_s=np.arange(0.0, DURATION_S + TIME_STEP_S, TIME_STEP_S)
)
two_body_runtime_s = time.perf_counter() - tic


# Encke's method with drag, stopped at reentry altitude
tic = time.perf_counter()
encke = enckes_method(
    state_0,
    t_span,
    TIME_STEP_S,
    perturbations=[drag],
    rectification_threshold=ENCKE_RECTIFICATION_TOLERANCE,
    events=[altitude_event(REENTRY_ALTITUDE_KM, direction="below")],
)
encke_runtime_s = time.perf_counter() - tic


# Cowell's method with drag, stopped at reentry altitude
tic = time.perf_counter()
cowell = cowells_method(
    state_0,
    t_span,
    TIME_STEP_S,
    perturbations=[drag],
    events=[altitude_crossing_event(REENTRY_ALTITUDE_KM, direction="below")],
)
cowell_runtime_s = time.perf_counter() - tic


# Results
def orbit_history(result):
    """Return perigee/apogee altitudes and angular elements over time."""
    rows = []
    for state in result.states:
        c = rv_to_coe(state[:3], state[3:])
        rows.append(
            (c.rp - R_EARTH_KM, c.ra - R_EARTH_KM, c.raan_rad, c.i_rad, c.argp_rad)
        )
    rows = np.array(rows)
    rows[:, 2:] = np.degrees(rows[:, 2:])
    return rows


two_body_hist = orbit_history(two_body)
encke_hist = orbit_history(encke)
cowell_hist = orbit_history(cowell)
t_two_body_days = two_body.times_s / 86400.0
t_encke_days = encke.times_s / 86400.0
t_cowell_days = cowell.times_s / 86400.0

encke_reentered = encke.raw_solution["terminated_by_event"]
cowell_reentered = cowell.raw_solution.status == 1
print(f"Two-body runtime : {two_body_runtime_s:.2f} s")
print(f"Encke runtime    : {encke_runtime_s:.2f} s")
if encke_reentered:
    print(f"Encke reached {REENTRY_ALTITUDE_KM:.0f} km at {t_encke_days[-1]:.3f} days")
else:
    print(f"Encke did not reach {REENTRY_ALTITUDE_KM:.0f} km in {DURATION_S / 86400:.0f} days")
print(f"Cowell runtime    : {cowell_runtime_s:.2f} s")
if cowell_reentered:
    print(f"Cowell reached {REENTRY_ALTITUDE_KM:.0f} km at {t_cowell_days[-1]:.3f} days")
else:
    print(f"Cowell did not reach {REENTRY_ALTITUDE_KM:.0f} km in {DURATION_S / 86400:.0f} days")
print(f"Final Encke perigee/apogee altitude: "
      f"{encke_hist[-1, 0]:.2f} / {encke_hist[-1, 1]:.2f} km")
print(f"Final Cowell perigee/apogee altitude: "
      f"{cowell_hist[-1, 0]:.2f} / {cowell_hist[-1, 1]:.2f} km")

fig, ax = plt.subplots()
ax.plot(t_two_body_days, two_body_hist[:, 0], "--", label="Two-body perigee")
ax.plot(t_two_body_days, two_body_hist[:, 1], "--", label="Two-body apogee")
ax.plot(t_encke_days, encke_hist[:, 0], label="Encke perigee")
ax.plot(t_encke_days, encke_hist[:, 1], label="Encke apogee")
ax.plot(t_cowell_days, cowell_hist[:, 0], label="Cowell perigee")
ax.plot(t_cowell_days, cowell_hist[:, 1], label="Cowell apogee")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Altitude (km)")
ax.set_title("Perigee and apogee altitude")
ax.grid(True)
ax.legend()

fig, axes = plt.subplots(3, 1, sharex=True)
for ax, col, label in zip(axes, (2, 3, 4), ("RAAN (deg)", "Inclination (deg)", "Arg. of perigee (deg)")):
    ax.plot(t_two_body_days, two_body_hist[:, col], "--", label="Two-body")
    ax.plot(t_encke_days, encke_hist[:, col], label="Encke")
    ax.plot(t_cowell_days, cowell_hist[:, col], label="Cowell")
    ax.set_ylabel(label)
    ax.grid(True)
axes[0].legend()
axes[-1].set_xlabel("Time (days)")
fig.suptitle("Angular elements")

plt.show()

