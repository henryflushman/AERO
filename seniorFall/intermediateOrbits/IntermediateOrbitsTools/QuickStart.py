"""Small end-to-end example for the IntermediateOrbitsTools package."""

from datetime import datetime, timezone

import numpy as np

from OrbitalElements import ClassicalOrbitalElements, rv_to_coe
from Propagators import propagate_j2
from CesiumOrbit import CesiumScene, times_from_seconds


# 1) Define an orbit from COEs.
coes = ClassicalOrbitalElements.from_degrees(
    a_km=6878.1363,
    e=0.001,
    i_deg=51.6,
    raan_deg=25.0,
    argp_deg=40.0,
    nu_deg=10.0,
)
state0 = coes.to_state()
print("Initial r [km]   :", state0.r_km)
print("Initial v [km/s] :", state0.v_km_s)

# 2) Propagate six hours with J2.
times_s = np.linspace(0.0, 6.0 * 3600.0, 721)
result = propagate_j2(
    state0,
    (times_s[0], times_s[-1]),
    times_s=times_s,
)

# 3) Convert the final Cartesian state back to COEs.
final_coes = rv_to_coe(result.final_state.r_km, result.final_state.v_km_s)
print("Final COEs:", final_coes.degrees)

# 4) Write a Cesium visualization.
epoch = datetime(2026, 8, 26, 0, 0, tzinfo=timezone.utc)
scene = CesiumScene("J2 Orbit Example", multiplier=120.0)
scene.add_trajectory(
    times_from_seconds(epoch, result.times_s),
    result.positions_km,
    entity_id="satellite",
    name="J2 propagated orbit",
    reference_frame="INERTIAL",
)
scene.save_czml("j2_orbit.czml")
scene.save_html("j2_orbit.html")

print("Wrote j2_orbit.czml and j2_orbit.html")
