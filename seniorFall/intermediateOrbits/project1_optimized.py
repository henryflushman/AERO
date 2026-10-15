"""V-bar rendezvous: phasing -> football -> hops -> V-bar approach.

LVLH frame: x = radial (R-bar), y = along-track (V-bar), z = cross-track.
Units: km, km/s, s. The ECI plot shows the two-body target orbit and a scaled
chaser path reconstructed from the mission's LVLH relative positions.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from IntermediateOrbitsTools import (
    ClassicalOrbitalElements,
    R_EARTH_KM,
    TLERecord,
    lvlh_to_eci,
    propagate_two_body,
    relative_state_lvlh_to_eci,
)
from IntermediateOrbitsTools.RelativeMotion import cwhill_state_transition_matrix


# === CONFIG =======================

TARGET_ORBIT = ClassicalOrbitalElements(
    [42164.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    names=["a", "e", "i", "RAAN", "argp", "true_anomaly"],
    angle_unit="deg",
)
N_TARGET = TARGET_ORBIT.mean_motion_rad_s
PERIOD_S = TARGET_ORBIT.period_s

START_Y_KM = -100.0  # chaser starts behind the target on the V-bar
PHASING_END_Y_KM = -40.0  # V-bar position after revolutions are complete (football apse)

FOOTBALL_AXES_KM = np.array([40.0, 20.0])  # semi-axes [along-track, radial]; CW forces 2:1

HOP_TARGETS_KM = [1.0, 0.3, 0.02]
VBAR_RATE_KM_S = 1e-5
APPROACH_STEP_S = 200.0  # time between corrective burns during the V-bar approach

MAX_REVS = 5  # maximum number of revolutions
MISSION_DURATION_S = 10.0 * 24.0 * 3600.0
MAX_APPROACH_DURATION_S = 3.0 * 3600.0
RELATIVE_SCALE = 1.0  # visual scale applied to the relative ECI offset only
ECI_PLOT_STEP_S = 300.0

# ==================================


def load_belinter_sat1_tle():
    """Load BELINTERSAT 1's TLE."""
    tle_path = Path(__file__).resolve().parent / "project1" / "BELINTERSAT1.csv"
    with tle_path.open(newline="", encoding="utf-8") as tle_file:
        rows = list(csv.DictReader(tle_file))

    if len(rows) != 1:
        raise ValueError(f"Expected one BELINTERSAT 1 TLE row in {tle_path}.")
    row = rows[0]
    if row.get("tle_line0", "").strip() != "0 BELINTERSAT 1":
        raise ValueError(f"Unexpected satellite in {tle_path}.")

    return TLERecord(row["tle_line1"], row["tle_line2"], name="BELINTERSAT 1")


def load_belinter_sat1_state():
    """Convert BELINTERSAT 1's epoch TLE state to osculating elements and ECI."""
    return load_belinter_sat1_tle().coes_at_epoch().to_state()


def propagate_cw(state, elapsed_s):
    """Propagate a relative LVLH state with the closed-form CW solution."""
    return cwhill_state_transition_matrix(elapsed_s, N_TARGET) @ state

def football_burn(state):
    """
    Single burn at the V-bar apse that enters the target-centered football.

    A drift-free CW ellipse has radial semi-axis A and along-track semi-axis 2A,
    so A = FOOTBALL_AXES_KM[1]. At the apse (x = 0) the required velocity is
    xdot = -A*n (moving toward -x) and ydot = -2*n*x = 0. Because the combined
    burn also removes the phasing drift, no separate braking burn is needed.
    """
    a_value = abs(state[1])
    if not np.isclose(a_value, FOOTBALL_AXES_KM[0]):
        raise ValueError("Burn point must be at the football apse (|y| = 2A).")
    new_state = state.copy()
    new_state[3:] = [-a_value * N_TARGET / 2.0, 0.0, 0.0]
    return new_state, new_state[3:] - state[3:]


def vbar_hop(state, target_y_km, revolutions=1):
    """Two-impulse V-bar to V-bar hop that ends at rest at (0, target_y_km).

    Requires x = 0. With xdot0 = 0 the orbit is drift-free in x, and
    y(T) = y0 - 3*T*ydot0, so ydot0 = (y0 - y_target)/(3*T). The second burn
    cancels the residual ydot so the chaser pauses on the V-bar.
    """
    if abs(state[0]) > 1e-9:
        raise ValueError("Hop must start on the V-bar (x = 0).")
    duration = revolutions * PERIOD_S
    ydot0 = (state[1] - target_y_km) / (3.0 * duration)
    after_burn1 = state.copy()
    after_burn1[3:] = [0.0, ydot0, 0.0]
    burn1 = after_burn1[3:] - state[3:]
    arrived = propagate_cw(after_burn1, duration)
    burn2 = -arrived[3:]
    arrived[3:] = 0.0
    return arrived, burn1, burn2, duration


def vbar_approach(state, rate_km_s=VBAR_RATE_KM_S, step_s=APPROACH_STEP_S):
    """Approach y = 0 at a positive constant speed using CW-targeted burns.

    Every step solves for the velocity that carries the chaser to the next
    point on a linear V-bar path from its current position to y = 0.
    Returns the final state, the burn list, and the position history.
    """
    y_start = state[1]
    if not np.isfinite(rate_km_s) or rate_km_s <= 0.0:
        raise ValueError("V-bar approach rate must be a positive speed.")
    total_s = abs(y_start) / rate_km_s
    steps = max(1, int(np.ceil(total_s / step_s)))
    dt = total_s / steps
    phi = cwhill_state_transition_matrix(dt, N_TARGET)
    phi_rr, phi_rv = phi[:3, :3], phi[:3, 3:]

    burns, pre_burn_states, history = [], [], [state[:3].copy()]
    for k in range(steps):
        goal = np.array([0.0, y_start * (1.0 - (k + 1) / steps), 0.0])
        pre_burn_states.append(state.copy())
        v_required = np.linalg.solve(phi_rv, goal - phi_rr @ state[:3])
        burns.append(v_required - state[3:])
        state = propagate_cw(np.concatenate([state[:3], v_required]), dt)
        history.append(state[:3].copy())
    # Final burn to stop at contact
    pre_burn_states.append(state.copy())
    burns.append(-state[3:])
    state = state.copy()
    state[3:] = 0.0
    return state, burns, np.array(history), total_s, pre_burn_states


def sample(state, duration_s, points=400):
    """Sample a coasting arc; returns times (s) and states."""
    times = np.linspace(0.0, duration_s, points)
    return times, np.array([propagate_cw(state, t) for t in times])


def run_mission(rev_numbers: list, extra_time=0.0):
    """Fly the mission, allocating extra seconds to approach then staging holds."""
    if (
        not np.isfinite(extra_time)
        or extra_time < 0.0
        or extra_time > MISSION_DURATION_S
    ):
        raise ValueError("extra_time must be finite and no longer than 10 days.")

    segments, burns = [], []  # burns: (label, time, delta-v LVLH, pre-burn state LVLH)

    def add_segment(label, t0, times, states):
        segments.append((label, t0, t0 + times, states))
        return t0 + times[-1]

    def add_burn(label, time_s, delta_v, state_before):
        burns.append((label, time_s, delta_v.copy(), state_before.copy()))

    base_approach_s = abs(HOP_TARGETS_KM[-1]) / VBAR_RATE_KM_S
    approach_extra_s = min(
        extra_time, max(0.0, MAX_APPROACH_DURATION_S - base_approach_s)
    )
    approach_duration_s = base_approach_s + approach_extra_s
    approach_rate = abs(HOP_TARGETS_KM[-1]) / approach_duration_s
    hold_duration_s = (extra_time - approach_extra_s) / len(HOP_TARGETS_KM)

    state0 = np.array([0.0, START_Y_KM, 0.0, 0.0, 0.0, 0.0])
    arrived, b1, b2, duration = vbar_hop(state0, PHASING_END_Y_KM, rev_numbers[0])
    post_burn = state0.copy()
    post_burn[3:] += b1
    add_burn(f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m (burn 1)", 0.0, b1, state0)
    t = add_segment(f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m", 0.0, *sample(post_burn, duration))
    state_before_b2 = arrived.copy()
    state_before_b2[3:] -= b2
    add_burn(f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m (burn 2)", t, b2, state_before_b2)
    state = arrived

    football_state, dv = football_burn(state)
    add_burn("Football entry", t, dv, state)
    t = add_segment("Football (half rev)", t, *sample(football_state, PERIOD_S / 2.0))
    state = propagate_cw(football_state, PERIOD_S / 2.0)
    
    if np.isclose(state[0], 0.0):
        state[0] = 0.0
    else:
        raise ValueError("[End of Football Orbit] Rbar nominal value must be close to 0.0")

    for target_y, revs in zip(HOP_TARGETS_KM, rev_numbers[1:]):
        label = f"Hop to {target_y * 1e3:.0f} m"
        arrived, b1, b2, duration = vbar_hop(state, target_y, revs)
        post_burn = state.copy()
        post_burn[3:] += b1
        add_burn(label + " (burn 1)", t, b1, state)
        t = add_segment(label, t, *sample(post_burn, duration))
        state_before_b2 = arrived.copy()
        state_before_b2[3:] -= b2
        add_burn(label + " (burn 2)", t, b2, state_before_b2)
        state = arrived
        if hold_duration_s > 0.0:
            hold_label = f"Hold at {target_y * 1e3:.0f} m"
            hold_states = np.vstack((state, state))
            t = add_segment(
                hold_label,
                t,
                np.array([0.0, hold_duration_s]),
                hold_states,
            )

    state, approach_burns, history, approach_s, approach_pre_burn_states = vbar_approach(
        state, rate_km_s=approach_rate
    )
    dt = approach_s / (len(history) - 1)
    for k, (burn, pre_burn_state) in enumerate(
        zip(approach_burns, approach_pre_burn_states)
    ):
        add_burn("V-bar approach", t + k * dt, burn, pre_burn_state)
    tf = add_segment("V-bar approach", t, np.linspace(0.0, approach_s, len(history)), history)
    return segments, burns, tf


def plot_mission(segments):
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    ax_traj, ax_range, ax_small_hops, ax_approach = axes.flat
    for label, _, times, states in segments:
        ax_traj.plot(states[:, 1], states[:, 0], label=label)
        ax_range.plot(times / 3600.0, np.linalg.norm(states[:, :3], axis=1), label=label)
    ax_traj.plot(0.0, 0.0, "k*", markersize=12, label="Target")
    ax_traj.set(xlabel="V-bar, along-track y (km)", ylabel="R-bar, radial x (km)",
                title="LVLH trajectory")
    ax_range.set(xlabel="Time (hr)", ylabel="Range (km)", yscale="symlog",
                 title="Range to target")

    zoom_segments = {
        "Hop to 300 m": ax_small_hops,
        "Hop to 20 m": ax_small_hops,
        "V-bar approach": ax_approach,
    }
    for label, _, _, states in segments:
        ax = zoom_segments.get(label)
        if ax is not None:
            ax.plot(states[:, 1], states[:, 0], label=label)
            ax.plot(states[-1, 1], states[-1, 0], "o", markersize=4)

    ax_small_hops.plot(0.0, 0.0, "k*", markersize=10, label="Target")
    ax_small_hops.set(
        xlabel="V-bar, along-track y (km)",
        ylabel="R-bar, radial x (km)",
        title="Zoom: 300 m and 20 m hops",
    )
    ax_approach.plot(0.0, 0.0, "k*", markersize=10, label="Target")
    ax_approach.set(
        xlabel="V-bar, along-track y (km)",
        ylabel="R-bar, radial x (km)",
        title="Zoom: V-bar approach",
    )
    for ax in (ax_small_hops, ax_approach):
        ax.relim()
        ax.autoscale_view()
        ax.set_aspect("equal", adjustable="datalim")

    for ax in axes.flat:
        ax.grid(True, alpha=0.3)
    ax_traj.legend(fontsize=8)
    ax_small_hops.legend(fontsize=8)
    ax_approach.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_eci_mission(segments, target_initial_state):
    """Create one ECI target/chaser plot for each mission segment."""
    figures = []
    for label, _, times, states in segments:
        if label.startswith("Hold at "):
            duration_s = times[-1] - times[0]
            sample_count = max(2, int(np.ceil(duration_s / ECI_PLOT_STEP_S)) + 1)
            plot_times = np.linspace(times[0], times[-1], sample_count)
            relative_positions = np.repeat(states[:1, :3], sample_count, axis=0)
        else:
            plot_times = times
            relative_positions = states[:, :3]

        target_history = propagate_two_body(
            target_initial_state,
            (float(plot_times[0]), float(plot_times[-1])),
            times_s=plot_times,
            method="DOP853",
        )
        target_positions = target_history.positions_km
        scaled_relative_eci = np.array(
            [
                lvlh_to_eci(
                    relative_position * RELATIVE_SCALE,
                    target_state[:3],
                    target_state[3:],
                )
                for relative_position, target_state in zip(
                    relative_positions, target_history.states
                )
            ]
        )
        chaser_positions = target_positions + scaled_relative_eci

        fig = plt.figure(figsize=(11, 9))
        ax = fig.add_subplot(111, projection="3d")
        ax.plot(
            *target_positions.T,
            label="BELINTERSAT 1 (2BP)",
            color="tab:blue",
            linewidth=1.5,
        )
        ax.plot(
            *chaser_positions.T,
            label=f"Chaser (relative scale {RELATIVE_SCALE:g}x)",
            color="tab:orange",
            linewidth=1.0,
        )

        sphere_angles = np.linspace(0.0, 2.0 * np.pi, 36)
        sphere_latitudes = np.linspace(-np.pi / 2.0, np.pi / 2.0, 18)
        earth_x = R_EARTH_KM * np.outer(np.cos(sphere_angles), np.cos(sphere_latitudes))
        earth_y = R_EARTH_KM * np.outer(np.sin(sphere_angles), np.cos(sphere_latitudes))
        earth_z = R_EARTH_KM * np.outer(np.ones_like(sphere_angles), np.sin(sphere_latitudes))
        ax.plot_surface(earth_x, earth_y, earth_z, color="steelblue", alpha=0.18, linewidth=0)

        plot_points = np.vstack((target_positions, chaser_positions))
        center = (plot_points.max(axis=0) + plot_points.min(axis=0)) / 2.0
        radius = np.max(plot_points.max(axis=0) - plot_points.min(axis=0)) / 2.0
        for axis, midpoint in zip((ax.set_xlim, ax.set_ylim, ax.set_zlim), center):
            axis(midpoint - radius, midpoint + radius)

        ax.set_box_aspect((1.0, 1.0, 1.0))
        ax.set(
            xlabel="ECI x (km)",
            ylabel="ECI y (km)",
            zlabel="ECI z (km)",
            title=f"ECI Trajectory: {label}",
        )
        ax.legend()
        fig.tight_layout()
        figures.append(fig)

    return figures


def maneuver_state_table(burns, target_initial_state):
    """Build ECI target/chaser states and LVLH/ECI delta-v for each maneuver."""
    maneuver_times = np.unique([burn_time for _, burn_time, _, _ in burns])
    target_history = propagate_two_body(
        target_initial_state,
        (float(maneuver_times[0]), float(maneuver_times[-1])),
        times_s=maneuver_times,
        method="DOP853",
    )
    target_states_by_time = {
        float(time_s): state
        for time_s, state in zip(maneuver_times, target_history.states)
    }

    rows = []
    for name, time_s, delta_v_lvlh, chaser_state_lvlh in burns:
        target_state_eci = target_states_by_time[float(time_s)]
        chaser_r_eci, chaser_v_eci = relative_state_lvlh_to_eci(
            target_state_eci[:3],
            target_state_eci[3:],
            chaser_state_lvlh[:3],
            chaser_state_lvlh[3:],
        )
        delta_v_eci = lvlh_to_eci(
            delta_v_lvlh,
            target_state_eci[:3],
            target_state_eci[3:],
        )
        rows.append(
            {
                "Maneuver": name,
                "Time (s)": time_s,
                "Time (hr)": time_s / 3600.0,
                "Target ECI rx (km)": target_state_eci[0],
                "Target ECI ry (km)": target_state_eci[1],
                "Target ECI rz (km)": target_state_eci[2],
                "Target ECI vx (km/s)": target_state_eci[3],
                "Target ECI vy (km/s)": target_state_eci[4],
                "Target ECI vz (km/s)": target_state_eci[5],
                "Chaser ECI rx (km)": chaser_r_eci[0],
                "Chaser ECI ry (km)": chaser_r_eci[1],
                "Chaser ECI rz (km)": chaser_r_eci[2],
                "Chaser ECI vx (km/s)": chaser_v_eci[0],
                "Chaser ECI vy (km/s)": chaser_v_eci[1],
                "Chaser ECI vz (km/s)": chaser_v_eci[2],
                "Chaser LVLH x (km)": chaser_state_lvlh[0],
                "Chaser LVLH y (km)": chaser_state_lvlh[1],
                "Chaser LVLH z (km)": chaser_state_lvlh[2],
                "Chaser LVLH vx (km/s)": chaser_state_lvlh[3],
                "Chaser LVLH vy (km/s)": chaser_state_lvlh[4],
                "Chaser LVLH vz (km/s)": chaser_state_lvlh[5],
                "Burn dV LVLH x (km/s)": delta_v_lvlh[0],
                "Burn dV LVLH y (km/s)": delta_v_lvlh[1],
                "Burn dV LVLH z (km/s)": delta_v_lvlh[2],
                "Burn dV ECI x (km/s)": delta_v_eci[0],
                "Burn dV ECI y (km/s)": delta_v_eci[1],
                "Burn dV ECI z (km/s)": delta_v_eci[2],
            }
        )
    return pd.DataFrame(rows)


def optimize_mission():
    """Select the lowest-delta-v revolution plan that fits within 10 days."""
    all_rev_numbers = []
    all_dvs = []
    all_times = []

    for i in range(MAX_REVS):
        for j in range(MAX_REVS):
            for k in range(MAX_REVS):
                for l in range(MAX_REVS):
                    candidate_revs = [i + 1, j + 1, k + 1, l + 1]
                    _, burns, tf = run_mission(candidate_revs)
                    if tf > MISSION_DURATION_S:
                        continue
                    all_rev_numbers.append(candidate_revs)
                    total_dv = sum(
                        np.linalg.norm(delta_v)
                        for _, _, delta_v, _ in burns
                    )
                    all_dvs.append(total_dv)
                    all_times.append(tf)

    if not all_dvs:
        raise RuntimeError("No revolution plan fits within the 10-day duration.")

    lowest_dv_index = np.argmin(all_dvs)
    rev_num = all_rev_numbers[lowest_dv_index]
    extra_time = MISSION_DURATION_S - all_times[lowest_dv_index]
    segments, burns, tf = run_mission(rev_num, extra_time=extra_time)
    if not np.isclose(tf, MISSION_DURATION_S):
        raise RuntimeError("Selected mission did not reach the required 10-day duration.")
    return segments, burns, tf, rev_num


if __name__ == "__main__":
    segments, burns, tf, rev_num = optimize_mission()

    table = pd.DataFrame(
        [
            (name, t / 3600.0, *(dv * 1e3), np.linalg.norm(dv) * 1e3)
            for name, t, dv, _ in burns
        ],
        columns=["Burn", "t (hr)", "dVx (m/s)", "dVy (m/s)", "dVz (m/s)", "|dV| (m/s)"],
    )
    print(table.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    print("\nDelta-v by phase (m/s):")
    phase = table["Burn"].str.replace(r" \(burn \d\)", "", regex=True)
    print(table.groupby(phase, sort=False)["|dV| (m/s)"].sum().to_string(float_format=lambda v: f"{v:.5f}"))
    print(f"Mission total dV: {table['|dV| (m/s)'].sum():.4f} m/s")
    print(f"Mission duration: {tf / 3600.0:.6f} hr ({tf / 86400.0:.6f} days)")
    print("Final position (km):", segments[-1][3][-1, :3])
    print("Revolution numbers:", rev_num)
    target_initial_state = load_belinter_sat1_state()
    maneuver_table = maneuver_state_table(burns, target_initial_state)
    print("\nTarget/chaser ECI and chaser LVLH states immediately before each burn, plus delta-v:")
    print(maneuver_table.to_string(index=False, float_format=lambda value: f"{value:.8f}"))
    maneuver_csv = Path(__file__).resolve().parent / "project1" / "mission_maneuver_states.csv"
    maneuver_table.to_csv(maneuver_csv, index=False)
    print(f"\nManeuver state and delta-v data saved to: {maneuver_csv}")

    plot_mission(segments)
    plot_eci_mission(segments, target_initial_state)
    plt.show()
