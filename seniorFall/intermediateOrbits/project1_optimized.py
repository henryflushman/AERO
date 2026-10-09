"""V-bar rendezvous: phasing -> football -> hops (1 km, 300 m, 20 m) -> V-bar approach.

LVLH frame: x = radial (R-bar), y = along-track (V-bar), z = cross-track.
Units: km, km/s, s.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from IntermediateOrbitsTools import ClassicalOrbitalElements
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
PHASING_REVS = 1  # number of revolutions made before reaching -25 km on V-bar

FOOTBALL_AXES_KM = np.array([40.0, 20.0])  # semi-axes [along-track, radial]; CW forces 2:1

HOP_TARGETS_KM = [1.0, 0.3, 0.02]
HOP_REVS = [1, 1, 1]  # revolutions per hop (to be changed later)
VBAR_RATE_KM_S = 1e-5
APPROACH_STEP_S = 200.0  # time between corrective burns during the V-bar approach

# ==================================


def propagate_cw(state, elapsed_s):
    """Propagate a relative LVLH state with the closed-form CW solution."""
    return cwhill_state_transition_matrix(elapsed_s, N_TARGET) @ state


def initial_phasing_state():
    """Start on the V-bar with an along-track velocity that drifts to PHASING_END_Y_KM.

    With x0 = xdot0 = 0, CW gives y(T) = y0 - 3*T*ydot0 after one period, so the
    chaser must be in a slightly lower (faster) orbit, i.e. ydot0 < 0.
    """
    duration = PHASING_REVS * PERIOD_S
    ydot0 = (START_Y_KM - PHASING_END_Y_KM) / (3.0 * duration)
    return np.array([0.0, START_Y_KM, 0.0, 0.0, ydot0, 0.0])


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
    """Constant-rate V-bar approach to y = 0 using CW-targeted burns each step.

    Every step solves for the velocity that carries the chaser from its current
    position to the next V-bar point (0, y - rate*step) in step_s, then re-burns.
    Returns the final state, the burn list, and the position history.
    """
    y_start = state[1]
    total_s = y_start / rate_km_s
    steps = max(1, int(np.ceil(total_s / step_s)))
    dt = total_s / steps
    phi = cwhill_state_transition_matrix(dt, N_TARGET)
    phi_rr, phi_rv = phi[:3, :3], phi[:3, 3:]

    burns, history = [], [state[:3].copy()]
    for k in range(steps):
        goal = np.array([0.0, y_start * (1.0 - (k + 1) / steps), 0.0])
        v_required = np.linalg.solve(phi_rv, goal - phi_rr @ state[:3])
        burns.append(v_required - state[3:])
        state = propagate_cw(np.concatenate([state[:3], v_required]), dt)
        history.append(state[:3].copy())
    # Final burn to stop at contact
    burns.append(-state[3:])
    state = state.copy()
    state[3:] = 0.0
    return state, burns, np.array(history), total_s


def sample(state, duration_s, points=400):
    """Sample a coasting arc; returns times (s) and states."""
    times = np.linspace(0.0, duration_s, points)
    return times, np.array([propagate_cw(state, t) for t in times])


def run_mission():
    """Fly the full mission; returns arc segments and a burn log."""
    segments, burns = [], []  # segments: (label, t0, times, states)

    def add_segment(label, t0, times, states):
        segments.append((label, t0, t0 + times, states))
        return t0 + times[-1]

    state0 = np.array([0.0, START_Y_KM, 0.0, 0.0, 0.0, 0.0])
    arrived, b1, b2, duration = vbar_hop(state0, PHASING_END_Y_KM, PHASING_REVS)
    post_burn = state0.copy()
    post_burn[3:] += b1
    burns.append((f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m (burn 1)", 0.0, b1))
    t = add_segment(f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m", 0.0, *sample(post_burn, duration))
    burns.append((f"Hop to {PHASING_END_Y_KM * 1e3:.0f} m (burn 2)", t, b2))
    state = arrived

    football_state, dv = football_burn(state)
    burns.append(("Football entry", t, dv))
    t = add_segment("Football (half rev)", t, *sample(football_state, PERIOD_S / 2.0))
    state = propagate_cw(football_state, PERIOD_S / 2.0)
    
    if np.isclose(state[0], 0.0):
        state[0] = 0.0
    else:
        raise ValueError("[End of Football Orbit] Rbar nominal value must be close to 0.0")

    for target_y, revs in zip(HOP_TARGETS_KM, HOP_REVS):
        label = f"Hop to {target_y * 1e3:.0f} m"
        arrived, b1, b2, duration = vbar_hop(state, target_y, revs)
        post_burn = state.copy()
        post_burn[3:] += b1
        burns.append((label + " (burn 1)", t, b1))
        t = add_segment(label, t, *sample(post_burn, duration))
        burns.append((label + " (burn 2)", t, b2))
        state = arrived

    state, approach_burns, history, approach_s = vbar_approach(state)
    dt = approach_s / (len(history) - 1)
    for k, b in enumerate(approach_burns):
        burns.append(("V-bar approach", t + k * dt, b))
    add_segment("V-bar approach", t, np.linspace(0.0, approach_s, len(history)), history)
    return segments, burns


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


if __name__ == "__main__":
    segments, burns = run_mission()

    table = pd.DataFrame(
        [(name, t / 3600.0, *(dv * 1e3), np.linalg.norm(dv) * 1e3) for name, t, dv in burns],
        columns=["Burn", "t (hr)", "dVx (m/s)", "dVy (m/s)", "dVz (m/s)", "|dV| (m/s)"],
    )
    print(table.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    print("\nDelta-v by phase (m/s):")
    phase = table["Burn"].str.replace(r" \(burn \d\)", "", regex=True)
    print(table.groupby(phase, sort=False)["|dV| (m/s)"].sum().to_string(float_format=lambda v: f"{v:.5f}"))
    print(f"Mission total dV: {table['|dV| (m/s)'].sum():.4f} m/s")
    print("Final position (km):", segments[-1][3][-1, :3])

    plot_mission(segments)
    plt.show()
