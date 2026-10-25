"""Generate an interactive Cesium scene for the optimized rendezvous mission."""

import os
from pathlib import Path

import numpy as np

from IntermediateOrbitsTools import (
    CesiumScene,
    lvlh_to_eci,
    propagate_two_body,
    times_from_seconds,
)
from project1_optimized import (
    ECI_PLOT_STEP_S,
    RELATIVE_SCALE,
    load_belinter_sat1_state,
    load_belinter_sat1_tle,
    optimize_mission,
)


OUTPUT_DIR = Path(__file__).resolve().parent / "project1"
SCENE_NAME = "BELINTERSAT 1 Rendezvous Mission"
SEGMENT_COLORS = [
    (255, 120, 40, 255),
    (235, 70, 70, 255),
    (190, 90, 230, 255),
    (40, 180, 120, 255),
    (245, 190, 40, 255),
]
TARGET_COLOR = (40, 150, 255, 255)


def main():
    segments, _, _, rev_numbers = optimize_mission()
    target_initial_state = load_belinter_sat1_state()
    epoch = load_belinter_sat1_tle().epoch
    scene = CesiumScene(SCENE_NAME, multiplier=600.0)

    for index, (label, _, segment_times, segment_states) in enumerate(segments):
        if label.startswith("Hold at "):
            duration_s = segment_times[-1] - segment_times[0]
            sample_count = max(2, int(np.ceil(duration_s / ECI_PLOT_STEP_S)) + 1)
            times_s = np.linspace(segment_times[0], segment_times[-1], sample_count)
            relative_positions = np.repeat(
                segment_states[:1, :3], sample_count, axis=0
            )
        else:
            times_s = segment_times
            relative_positions = segment_states[:, :3]

        target_history = propagate_two_body(
            target_initial_state,
            (float(times_s[0]), float(times_s[-1])),
            times_s=times_s,
            method="DOP853",
        )
        chaser_positions = target_history.positions_km + np.array(
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
        absolute_times = times_from_seconds(epoch, times_s)
        scene.add_trajectory(
            absolute_times,
            target_history.positions_km,
            entity_id=f"target-{index}",
            name=f"BELINTERSAT 1 - {label}",
            color=TARGET_COLOR,
            path_width=2.0,
        )
        scene.add_trajectory(
            absolute_times,
            chaser_positions,
            entity_id=f"chaser-{index}",
            name=f"Chaser - {label} ({RELATIVE_SCALE:g}x relative scale)",
            color=SEGMENT_COLORS[index % len(SEGMENT_COLORS)],
            path_width=3.0,
        )

    czml_path = scene.save_czml(OUTPUT_DIR / "mission_cesium.czml")
    html_path = scene.save_html(
        OUTPUT_DIR / "mission_cesium.html",
        title=f"{SCENE_NAME} (revolutions {rev_numbers})",
        ion_token=os.environ.get("CESIUM_ION_TOKEN"),
    )
    print(f"Saved Cesium CZML: {czml_path}")
    print(f"Saved Cesium viewer: {html_path}")


if __name__ == "__main__":
    main()
