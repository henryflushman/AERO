"""
AERO 321 – Lidar Wall Scan / Surface Plot

Talking to Arduino:
    - Python sends:   "START\\n"
    - Arduino replies with many lines: "angle_deg, distance_mm"
    - Arduino ends each sweep with: "SWEEP_DONE"

What this script does:
    1) Grab several sweeps at different vertical heights.
    2) Apply an angle-dependent bias correction to the raw distances.
    3) Average multiple sweeps per height level.
    4) Build a 3D surface of the wall and fit a best-fit plane.
    5) Color the surface by corrected distance.
"""

# ---------------------------- IMPORTS ---------------------------- #

import os
import time
import math
import serial
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (needed for 3D)

os.system('cls' if os.name == 'nt' else 'clear')


# ---------------------------- USER INPUT ------------------------- #

# Serial port settings (change PORT to match your Arduino)
PORT     = "COM6"          # ex: "COM3" on Windows, "/dev/ttyACM0" on Linux
BAUDRATE = 115200

# Vertical sampling
HEIGHT_LEVELS      = 12     # how many vertical levels in the scan
SWEEPS_PER_HEIGHT  = 3      # how many sweeps to average at each height (2–5 is reasonable)
TOTAL_SWEEPS       = HEIGHT_LEVELS * SWEEPS_PER_HEIGHT

# Vertical spacing between scans
Z_STEP_INCH = 0.5
INCH_TO_MM  = 25.4
Z_STEP_MM   = Z_STEP_INCH * INCH_TO_MM

# Range limits (used for clamping out nonsense readings)
MIN_DIST_MM = 50.0
MAX_DIST_MM = 3000.0

# Serial protocol strings
START_COMMAND    = b"START\n"
SWEEP_DONE_TOKEN = "SWEEP_DONE"


# ------------------------ BIAS CALIBRATION ----------------------- #
# Mean residuals from JMP (mm), one per scan angle (deg)
# angle : mean residual (mm)
#  -45  : -72.536
#  -30  : -31.486
#  -15  : -17.759
#    0  : -15.310
#   15  : -19.679
#   30  : -32.818
#   45  : -75.529

CAL_ANGLES_DEG = np.array([-45, -30, -15, 0, 15, 30, 45], dtype=float)
CAL_BIASES_MM  = np.array(
    [-72.536, -31.486, -17.759, -15.310, -19.679, -32.818, -75.529],
    dtype=float
)


def angle_bias(theta_deg):
    """
    Returns the bias (mm) for any input angle (scalar or array),
    using linear interpolation between the calibrated angles.
    """
    return np.interp(theta_deg, CAL_ANGLES_DEG, CAL_BIASES_MM)


# --------------------------- UTILITIES --------------------------- #

def clamp_distance(d):
    """
    Clamp distances into [MIN_DIST_MM, MAX_DIST_MM].

    Notes:
        - If d is <= 0 or not finite, treat it as MAX_DIST_MM
          (essentially a "no return / too far").
    """
    if not math.isfinite(d) or d <= 0:
        return MAX_DIST_MM
    if d < MIN_DIST_MM:
        return MIN_DIST_MM
    if d > MAX_DIST_MM:
        return MAX_DIST_MM
    return d


# ----------------------- SERIAL COMMUNICATION -------------------- #

def get_sweep(ser):
    """
    Trigger a single sweep on the Arduino and return the raw data.

    Returns:
        list of (angle_deg, distance_mm) using the *raw* sensor distance
        (no clamping or bias correction yet).
    """
    # clear anything left in the buffer
    ser.reset_input_buffer()

    # tell Arduino to run one sweep
    ser.write(START_COMMAND)
    ser.flush()

    sweep = []

    while True:
        line = ser.readline().decode("ascii", errors="ignore").strip()
        if not line:
            continue

        if line == SWEEP_DONE_TOKEN:
            # end of sweep
            break

        # expect "angle, distance"
        if "," not in line:
            print(f"  [info] non-data line: {line}")
            continue

        try:
            angle_str, dist_str = line.split(",", 1)
            angle_deg = float(angle_str.strip())
            dist_mm   = float(dist_str.strip())
        except ValueError:
            print(f"  [warn] could not parse line: {line}")
            continue

        sweep.append((angle_deg, dist_mm))

    return sweep


def collect_sweeps():
    """
    Interactively collect TOTAL_SWEEPS sweeps with
    [k]eep / [r]edo / [q]uit logic.

    Each stored sweep is saved as an array with shape (N_points, 2):
        [angle_deg, corrected_distance_mm]
    """
    print(f"Opening serial port {PORT} at {BAUDRATE} baud...")
    ser = serial.Serial(PORT, BAUDRATE, timeout=1)
    time.sleep(2.0)   # give the board time to reset
    print("Serial connected.\n")

    sweeps = []

    try:
        while len(sweeps) < TOTAL_SWEEPS:
            sweep_idx = len(sweeps) + 1
            print(f"-------------------- Sweep {sweep_idx} / {TOTAL_SWEEPS} --------------------")
            input("Press <Enter> to start this sweep...")

            raw_sweep = get_sweep(ser)

            if not raw_sweep:
                print("  [warn] Empty sweep received, trying again.\n")
                continue

            # pull out angles and raw distances
            angles = np.array([a for (a, d) in raw_sweep], dtype=float)
            dists_raw = np.array([d for (a, d) in raw_sweep], dtype=float)

            # clamp obviously bad values
            dists_clamped = np.array([clamp_distance(d) for d in dists_raw], dtype=float)

            # bias correction based on scan angle
            biases = angle_bias(angles)
            dists_corrected = dists_clamped - biases

            # final clamp for safety
            dists_corrected = np.clip(dists_corrected, MIN_DIST_MM, MAX_DIST_MM)

            # quick summary
            print(f"  Points: {len(raw_sweep)}")
            print(f"  Angle range: {angles.min():.1f} deg to {angles.max():.1f} deg")
            print(f"  Raw distances:      {dists_raw.min():.1f} mm to {dists_raw.max():.1f} mm")
            print(f"  Corrected distances:{dists_corrected.min():.1f} mm to {dists_corrected.max():.1f} mm")

            # decide what to do with this sweep
            while True:
                choice = input("Keep this sweep?  [k]eep / [r]edo / [q]uit : ").strip().lower()
                if choice in ("k", "r", "q"):
                    break
                print("  Please type k, r, or q.")

            if choice == "q":
                print("Quitting sweep collection early.\n")
                break
            elif choice == "r":
                print("Discarding this sweep.\n")
                continue
            else:
                # store angle + corrected distance
                sweep_array = np.column_stack((angles, dists_corrected))
                sweeps.append(sweep_array)
                print("Sweep saved.\n")

    finally:
        ser.close()
        print("Serial port closed.\n")

    return sweeps


# --------------------------- SURFACE PLOT ------------------------- #

def build_and_plot_surface(sweeps):
    """
    Take the collected sweeps and produce a surface plot.

    Steps:
        1) Group sweeps into blocks of size SWEEPS_PER_HEIGHT.
        2) Average corrected distances inside each block (one "height level").
        3) Convert polar data (angle, r) -> Cartesian (X, Y, Z).
        4) Fit best-fit plane X = a*Y + b*Z + c.
        5) Plot measured surface + plane, colored by corrected distance.
    """
    if not sweeps:
        print("No sweeps to plot.")
        return

    if len(sweeps) < HEIGHT_LEVELS:
        print(f"[warn] Only {len(sweeps)} sweeps collected, HEIGHT_LEVELS = {HEIGHT_LEVELS}")

    # how many full groups do we have?
    num_groups = len(sweeps) // SWEEPS_PER_HEIGHT
    if num_groups == 0:
        print("[error] Not enough sweeps for even one full height level.")
        return

    if num_groups < HEIGHT_LEVELS:
        print(
            f"[warn] With SWEEPS_PER_HEIGHT = {SWEEPS_PER_HEIGHT}, "
            f"only {num_groups} full height levels are available.\n"
            f"Using {num_groups} height levels."
        )
        height_levels = num_groups
    else:
        height_levels = HEIGHT_LEVELS

    # angle grid from the first sweep (we'll truncate others to match)
    base_angles = sweeps[0][:, 0]
    num_points = len(base_angles)

    # averaged corrected distances per height
    dists_avg = np.zeros((height_levels, num_points), dtype=float)

    for h in range(height_levels):
        group = sweeps[h * SWEEPS_PER_HEIGHT:(h + 1) * SWEEPS_PER_HEIGHT]

        # use the smallest point count in the group
        min_points = min(s.shape[0] for s in group)
        angles_group = group[0][:min_points, 0]

        # stack corrected distances from each sweep and average
        vals = np.stack([g[:min_points, 1] for g in group], axis=0)  # (reps, min_points)
        d_mean = np.mean(vals, axis=0)

        dists_avg[h, :min_points] = d_mean

        if h == 0:
            base_angles = angles_group
            num_points = min_points

    # trim to consistent length
    dists_avg = dists_avg[:, :num_points]
    angles_deg = base_angles[:num_points]
    angles_rad = np.deg2rad(angles_deg)

    # polar -> Cartesian
    cos_t = np.cos(angles_rad)[np.newaxis, :]  # shape (1, N)
    sin_t = np.sin(angles_rad)[np.newaxis, :]

    X = dists_avg * cos_t  # (H, N)   distance "straight out" from sensor
    Y = dists_avg * sin_t  # (H, N)   sideways along wall

    # Z coordinate for each height level
    level_indices = np.arange(height_levels)[:, None]  # (H, 1)
    Z_levels = level_indices * Z_STEP_MM
    Z = np.repeat(Z_levels, num_points, axis=1)        # (H, N)

    # --------------------- Plane Fit: X = a*Y + b*Z + c --------------------- #
    mask = np.isfinite(X)
    X_flat = X[mask]
    Y_flat = Y[mask]
    Z_flat = Z[mask]

    A = np.column_stack((Y_flat, Z_flat, np.ones_like(Y_flat)))
    coef, _, _, _ = np.linalg.lstsq(A, X_flat, rcond=None)
    a, b, c = coef
    X_plane = a * Y + b * Z + c

    # --------------------- Coloring by corrected distance ------------------- #
    dist_min = float(dists_avg.min())
    dist_max = float(dists_avg.max())

    norm = mcolors.Normalize(vmin=dist_min, vmax=dist_max)
    cmap = cm.get_cmap("viridis")

    facecolors = cmap(norm(dists_avg))  # shape (H, N, 4)

    # --------------------- Plotting --------------------- #
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection="3d")

    # measured surface (colored by corrected distance)
    ax.plot_surface(
        X,
        Y,
        Z,
        facecolors=facecolors,
        edgecolor="none",
        antialiased=True,
    )

    # best-fit plane (grey, semi-transparent)
    ax.plot_surface(
        X_plane,
        Y,
        Z,
        color="0.5",
        alpha=0.25,
        edgecolor="none",
    )

    # axis labels
    ax.set_xlabel("X (mm)")   # sensor -> wall
    ax.set_ylabel("Y (mm)")   # sideways
    ax.set_zlabel("Z (mm)")   # height

    # X limits with some padding so small bumps don't look huge
    x_min = float(np.nanmin(X))
    x_max = float(np.nanmax(X))
    x_span = x_max - x_min
    x_pad = max(10.0, 0.2 * x_span)
    ax.set_xlim(x_min - x_pad, x_max + x_pad)

    # Y limits with padding so the ends aren't cramped
    y_min = float(np.nanmin(Y))
    y_max = float(np.nanmax(Y))
    y_span = y_max - y_min
    y_pad = max(10.0, 0.1 * y_span)
    ax.set_ylim(y_min - y_pad, y_max + y_pad)

    ax.set_title(
        f"Bias-Corrected Lidar Surface\n"
        f"{height_levels} height levels, "
        f"{SWEEPS_PER_HEIGHT} sweeps averaged per level"
    )

    # colorbar tied to corrected distance
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.5, aspect=10)
    cbar.set_label("Corrected distance (mm)")

    plt.tight_layout()
    plt.show()


# ------------------------------ MAIN ------------------------------ #

def main():
    sweeps = collect_sweeps()

    if not sweeps:
        print("No sweeps collected. Exiting.")
        return

    build_and_plot_surface(sweeps)


if __name__ == "__main__":
    main()