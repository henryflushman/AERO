from pathlib import Path
import numpy as np
import scipy.io as sio

def compute_from_mat(mat_path: Path):
    """Return (filename, cp[24], Cd, Cl, Po, Ps). None if variable P/V missing."""
    d = sio.loadmat(str(mat_path), squeeze_me=True, struct_as_record=False, simplify_cells=True)
    P = d.get("P"); V = d.get("V")
    if P is None or V is None:
        return None
    
    Po = np.mean(P[:, 0])     # adjust indices to your schema
    Ps = np.mean(P[:, 1])
    U  = np.mean(V)
    rho = 1.225
    kin_visc = 1.46e-5
    Re = (rho * U * 0.15875) / kin_visc  # not used later, but OK to keep

    cp    = np.zeros(24)
    P_cyl = np.zeros(24)
    D     = np.zeros(24)
    L     = np.zeros(24)

    for i in range(3, 27):  # i=3..26 (24 taps)
        if np.shape(P) == (100, 26):
            i -= 1
        radians     = (i - 3) * np.pi / 12  # 15° steps
        P_cyl[i-3]  = np.mean(P[:, i])
        cp[i-3]     = (P_cyl[i-3] - Ps) / (Po - Ps)

        F_r = P_cyl[i-3] * 0.15875 * np.pi / 24  # scalar for this tap
        D[i-3] = -F_r * np.cos(radians)
        L[i-3] =  F_r * np.sin(radians)

    D_tot = D.sum()
    L_tot = L.sum()
    Cd = D_tot / ((Po - Ps) * 0.15875)
    Cl = L_tot / ((Po - Ps) * 0.15875)

    return mat_path.name, cp, Cd, Cl, Po, Ps




def download_mat(path_pattern: str):
    """Return first match for path_pattern, or None if not found."""
    matches = list(Path.cwd().glob(f"**/{path_pattern}"))
    return matches[0] if matches else None


import matplotlib.pyplot as plt
import numpy as np

# Assuming `cp` and your loop's angle calculation exist from compute_from_mat()
# If you used radians inside your loop:
theta = np.arange(len(cp)) * 15  # 24 taps × 15° increments = 360°

plt.figure(figsize=(7,5))
plt.plot(theta, cp, 'o-', linewidth=2, markersize=6)
plt.xlabel('θ (degrees)', fontsize=12)
plt.ylabel('$C_p$', fontsize=12)
plt.title('Pressure Coefficient Distribution Around Cylinder', fontsize=14)
plt.grid(True)
plt.gca().invert_yaxis()  # Optional: aerodynamic convention (lower Cp = higher pressure)
plt.tight_layout()
plt.show()