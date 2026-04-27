import libFRC_v1_1 as frc
import numpy as np
import math as m

# Geometry of laminate
Lx = 1
Ly = 2
t = 3e-3

# Strengths for carbon-epoxy
sig_xt = 2e9
sig_xc= 1.1e9
sig_yt = 0.08e9
sig_yc = 0.28e9
sig_tau = 0.12e9

# Carbon-epoxy material characteristics
E1 = 140e9
E2 = 10e9
G12 = 7e9
nu12 = .3

# Ply angles
theta = [0.0, 30.0, -30.0, 45.0]

# Stress transformation matrix
def Ts(theta_deg: float) -> np.ndarray:
    th = np.deg2rad(theta_deg)
    c = np.cos(th)
    s = np.sin(th)
    return np.array([
        [c*c, s*s, 2*c*s],
        [s*s, c*c, -2*c*s],
        [-c*s, c*s, c*c - s*s]
    ], dtype=float)
    
# Tsai-Hill Failure criteria function
def tsaiHill(s11, s22, s12, Y1T, Y1C, Y2T, Y2C, Y12):
    X1 = Y1T if s11 >= 0.0 else Y1C
    X2 = Y2T if s22 >= 0.0 else Y2C
    return (s11/X1)**2 + (s22/X2)**2 + (s12/Y12)**2 - (s11*s22)/(X1**2) - 1.0

def failNM(N, M):
    # ---- 1. Build N & M
    N = np.asarray(N, dtype=float).reshape(3)
    M = np.asarray(M, dtype=float)
    
    # ---- 2. Define laminate layers
    tp = t/8
    z = np.array([-t/2 + i*tp for i in range(5)], dtype=float)
    # Ply
    ply = []
    for i, ang in enumerate(theta):
        ply.append(frc.Ply(theta=ang, zL=z[i], zH=z[i+1], E1=E1, E2=E2, G12=G12, nu12=nu12))
    # Qbar
    Qb = [p.Qbar for p in ply]
    
    # ---- 3. Build A, B, and D
    A = np.zeros((3, 3), dtype=float)
    D = np.zeros((3, 3), dtype=float)
    
    for i in range(4):
        A += Qb[i] * 2 * (z[i+1] - z[i])
        D += Qb[i] * (2/3) * (z[i+1]**3 - z[i]**3)
    
    B = np.zeros_like(A) # Symmetric
    
    # ---- 4. Solve midplane strains/curvature
    eps0 = np.linalg.solve(A, N) # [eps_x0, eps_y0, eps_xy0]
    kappa = np.linalg.solve(D, M) # [kappa_x, kappa_y, kappa_xy]
    
    # ---- 5. Ply midplae strains
    zmid = np.array([(z[i] + z[i+1]) / 2 for i in range(4)], dtype=float)
    eps_ply = [eps0 - zmid[i]*kappa for i in range(4)]
    
    # ---- 6. Global ply stress
    sigGlobal = [Qb[i] @ eps_ply[i] for i in range(4)]
    
    # ---- 7. Transform to local and evaluate failure conditions
    f = np.zeros(4, dtype=float)
    for i, ang in enumerate(theta):
        sigLocal = Ts(ang) @ sigGlobal[i]
        f[i] = tsaiHill(*sigLocal, sig_xt, sig_xc, sig_yt, sig_yc, sig_tau)
    
    # Because symmetric, top-half plies should have the same magnitudes (M=0)
    return f

def solvePMax(
    P_upper=1e7,
    Ndir=(1.0/Ly, 1.0/Lx, 0),
    Mdir=(0.0, 0.0, 0.0),
    N0=(0.0,0.0,0.0),
    M0=(0.0,0.0,0.0),
    tol=1e-6,
    maxIter=100
):
    
    Ndir = np.asarray(Ndir, dtype=float).reshape(3)
    Mdir = np.asarray(Mdir, dtype=float).reshape(3)
    N0 = np.asarray(N0, dtype=float).reshape(3)
    M0 = np.asarray(M0, dtype=float).reshape(3)
    
    def g(P):
        N = N0 + P * Ndir
        M = M0 + P * Mdir
        return float(np.max(failNM(N, M)))
    
    P_lower = 0
    
    g_lower = g(P_lower)
    g_upper = g(P_upper)
    
    bracketingIter = 0
    while g_upper < 0.0 and bracketingIter < maxIter:
        P_upper *= 2
        g_upper = g(P_upper)
        bracketingIter += 1
        
    if g_upper < 0.0:
        raise RuntimeError("Failed bracketing procedure. Increase P_upper")
    
    # Bisection
    for _ in range(maxIter):
        P_mid = 0.5*(P_lower + P_upper)
        g_mid = g(P_mid)
        
        if abs(g_mid) < tol:
            return P_mid
        
        if g_mid < 0.0:
            P_lower = P_mid
        else:
            P_upper = P_mid
            
    return 0.5*(P_lower + P_upper)

def main():
    # === EXERCISE 6 ============================================================
    # Uniform biaxial in-plane loading, no moments
    Pmax = solvePMax(
        P_upper=1e6,
        Ndir=(1.0/Ly, 1.0/Lx, 0.0),
        Mdir=(0.0, 0.0, 0.0)
    )
    print(f"\n\n\nExercise 6:\nMaximum load: {Pmax:.3e} N\n\n")
    
    # === EXERCISE 7 ============================================================
    # Uniform biaxial moments, no in-plane loading
    Pmax_moment = solvePMax(
        P_upper=1e6,
        Ndir=(0.0, 0.0, 0.0),
        Mdir=(1.0/Ly, 1.0/Lx, 0.0)
    )
    print(f"Exercise 7:\nMaximum moment: {Pmax_moment:.3e} N-m")
    
    
if __name__ == "__main__":
    main()