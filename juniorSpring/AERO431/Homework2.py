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
# ║   Course      :  AERO431 - Aerospace Structural Analysis II ║
# ║   Assignment  :  Homework 2                                 ║
# ║   Date        :  April 27, 2026                             ║
# ╚═════════════════════════════════════════════════════════════╝

import sys
from pathlib import Path

libFRC_path = Path(__file__).parent.parent.parent / "juniorWinter" / "AERO331"
sys.path.insert(0, str(libFRC_path))

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from juniorWinter.AERO331.libFRC_v1_1 import Laminate

# Sympy symbols
x = sp.Symbol('x')
y = sp.Symbol('y')

# =============================================================
#  PROBLEM 1 – Simply Supported Beam, Point Load at 3L/4
# =============================================================

L  = 1.0
P  = 1.0
EI = 1.0 / 12.0
x_load = 3 * L / 4

phi_all = [
    sp.sin(sp.pi * x / L),          # type: ignore
    x * sp.sin(sp.pi * x / L),      # type: ignore
    x**2 * sp.sin(sp.pi * x / L),   # type: ignore
    x**3 * sp.sin(sp.pi * x / L),   # type: ignore
]

def a_bilinear_ss(w, dw):
    return float(sp.integrate(EI * sp.diff(w, x, 2) * sp.diff(dw, x, 2), (x, 0, L)))    # type: ignore

def l_linear_ss(dw):
    return float(P * dw.subs(x, x_load))

# Exact (Castigliano)
a_c, b_c = 3*L/4, L/4
exact = P * a_c**2 * b_c**2 / (3 * L * EI)

print("=" * 50)
print("PROBLEM 1 — Simply Supported Beam")
print(f"Exact: δ(3L/4) = {exact:.6f} m")
print("=" * 50)

deflections = {}
w_exprs = {}

for n in [1, 2, 3, 4]:
    phi = phi_all[:n]
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            K[i, j] = a_bilinear_ss(phi[j], phi[i])
            K[j, i] = K[i, j]
    F = np.array([l_linear_ss(phi[i]) for i in range(n)])
    c = np.linalg.solve(K, F)
    w = sum(c[i] * phi[i] for i in range(n))
    w_exprs[n] = w
    delta = float(sp.lambdify(x, w, 'numpy')(x_load))
    deflections[n] = delta
    err = 100 * abs(delta - exact) / exact
    print(f"n={n}: δ = {delta:.6f} m  (error = {err:.2f}%)")

# Plot – deflection shapes
xs = np.linspace(0, L, 300)

def exact_shape(xv):
    b = L - a_c
    return np.where(
        xv <= a_c,
        P*b*xv*(L**2 - b**2 - xv**2) / (6*L*EI),
        P*a_c*(L-xv)*(2*L*xv - xv**2 - a_c**2) / (6*L*EI)
    )

plt.figure(figsize=(7, 4))
for n in [1, 2, 3, 4]:
    ws = sp.lambdify(x, w_exprs[n], 'numpy')(xs)
    plt.plot(xs, 1000*ws, label=f'n={n} ({deflections[n]*1000:.3f} mm)')
plt.plot(xs, 1000*exact_shape(xs), 'k--', label=f'Exact ({exact*1000:.3f} mm)')
plt.gca().invert_yaxis()
plt.xlabel('x (m)')
plt.ylabel('w (mm)')
plt.title('Problem 1: Deflection shapes, load at x = 3L/4')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# =============================================================
#  PROBLEM 2 – Cantilever Wing, Distributed Lift Load
# =============================================================

l   = 1.0
c   = 0.3
t   = 0.01
E   = 70e9          # aluminium (assumed)
EI2 = E * c * t**3 / 12.0

rho = 1.2
S   = l * c
CL  = 1.1

print("\n" + "=" * 50)
print("PROBLEM 2 — Cantilever Wing")
print(f"EI = {EI2:.4f} N·m²  (E = 70 GPa assumed, aluminium)")
print("=" * 50)

def make_basis(n):
    # phi_k = x^(k+2), k=0,1,...  satisfies w(0)=w'(0)=0
    return [x**(k+2) for k in range(n)]

def a_bilinear_c(wi, wj):
    return float(sp.integrate(EI2 * sp.diff(wi, x, 2) * sp.diff(wj, x, 2), (x, 0, l))) # type: ignore

def l_linear_c(dw, v):
    Lx = 0.5 * rho * v**2 * S * CL * (1 - x**2 / l**2)  # type: ignore
    return float(sp.integrate(Lx * dw, (x, 0, l)))      # type: ignore

def solve_ritz_c(n, v):
    phi = make_basis(n)
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            K[i, j] = a_bilinear_c(phi[i], phi[j])
            K[j, i] = K[i, j]
    F = np.array([l_linear_c(phi[i], v) for i in range(n)])
    c_vec = np.linalg.solve(K, F)
    tip = sum(float(c_vec[i] * l**(i+2)) for i in range(n))
    return tip, c_vec, phi

# Part (i): find N at v=20 m/s
v0  = 20.0
tol = 1e-4   # 0.1 mm

print(f"\nPart (i): convergence at v = {v0} m/s")
print(f"{'n':>3}  {'tip (mm)':>10}  {'|Δtip| (mm)':>12}")
print("-" * 30)

tips = {}
tips[1], _, _ = solve_ritz_c(1, v0)
print(f"  1  {tips[1]*1000:>10.4f}  {'—':>12}")

N = None
for n in range(2, 20):
    tips[n], _, _ = solve_ritz_c(n, v0)
    delta = abs(tips[n] - tips[n-1])
    print(f"{n:>3}  {tips[n]*1000:>10.4f}  {delta*1000:>12.4f}{'  ✓' if delta < tol else ''}")
    if delta < tol:
        N = n
        break

print(f"\nConverged at N = {N},  tip = {tips[N]*1000:.4f} mm")

# Part (ii): sweep velocities
velocities = np.arange(2, 32, 2)
tip_mm = []
print(f"\nPart (ii): tip deflection vs velocity (N={N} bases)")
print(f"{'v (m/s)':>8}  {'tip (mm)':>10}")
print("-" * 22)
for v in velocities:
    tip_v, _, _ = solve_ritz_c(N, float(v))
    tip_mm.append(tip_v * 1000)
    print(f"{v:>8.0f}  {tip_v*1000:>10.4f}")

plt.figure(figsize=(6, 4))
plt.plot(velocities, tip_mm, 'o-')
plt.xlabel('Velocity v (m/s)')
plt.ylabel('Tip deflection (mm)')
plt.title(f'Problem 2: Tip deflection vs velocity (N={N} bases)')
plt.grid(True)
plt.tight_layout()
plt.show()

# =============================================================
#  PROBLEM 3 – Cantilever Plate, Distributed Lift Load Along x
# =============================================================
# Plate fixed at (x,y,z) = (0,0,0)

nu = 0.3
D = E * t**3 / (12.0 * (1 - nu**2))
v_plate = 20.0

print("\n" + "=" * 50)
print("PROBLEM 3 — Cantilever Plate")
print(f"D = {D:.4f} N·m  (E=70 GPa, nu=0.3)")
print("=" * 50)

V1 = [x**2]
V2 = [x**2, x**3]
V3 = [x**2, x**2*y**2, x**3, x**3*y**2]   # type:ignore
V4 = [x**2, x**2*y**2, x**2*y**4,   # type:ignore
      x**3, x**3*y**2, x**3*y**4,   # type:ignore
      x**4, x**4*y**2, x**4*y**4]   # type:ignore

spaces = [
    ('V1', V1),
    ('V2', V2),
    ('V3', V3),
    ('V4', V4)
]

def plate_bilinear(wi, wj):
    """Kirchoff plate: bilinear form over [0,1] x [-c/2, c/2]"""
    wi_xx = sp.diff(wi, x, 2)
    wi_yy = sp.diff(wi, y, 2)
    wi_xy = sp.diff(wi, x, y)
    wj_xx = sp.diff(wj, x, 2)
    wj_yy = sp.diff(wj, y, 2)
    wj_xy = sp.diff(wj, x, y)
    term1 = D * (wi_xx + wi_yy) * (wj_xx + wj_yy)   # type: ignore
    term2 = D * (1 - nu) * (wi_xx*wj_yy + wi_yy*wj_xx - 2*wi_xy*wj_xy)  # type: ignore
    integrand = term1 - term2
    return float(sp.integrate(
        sp.integrate(
            integrand,
            (y, -c/2, c/2)
        ),
        (x, 0, 1)
    ))  # type: ignore

def plate_load(dw, v):
    """Functional Load: integral of L(x,y)*dw over plate area."""
    Lxy = (1/2) * rho * v**2 * l * CL * (1 - x**2 / l**2)   # type: ignore
    return float(sp.integrate(
        sp.integrate(
            Lxy * dw,
            (y, -c/2, c/2)
        ),
        (x, 0, 1)
    ))  # type: ignore

def solve_plate(basis, v):
    n = len(basis)
    K = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            K[i, j] = plate_bilinear(basis[i], basis[j])
            K[j, i] = K[i, j]
    F = np.array([plate_load(basis[i], v) for i in range(n)])
    c_vec = np.linalg.solve(K, F)
    w_tip = sum(float(c_vec[i]) * float(basis[i].subs([(x, l), (y, 0)])) for i in range(n)) # type: ignore
    return w_tip, c_vec

print(f"\nv = {v_plate} m/s, tip deflection at (1, 0, 0):")
print(f"{'Space':>4}    {'n':>3}    {'tip (mm)':>10}")
print("-" * 24)

plate_tips = {}
    
for name, basis in spaces:
    print(f"    Computing basis {name} ({len(basis)} basis functions)...", flush=True)
    tip, _ = solve_plate(basis, v_plate)
    plate_tips[name] = tip
    print(f"{name:>4}   {len(basis):>3} {tip*1000:>10.4f}")
    
print(f"\nBeam theory tip deflection (Problem 2, v=20): {tips[N]*1000:.4f} mm")
print("Plate model is stiffer due to Poisson coupling between x and y directions")
print("this makes plate_tip deflection < beam_tip deflection")


# =============================================================
#  PROBLEM 4 – Cantilever Laminate Plate, [θ/0/90]_s
# =============================================================
# Same plate as Q3
# Layup:
#   [θ, 0, 90, 90, 0, θ]

# --- Basis [9 functions] ---
basis_p4 = [
    x**2, x**2*y, x**2*y**2,    # type: ignore
    x**3, x**3*y, x**3*y**2,    # type: ignore
    x**4, x**4*y, x**4*y**2,    # type: ignore
]
nb = len(basis_p4)

phi_xx_p4 = [sp.diff(b, x, 2) for b in basis_p4]
phi_yy_p4 = [sp.diff(b, y, 2) for b in basis_p4]
phi_xy_p4 = [sp.diff(b, x, y) for b in basis_p4]
phi_y_p4  = [sp.diff(b, y)    for b in basis_p4]

def p_int(f):
    return float(sp.integrate(
        sp.integrate(
            f,
            (y, -c/2, c/2)
        ),
        (x, 0, 1)
    ))  # type: ignore
    
print("\n" + "=" * 50)
print("PROBLEM 4 - Cantilever Laminate Plate")
print("Base layup [θ/0/90]_s")
print("=" * 50)
print("Computing geometric derivative and integrals")

G_xxxx = np.zeros((nb, nb)); G_xxyy = np.zeros((nb, nb))
G_xxxy = np.zeros((nb, nb)); G_yyyy = np.zeros((nb, nb))
G_yyxy = np.zeros((nb, nb)); G_xyxy = np.zeros((nb, nb))


for i in range(nb):
    for j in range(nb):
        G_xxxx[i, j] = p_int(phi_xx_p4[i] * phi_xx_p4[j]) # type: ignore
        G_xxyy[i, j] = p_int(phi_xx_p4[i] * phi_yy_p4[j]) # type: ignore
        G_xxxy[i, j] = p_int(phi_xx_p4[i] * phi_xy_p4[j]) # type: ignore
        G_yyyy[i, j] = p_int(phi_yy_p4[i] * phi_yy_p4[j]) # type: ignore
        G_yyxy[i, j] = p_int(phi_yy_p4[i] * phi_xy_p4[j]) # type: ignore
        G_xyxy[i, j] = p_int(phi_xy_p4[i] * phi_xy_p4[j]) # type: ignore
        
v_p4 = 20.0
Lxy = 0.5 * rho * v_p4**2 * l * CL * (1 - x**2 / l**2) # type: ignore
F_p4 = np.array([p_int(Lxy * b) for b in basis_p4])

tip_vals = np.array([float(b.subs([(x, l), (y, 0)])) for b in basis_p4])
tip_dy = np.array([float(d.subs([(x, l), (y, 0)])) for d in phi_y_p4]) # type: ignore

def assemble_K_lam(D):
    D11,D12,D16 = D[0,0],D[0,1],D[0,2]
    D22,D26,D66 = D[1,1],D[1,2],D[2,2]
    return (D11 * G_xxxx
            + D12 * (G_xxyy + G_xxyy.T)
            + 2*D16 * (G_xxxy + G_xxxy.T)
            + D22 * G_yyyy
            + 2*D26 * (G_yyxy + G_yyxy.T)
            + 4*D66 * G_xyxy)
    
thetas = np.arange(-90, 91, 15)
deflections_p4 = []
twists_p4 = []

print(f"\n{'θ':>6}  {'w(l,0) mm':>12}   {'∂w/∂y(l,0)':>12}")
print("-" * 35)

for th in thetas:
    lam = Laminate(theta=[th, 0, 90, 90, 0, th], t=t)
    K = assemble_K_lam(lam.D)
    c_vec = np.linalg.solve(K, F_p4)
    w_tip = float(c_vec @ tip_vals)
    dw_tip = float(c_vec @ tip_dy)
    deflections_p4.append(w_tip * 1000)
    twists_p4.append(dw_tip)
    print(f"{th:>6} {w_tip*1000:>12.4f} {dw_tip:>12.6f}")
    
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

ax1.plot(thetas, deflections_p4, 'o-')
ax1.set_xlabel('θ')
ax1.set_ylabel('Tip deflection w(l,0)   [mm]')
ax1.set_title('Problem 4: Tip deflection vs θ')
ax1.grid(True)

ax2.plot(thetas, twists_p4, 's-')
ax2.set_xlabel('θ')
ax2.set_ylabel('Twist ∂w/∂y(l,0)    [m/m]')
ax2.set_title('Problem 4: Tip twist vs θ')

plt.tight_layout()
plt.show()