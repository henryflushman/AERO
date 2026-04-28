# =============================================================
# AERO 421 - Assignment 0
#   Review of Rigid Body Kinematics and Dynamics
#
# Written by Henry Flushman
# =============================================================

from sympy import symbols, Matrix
import sympy as sp

# ── Helpers ───────────────────────────────────────────────────
def c(ang): return sp.cos(ang)
def s(ang): return sp.sin(ang)

def banner(title):
    width = 60
    print("\n" + "═" * width)
    print(f"  {title}")
    print("═" * width)

def section(label):
    print(f"\n  ┌─ {label} {'─' * (52 - len(label))}")

def show(label, expr):
    print(f"\n  ┌─ {label} {'─' * (52 - len(label))}")
    print(f'  │')
    lines = sp.pretty(expr).splitlines()
    for line in lines:
        print(f"  │  {line}")
    print(f"  └{'─' * 55}")

# =============================================================
# === Problem 10, Chapter 1 ===================================
# =============================================================

banner("Problem 10 — Chapter 1  │  Euler Angle Kinematics")

phi, theta, psi       = symbols('phi theta psi',             real=True)
phi_dot, theta_dot, psi_dot = symbols('phi_dot theta_dot psi_dot', real=True)

# ── Rotation matrices ─────────────────────────────────────────
C_x = Matrix([
    [1,       0,        0      ],
    [0,       c(phi),  -s(phi) ], # type: ignore
    [0,       s(phi),   c(phi) ],
])
C_y = Matrix([
    [ c(theta), 0, s(theta)],
    [ 0,        1, 0       ],
    [-s(theta), 0, c(theta)], # type: ignore
])
C_z = Matrix([
    [c(psi), -s(psi), 0], # type: ignore
    [s(psi),  c(psi), 0],
    [0,       0,      1],
])

# ── Part (a): DCM C21 = Cx · Cz · Cy ─────────────────────────
section("Part (a) — Direction Cosine Matrix  C₂₁ = Cₓ · C_z · C_y")
C21 = sp.simplify(C_x @ C_z @ C_y)
show("C₂₁", C21)



# ── Part (d): Basis vectors of F1 ───────────────────────────────────────
section("Part (d) — Basis vectors of F₁")
x_F1 = Matrix([1, 0, 0])
y_F1 = Matrix([0, 1, 0])
z_F1 = Matrix([0, 0, 1])

# ── Angular velocity components ───────────────────────────────
w_phi   = phi_dot   * x_F1
w_psi   = psi_dot   * (C_x  @ z_F1)
w_theta = theta_dot * (C_x  @ C_z @ y_F1)

omega = sp.simplify(w_phi + w_theta + w_psi)
show("ω  (angular velocity vector)", omega)

# ── Kinematic matrix B  where  ω = B · q̇ ────────────────────
qdot = Matrix([phi_dot, theta_dot, psi_dot])
B    = sp.simplify(omega.jacobian(qdot))
show("B  (kinematic matrix,  ω = B·q̇)", B)

print()