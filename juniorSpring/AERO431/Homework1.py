# =============================================================
# AERO 421 - Aerospace Structural Analysis II
#   o Homework 1
#
# Written by Henry Flushman
# =============================================================

# === Imports =================================================
import numpy as np
import matplotlib.pyplot as plt
# =============================================================

# =============================================================
#                       PROBLEM 1
#        AERO 431 – Aerospace Structural Analysis II
# -------------------------------------------------------------
#        Displacement in a Stepped Bar using
#              Castigliano’s Second Theorem


# --- Given Values --------------------------------------------
L1 = 1.0
L2 = 1.0
r1 = 0.1
r2 = 0.05
E = 70e9
P = 10e3
# -------------------------------------------------------------

# --- Analysis ------------------------------------------------
A1 = np.pi * r1**2
A2 = np.pi * r2**2

# Displacement functions
def u(x):
    x = np.asarray(x)
    out = np.zeros_like(x, dtype=float)
    
    region1 = (x >= 0) & (x <= L1)
    region2 = (x > L1) & (x <= (L1 + L2))
    
    out[region1] = (P / (E * A1)) * x[region1]
    out[region2] = ((P * L1) / (E * A1)) + (P / (E * A2)) * (x[region2] - L1)
    
    return out

# Plotting values
x_vals = np.linspace(0, (L1 + L2), 500)
u_vals = u(x_vals)

# Important values
u_L1 = u(np.array([L1]))[0]
u_tip = u(np.array([L1 + L2]))[0]

print(f"A1 = {A1:.2e} m^2")
print(f"A2 = {A2:.2e} m^2")
print(f"u(L1) = {u_L1:.3e} m")
print(f"u(L1 + L2) = {u_tip:.3e} m")
# -------------------------------------------------------------

# --- Plotting-------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(x_vals, u_vals, linewidth=2)
plt.axvline(L1, linestyle='--')
plt.xlabel('x (m)')
plt.ylabel('u(x) (m)')
plt.title('Displacement along the stepped bar')
plt.grid(True)
plt.tight_layout()
plt.show()
# -------------------------------------------------------------


# ============================================================
#                        PROBLEM 2
# ------------------------------------------------------------
# AERO 431 - Aerospace Structural Analysis II
#
# SymPy derivation of beam deflection w(x) using
# Castigliano's Second Theorem (normal dummy-load method)
# ============================================================

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Symbolic variables
# ------------------------------------------------------------
s, x = sp.symbols('s x', real=True)
L, E, I = sp.symbols('L E I', positive=True, real=True)
P, Q, R = sp.symbols('P Q R', real=True)

# ------------------------------------------------------------
# Full reactions with real loads P, Q and dummy load R included
# ------------------------------------------------------------
A = 2*P/3 - Q/3 + R*(1 - x/L)
B = P/3 - 2*Q/3 + R*x/L

# Given numeric relation in the problem
#   P = 10 kN, Q = 20 kN => Q = 2P
subs_Q = {Q: 2*P}

A_sub = sp.simplify(A.subs(subs_Q))
B_sub = sp.simplify(B.subs(subs_Q))

# ------------------------------------------------------------
# Helper for symbolic cleanup
# ------------------------------------------------------------
def clean(expr):
    """Expand, simplify, and factor an expression."""
    return sp.factor(sp.simplify(sp.expand(expr)))

# ------------------------------------------------------------
# Helper: compute U and w(x) from a list of moment pieces
# ------------------------------------------------------------
def castigliano_deflection(moment_pieces):
    """
    moment_pieces: list of tuples
        [(M1, (s, a1, b1)), (M2, (s, a2, b2)), ...]

    Returns
    -------
    U_simplified : symbolic strain energy
    w_simplified : symbolic deflection = dU/dR at R=0
    """
    U = sp.S(0)
    for M_expr, limits in moment_pieces:
        integrand = sp.expand(M_expr**2) / (2 * E * I)
        U += sp.integrate(integrand, limits)

    U = clean(U)
    w = sp.diff(U, R)
    w = clean(w.subs(R, 0))

    return U, w

# ============================================================
# CASE 1: 0 <= x <= L/3
# Load order:
#   0 ---- x ---- L/3 ---- 2L/3 ---- L
# ============================================================

M1_1 = A_sub * s
M1_2 = A_sub * s - R * (s - x)
M1_3 = A_sub * s - R * (s - x) - P * (s - L/3)
M1_4 = A_sub * s - R * (s - x) - P * (s - L/3) + Q * (s - 2*L/3)

M1_1 = clean(M1_1.subs(subs_Q))
M1_2 = clean(M1_2.subs(subs_Q))
M1_3 = clean(M1_3.subs(subs_Q))
M1_4 = clean(M1_4.subs(subs_Q))

case1_pieces = [
    (M1_1, (s, 0, x)),
    (M1_2, (s, x, L/3)),
    (M1_3, (s, L/3, 2*L/3)),
    (M1_4, (s, 2*L/3, L)),
]

U1, w1 = castigliano_deflection(case1_pieces)

# ============================================================
# CASE 2: L/3 <= x <= 2L/3
# Load order:
#   0 ---- L/3 ---- x ---- 2L/3 ---- L
# ============================================================

M2_1 = A_sub * s
M2_2 = A_sub * s - P * (s - L/3)
M2_3 = A_sub * s - P * (s - L/3) - R * (s - x)
M2_4 = A_sub * s - P * (s - L/3) - R * (s - x) + Q * (s - 2*L/3)

M2_1 = clean(M2_1.subs(subs_Q))
M2_2 = clean(M2_2.subs(subs_Q))
M2_3 = clean(M2_3.subs(subs_Q))
M2_4 = clean(M2_4.subs(subs_Q))

case2_pieces = [
    (M2_1, (s, 0, L/3)),
    (M2_2, (s, L/3, x)),
    (M2_3, (s, x, 2*L/3)),
    (M2_4, (s, 2*L/3, L)),
]

U2, w2 = castigliano_deflection(case2_pieces)

# ============================================================
# CASE 3: 2L/3 <= x <= L
# Load order:
#   0 ---- L/3 ---- 2L/3 ---- x ---- L
# ============================================================

M3_1 = A_sub * s
M3_2 = A_sub * s - P * (s - L/3)
M3_3 = A_sub * s - P * (s - L/3) + Q * (s - 2*L/3)
M3_4 = A_sub * s - P * (s - L/3) + Q * (s - 2*L/3) - R * (s - x)

M3_1 = clean(M3_1.subs(subs_Q))
M3_2 = clean(M3_2.subs(subs_Q))
M3_3 = clean(M3_3.subs(subs_Q))
M3_4 = clean(M3_4.subs(subs_Q))

case3_pieces = [
    (M3_1, (s, 0, L/3)),
    (M3_2, (s, L/3, 2*L/3)),
    (M3_3, (s, 2*L/3, x)),
    (M3_4, (s, x, L)),
]

U3, w3 = castigliano_deflection(case3_pieces)

# ------------------------------------------------------------
# Build final symbolic piecewise deflection function
# ------------------------------------------------------------
w_piecewise = sp.Piecewise(
    (w1, sp.And(x >= 0, x <= L/3)),
    (w2, sp.And(x >= L/3, x <= 2*L/3)),
    (w3, sp.And(x >= 2*L/3, x <= L)),
)

# ------------------------------------------------------------
# Print symbolic results
# ------------------------------------------------------------
print("=" * 70)
print("FULL REACTIONS (WITH P, Q, R)")
print("=" * 70)
print("\nA =")
sp.pprint(clean(A))
print("\nB =")
sp.pprint(clean(B))

print("\n" + "=" * 70)
print("REACTIONS AFTER USING Q = 2P")
print("=" * 70)
print("\nA =")
sp.pprint(A_sub)
print("\nB =")
sp.pprint(B_sub)

print("\n" + "=" * 70)
print("CASE 1: 0 <= x <= L/3")
print("=" * 70)
print("\nU1 =")
sp.pprint(U1)
print("\nw1(x) =")
sp.pprint(sp.collect(sp.expand(w1), x))

print("\n" + "=" * 70)
print("CASE 2: L/3 <= x <= 2L/3")
print("=" * 70)
print("\nU2 =")
sp.pprint(U2)
print("\nw2(x) =")
sp.pprint(sp.collect(sp.expand(w2), x))

print("\n" + "=" * 70)
print("CASE 3: 2L/3 <= x <= L")
print("=" * 70)
print("\nU3 =")
sp.pprint(U3)
print("\nw3(x) =")
sp.pprint(sp.collect(sp.expand(w3), x))

print("\n" + "=" * 70)
print("FINAL PIECEWISE DEFLECTION w(x)")
print("=" * 70)
sp.pprint(w_piecewise)

# ------------------------------------------------------------
# Numeric values for plotting and evaluation
# ------------------------------------------------------------
L_val = 1.0
b_val = 0.1
h_val = 0.05
E_val = 70e9
P_val = 10e3
Q_val = 20e3
I_val = b_val * h_val**3 / 12.0
EI_val = E_val * I_val

numeric_values = {
    L: L_val,
    E: E_val,
    I: I_val,
    P: P_val,
    Q: Q_val,
}

# Numeric symbolic forms
w1_num = sp.simplify(sp.expand(w1.subs(numeric_values)))
w2_num = sp.simplify(sp.expand(w2.subs(numeric_values)))
w3_num = sp.simplify(sp.expand(w3.subs(numeric_values)))

print("\n" + "=" * 70)
print("NUMERICAL DEFLECTION EXPRESSIONS")
print("=" * 70)
print("\nCase 1:")
sp.pprint(sp.collect(w1_num, x))

print("\nCase 2:")
sp.pprint(sp.collect(w2_num, x))

print("\nCase 3:")
sp.pprint(sp.collect(w3_num, x))

# ------------------------------------------------------------
# Lambdify for plotting
# ------------------------------------------------------------
w1_func = sp.lambdify(x, w1_num, 'numpy')
w2_func = sp.lambdify(x, w2_num, 'numpy')
w3_func = sp.lambdify(x, w3_num, 'numpy')

# Real bending moment function (numeric)
def M_real(x_vals):
    x_vals = np.asarray(x_vals, dtype=float)
    M_vals = np.zeros_like(x_vals)

    r1 = (x_vals >= 0.0) & (x_vals <= L_val / 3.0)
    r2 = (x_vals > L_val / 3.0) & (x_vals <= 2.0 * L_val / 3.0)
    r3 = (x_vals > 2.0 * L_val / 3.0) & (x_vals <= L_val)

    M_vals[r1] = 0.0
    M_vals[r2] = -P_val * (x_vals[r2] - L_val / 3.0)
    M_vals[r3] = -P_val * (x_vals[r3] - L_val / 3.0) + Q_val * (x_vals[r3] - 2.0 * L_val / 3.0)

    return M_vals

# Piecewise numeric deflection function
def w_numeric(x_vals):
    x_vals = np.asarray(x_vals, dtype=float)
    w_vals = np.zeros_like(x_vals)

    c1 = (x_vals >= 0.0) & (x_vals <= L_val / 3.0)
    c2 = (x_vals > L_val / 3.0) & (x_vals <= 2.0 * L_val / 3.0)
    c3 = (x_vals > 2.0 * L_val / 3.0) & (x_vals <= L_val)

    w_vals[c1] = w1_func(x_vals[c1])
    w_vals[c2] = w2_func(x_vals[c2])
    w_vals[c3] = w3_func(x_vals[c3])

    return w_vals

# ------------------------------------------------------------
# Maximum bending moment and stress
# ------------------------------------------------------------
x_mmax = 2.0 * L_val / 3.0
Mmax = M_real(np.array([x_mmax]))[0]
sigma_max = abs(Mmax) * (h_val / 2.0) / I_val

# ------------------------------------------------------------
# Maximum deflection location from symbolic derivative
# Use middle region, since that is where the interior extremum occurs
# ------------------------------------------------------------
dw2_dx = sp.diff(w2, x)
critical_points = sp.solve(sp.Eq(dw2_dx, 0), x)

# Keep only real point in [L/3, 2L/3]
critical_points_filtered = []
for cp in critical_points:
    cp_num = sp.N(cp.subs({L: L_val}))
    if cp_num.is_real and (L_val/3.0 <= float(cp_num) <= 2.0*L_val/3.0):
        critical_points_filtered.append(cp)

x_defl_sym = sp.simplify(critical_points_filtered[0])
x_defl_val = float(x_defl_sym.subs({L: L_val}))
w_defl_val = float(w_numeric(np.array([x_defl_val]))[0])

print("\n" + "=" * 70)
print("KEY RESULTS")
print("=" * 70)
print(f"I = {I_val:.8e} m^4")
print(f"EI = {EI_val:.8e} N*m^2")
print(f"Maximum |M| occurs at x = {x_mmax:.8f} m")
print(f"M(2L/3) = {Mmax:.8e} N*m")
print(f"Maximum bending stress magnitude = {sigma_max:.8e} Pa")
print(f"Maximum bending stress magnitude = {sigma_max/1e6:.6f} MPa")

print("\nSymbolic max-deflection location:")
sp.pprint(x_defl_sym)

print(f"\nMaximum downward deflection occurs at x = {x_defl_val:.8f} m")
print(f"w(x) = {w_defl_val:.8e} m")
print(f"|w(x)| = {abs(w_defl_val)*1e3:.6f} mm")

# ------------------------------------------------------------
# Plotting
# ------------------------------------------------------------
x_plot = np.linspace(0.0, L_val, 800)
w_plot = w_numeric(x_plot)
M_plot = M_real(x_plot)

plt.figure(figsize=(8, 5))
plt.plot(x_plot, w_plot, linewidth=2)
plt.axvline(L_val / 3.0, linestyle='--', linewidth=1)
plt.axvline(2.0 * L_val / 3.0, linestyle='--', linewidth=1)
plt.axvline(x_defl_val, linestyle=':', linewidth=1)
plt.xlabel('x (m)')
plt.ylabel('w(x) (m)')
plt.title('Deflection Profile')
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(8, 5))
plt.plot(x_plot, M_plot, linewidth=2)
plt.axvline(L_val / 3.0, linestyle='--', linewidth=1)
plt.axvline(2.0 * L_val / 3.0, linestyle='--', linewidth=1)
plt.xlabel('x (m)')
plt.ylabel('M(x) (N*m)')
plt.title('Bending Moment Diagram')
plt.grid(True)
plt.tight_layout()
plt.show()