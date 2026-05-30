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
# ║   Assignment  :  Homework 3                                 ║
# ║   Date        :  May 11, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp


# === Problem 1 ===

# --- Givens ------------------------------
L = 1.0
P = -1.0
Ebh3 = 1.0

EI = Ebh3 / 12.0

a = 3.0 * L / 4.0
b = L - a

n_elements = [2, 4, 6, 8, 10]

# --- Hermite shape functions -------------

s, h = sp.symbols("s h", positive=True)

# w(s) = a0 + a1*s + a2*s^2 + a3*s^3
P_row = sp.Matrix([[1, s, s**2, s**3]])

# B matrix
# w(0) = a0
# w'(0) = a1
# w(h) = a0 + a1*h + a2*h^2 + a3*h^3
# w'(h) = a1 + 2*a2*h + 3*a3*h^2
B = sp.Matrix([
    [1,     0,      0,      0],
    [0,     1,      0,      0],
    [1,     h,      h**2,   h**3],
    [0,     1,      2*h,    3*h**2]
])


# Hermite shape functions
# N = PB^{-1}
N = sp.simplify(P_row * B.inv())

print("Hermite shape functions: \n      Derived from [N] = [P][B]^{-1}:")
for i in range(4):
    print(f"N{i+1} =", sp.simplify(N[0, i]))
    
    
# --- Stiffness Matrix --------------------

EI_sym = sp.symbols("EI", positive=True)

Ke_sym = sp.zeros(4,4)

for i in range(4):
    for j in range(4):
        Ni_dd = sp.diff(N[0, i], s, 2)
        Nj_dd = sp.diff(N[0,j], s, 2)
        
        # Element stiffness
        #   K_ij^e = intergal_0^h EI * N_i'' * N_j'' ds
        Ke_sym[i, j] = sp.simplify(sp.integrate(
            EI_sym * Ni_dd * Nj_dd, (s, 0, h)
        ))
        
print("\nElement stiffness matrix derived from:")
print("K_ij^e = integral_0^h EI * N_i'' * N_j'' ds")
sp.pprint(Ke_sym)

# Converting symbolic expressions into numerical functions
N_func = sp.lambdify((s, h), N, "numpy")
Ke_func = sp.lambdify((EI_sym, h), Ke_sym, "numpy")



# --- Exact solution ----------------------
def exact_deflection(x):
    x = np.asarray(x)
    w = np.zeros_like(x, dtype=float)
    
    left = x <= a
    right = x > a
    
    # for 0 <= x <= a
    w[left] = (
        P * b * x[left] / (6.0 * L * EI)
        * (L**2 - b**2 - x[left]**2)
    )
    
    # for a <= x <= L
    w[right] = (
        P * a * (L - x[right]) / (6.0 * L * EI)
        * (L**2 - a**2 - (L - x[right])**2)
    )
    return w


# --- FEM Solver --------------------------
def beam_fem(Ne):
    element_length = L / Ne
    nodes = np.linspace(0.0, L, Ne+1)
    
    # Each node has two DOFs:
    # w_i and theta_i
    n_dof = 2 * (Ne + 1)
    
    # Initialize stiffness matrix
    K = np.zeros((n_dof, n_dof))
    # Initialize load vector
    F = np.zeros(n_dof)
    
    # Sympy local stiffness
    Ke = np.array(Ke_func(EI, element_length), dtype=float)
    
    # Constructing stiffness matrix
    for e in range(Ne):
        element_dofs = np.array([
            2*e,
            2*e + 1,
            2*(e + 1),
            2*(e + 1) + 1
        ])
        
        for i in range(4):
            for j in range(4):
                K[element_dofs[i], element_dofs[j]] += Ke[i,j]

    # Constructing load vector
    if np.isclose(a, L):
        e = Ne - 1
    else:
        e = np.searchsorted(nodes, a, side="right") - 1
        
    x_left = nodes[e]
    x_local_load = a - x_left
    
    N_load = np.array(N_func(x_local_load, element_length), dtype=float).reshape(4)
    
    load_dofs = np.array([
        2*e,
        2*e+1,
        2*(e+1),
        2*(e+1)+1,
    ])
    
    F[load_dofs] += P * N_load
    
    
    # Boundary conditions
    constrained = [0, 2*Ne]
    free = [i for i in range(n_dof) if i not in constrained]
    
    W = np.zeros(n_dof)
    
    Kff = K[np.ix_(free, free)]
    Ff = F[free]
    
    W[free] = np.linalg.solve(Kff, Ff)
    
    return nodes, W, K, F

def fem_solution(x_plot, nodes, W):
    """
    w^e(x) = N1*w_i + N2*theta_i + N3*w_{i+1} + N4*theta_{i+1}
    """
    
    
    Ne = len(nodes) - 1
    element_length = nodes[1] - nodes[0]
    
    w_plot = np.zeros_like(x_plot, dtype=float)
    
    for k, x in enumerate(x_plot):
        if np.isclose(x, L):
            e = Ne - 1
        else:
            e = np.searchsorted(nodes, x, side="right") - 1
            
        x_left = nodes[e]
        x_local = x - x_left
        
        N_values = np.array(N_func(x_local, element_length), dtype=float).reshape(4)
        
        element_dofs = np.array([
            2*e,
            2*e + 1,
            2*(e+1),
            2*(e+1)+1
        ])
        
        w_plot[k] = N_values @ W[element_dofs]
        
    return w_plot
    

# Analysis
exact_at_load = exact_deflection(np.array([a]))[0]

fem_at_load = []

print("\nConvergence of deflection at load point x = 3L/4")
print("------------------------------------------------")
print(f"Exact deflection at load = {exact_at_load:.10f} m")
print()
print("Ne     FEM deflection      absolute error")
print("------------------------------------------")

for Ne in n_elements:
    nodes, W, K, F = beam_fem(Ne)
    
    w_load = fem_solution(np.array([a]), nodes, W)[0]
    fem_at_load.append(w_load)
    
    error = abs(w_load - exact_at_load)
    
    print(f"{Ne:<6d} {w_load:<18.10f} {error:.3e}")
    

# plotting

# convergence of FEM
plt.figure()
plt.plot(n_elements, fem_at_load, "o-", label="FEM")
plt.axhline(exact_at_load, linestyle="--", label="Exact")
plt.xlabel('Number of Elements')
plt.ylabel("Deflection at load point, w(3L/4) [m]")
plt.title("Convergence of FEM Deflection with Point Load")
plt.grid(True)
plt.legend()
plt.show()

# full solution @ Ne = 10
Ne = 10
nodes, W, K, F = beam_fem(Ne)

x_plot = np.linspace(0.0, L, 500)

w_fem = fem_solution(x_plot, nodes, W)
w_nodes = fem_solution(nodes, nodes, W)
w_exact = exact_deflection(x_plot)

plt.figure()
plt.plot(x_plot, w_fem, label="FEM, Ne = 10")
plt.plot(x_plot, w_exact, "--", label="Exact")
plt.plot(a, exact_at_load, "o", label="Load Point")
plt.plot(nodes, w_nodes, "o", color='black', label='Nodes')
plt.xlabel("x [m]")
plt.ylabel("Deflection w(x) [m]")
plt.title("Simply supported beam deflection under point load")
plt.grid(True)
plt.legend()
plt.show()


# === Problem 2 ===========================

# --- Givens ------------------------------
l = 1.0
chord = 0.3
thickness = 0.01

# Assumed aluminum
E_wing = 70e9

rho = 1.2
CL = 1.1
v_inf = 20.0
S = l * chord

# Rectangular MMOI
#   I = c*t^3 / 12
EI_wing = E_wing * chord * thickness**3 / 12.0

# Lift distribution:
#   L(x) = (1/2)*rho*v^2*S*CL*(1 - x^2/l^2)
q0 = 0.5 * rho * v_inf**2 * S * CL

print("\n" + "=" * 60)
print("PROBLEM 2 - Cantilever wing using Hermite FEM")
print(f"EI = {EI_wing:.6f} N*m^2")
print(f"q0 = {q0:.6f} N/m")
print("=" * 60)


# Distributed element load vector
x_left_sym = sp.symbols("x_left", real=True)
q0_sym = sp.symbols("q0", real=True)

# x position in each element:
#   x = x_left + s
x_global = x_left_sym + s

# distribution lift load
Lx_sym = q0_sym * (1 - x_global**2 / l**2)

# Element load vector:
#   F_i^e = integral_0^h L(x_left+s)*N_i(s) ds
Fe_sym = sp.zeros(4,1)
for i in range(4):
    Fe_sym[i, 0] = sp.simplify(
        sp.integrate(Lx_sym * N[0, i], (s, 0, h))
    )
    
print("\nElement load vector derived from:")
print("F_i^e = integral_0^h L(x_left+s)*N_i(s) ds")
sp.pprint(Fe_sym)

Fe_func = sp.lambdify((q0_sym, x_left_sym, h), Fe_sym, "numpy")

# FEM solver for cantilever wing
def wing_fem(Ne):
    element_length = l / Ne
    nodes = np.linspace(0.0, l, Ne + 1)
    
    n_dof = 2 * (Ne + 1)
    
    K = np.zeros((n_dof, n_dof))
    F = np.zeros(n_dof)
    
    Ke = np.array(Ke_func(EI_wing, element_length), dtype=float)
    
    for e in range(Ne):
        element_dofs = np.array([
            2*e,
            2*e + 1,
            2*(e+1),
            2*(e+1)+1
        ])
        for i in range(4):
            for j in range(4):
                K[element_dofs[i], element_dofs[j]] += Ke[i, j]
                
        x_left = nodes[e]
        
        Fe = np.array(
            Fe_func(q0, x_left, element_length),
            dtype=float
        ).reshape(4)
        
        F[element_dofs] += Fe
    
    constrained = [0, 1]
    free = [i for i in range(n_dof) if i not in constrained]
    
    W = np.zeros(n_dof)
    
    Kff = K[np.ix_(free, free)]
    Ff = F[free]

    W[free] = np.linalg.solve(Kff, Ff)
    
    return nodes, W, K, F

# FEM solution for cantilever wing
def wing_fem_solution(x_plot, nodes, W):
    Ne = len(nodes) - 1
    element_length = nodes[1] - nodes[0]
    
    w_plot = np.zeros_like(x_plot, dtype=float)
    
    for k, x_val in enumerate(x_plot):
        if np.isclose(x_val, l):
            e = Ne - 1
        else:
            e = np.searchsorted(nodes, x_val, side="right") - 1
        
        x_left = nodes[e]
        x_local = x_val - x_left
        
        N_values = np.array(
            N_func(x_local, element_length)
        ).reshape(4)
        
        element_dofs = np.array([
            2*e,
            2*e + 1,
            2*(e+1),
            2*(e+1) + 1
        ])
        
        w_plot[k] = N_values @ W[element_dofs]
    
    return w_plot

# =============================================================
# Rayleigh-Ritz comparison from Homework 2
# =============================================================

def solve_ritz_wing(n_basis):
    """
    Rayleigh-Ritz comparison using the same basis style from Homework 2:

        phi_k = x^(k+2)

    These satisfy the cantilever essential boundary conditions:
        w(0) = 0
        w'(0) = 0
    """

    x_rr = sp.symbols("x_rr")

    basis = [x_rr**(k + 2) for k in range(n_basis)]

    Lx_rr = q0 * (1 - x_rr**2 / l**2)

    K_rr = np.zeros((n_basis, n_basis))
    F_rr = np.zeros(n_basis)

    for i in range(n_basis):
        for j in range(n_basis):
            integrand_K = EI_wing * sp.diff(basis[i], x_rr, 2) * sp.diff(basis[j], x_rr, 2)
            K_rr[i, j] = float(sp.integrate(integrand_K, (x_rr, 0, l)))

        integrand_F = Lx_rr * basis[i]
        F_rr[i] = float(sp.integrate(integrand_F, (x_rr, 0, l)))

    coeffs = np.linalg.solve(K_rr, F_rr)

    w_expr = sum(coeffs[i] * basis[i] for i in range(n_basis))

    w_func = sp.lambdify(x_rr, w_expr, "numpy")

    tip = float(w_func(l))

    return tip, w_func


N_RITZ_BEST = 3
tip_ritz, w_ritz_func = solve_ritz_wing(N_RITZ_BEST)

# --- FEM convergence -----------------------------

fem_tip_deflections = []

print("\nTip deflection convergence")
print("--------------------------------------")
print("Ne       FEM tip deflection [mm]     Difference from Ritz [mm]")
print("--------------------------------------")

for Ne in n_elements:
    nodes, W, K, F = wing_fem(Ne)

    tip_fem = W[2*Ne]
    
    fem_tip_deflections.append(tip_fem)
    
    diff_from_ritz = tip_fem - tip_ritz
    
    print(f"{Ne:<6d}    {tip_fem*1000:<25.6f}   {diff_from_ritz*1000:.6f}")
    
    
# Plotting

plt.figure()
plt.plot(n_elements, np.array(fem_tip_deflections)*1000, "o-", color="black", label="Hermite FEM")
plt.axhline(tip_ritz*1000, linestyle="--", color="red", label=f"Reayleigh-Ritz, N={N_RITZ_BEST}")
plt.xlabel("Number of Elements")
plt.ylabel("Tip Deflection w(l) [mm]")
plt.title("Problem 2: Tip Deflection Convergence")
plt.grid(True)
plt.legend()
plt.ticklabel_format(axis="y", style="plain", useOffset=False)

tip_values_mm = np.array(fem_tip_deflections)*1000
tip_ref_mm = tip_ritz * 1000
ymin = min(np.min(tip_values_mm), tip_ref_mm)
ymax = max(np.max(tip_values_mm), tip_ref_mm)
margin = 0.05 * abs(ymax)
plt.ylim(ymin - margin, ymax + margin)

plt.show()

# Full solution
Ne = 10
nodes, W, K, F = wing_fem(Ne)

x_plot = np.linspace(0.0, l, 500)

w_fem = wing_fem_solution(x_plot, nodes, W)
w_ritz = w_ritz_func(x_plot)

plt.figure()
plt.plot(x_plot, w_fem*1000, label="Hermite FEM, Ne = 10")
plt.plot(x_plot, w_ritz*1000, "--", label=f"Rayleigh-Ritz, N={N_RITZ_BEST}")
plt.plot(l, W[2*Ne]*1000, "o", label="FEM tip")
plt.xlabel("x [m]")
plt.ylabel("Deflection w(x) [mm]")
plt.title("Problem 2: Cantilever Wing Deflection")
plt.grid(True)
plt.legend()
plt.show()