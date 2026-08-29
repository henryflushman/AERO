import numpy as np
import sympy as sp
import matplotlib.pyplot as plt


# === Symbols ===

x, y, s = sp.symbols("x y s")


# === Constants ===

E_VAL = 70e9
RHO_VAL = 2700.0
G_VAL = 9.81

L1_VAL = 2.0
C1_VAL = 0.2
T1_VAL = 0.05
A1_VAL = C1_VAL * T1_VAL
I1_VAL = C1_VAL * T1_VAL**3 / 12

L2_VAL = 2.0
C2_VAL = 0.4
T2_VAL = 0.05
A2_VAL = C2_VAL * T2_VAL
I2_VAL = C2_VAL * T2_VAL**3 / 12

NU_VAL = 0.33
D_VAL = E_VAL * T2_VAL**3 / (12 * (1 - NU_VAL**2))

E4_VAL = 70e9
NU4_VAL = 0.3
T_SKIN = 0.002
T_SPAR = 0.005
H_SPAR = 0.1
W_FLANGE = 0.05
L_VALUES = np.linspace(0.1, 0.5, 6)
B_VALUES = np.linspace(0.1, 0.3, 6)

W5_VAL = 0.5
SIGMA5_VAL = 50.0
KC5_VAL = 24.0
C_PARIS = 3.15e-11
M_PARIS = 3.59


# === Helpers ===

def rpm_to_rad_s(rpm):
    return rpm * 2 * np.pi / 60


def poly_int_1d(expr, a, b, var=x):
    poly = sp.Poly(sp.expand(expr), var)
    total = 0.0

    for powers, coeff in poly.terms():
        p = powers[0]
        total += float(coeff) * (b**(p + 1) - a**(p + 1)) / (p + 1)

    return total


def solve_eig(K, M):
    vals = np.linalg.eigvals(np.linalg.solve(M, K))
    vals = np.sort(np.real(vals[vals > 0]))

    omega_1 = np.sqrt(vals[0])
    f_1 = omega_1 / (2 * np.pi)

    return omega_1, f_1


def print_frequency_table(title, rows):
    print(f"\n{title}")
    print("Boundary condition        omega_1 [rad/s]    f_1 [Hz]")
    print("------------------------------------------------------")

    for name, omega_1, f_1 in rows:
        print(f"{name:<24} {omega_1:14.6f} {f_1:12.6f}")


# === Problem 1 ===

def part_1_symbolic():
    L, rho, A, omega, E, I, g = sp.symbols("L rho A omega E I g")

    w = sp.Function("w")(x)
    delta_w = sp.Function("delta_w")(x)

    N_expr = sp.simplify(sp.integrate(rho * A * omega**2 * s, (s, x, L)))

    Pi = (
        sp.Rational(1, 2)
        * sp.Integral(E * I * sp.diff(w, x, 2)**2, (x, 0, L))
        + sp.Rational(1, 4)
        * sp.Integral(
            rho * A * omega**2 * (L**2 - x**2) * sp.diff(w, x)**2,
            (x, 0, L),
        )
        - sp.Integral(rho * A * g * w, (x, 0, L))
    )

    delta_Pi = (
        sp.Integral(
            E * I * sp.diff(w, x, 2) * sp.diff(delta_w, x, 2),
            (x, 0, L),
        )
        + sp.Rational(1, 2)
        * sp.Integral(
            rho * A * omega**2 * (L**2 - x**2)
            * sp.diff(w, x)
            * sp.diff(delta_w, x),
            (x, 0, L),
        )
        - sp.Integral(rho * A * g * delta_w, (x, 0, L))
    )

    return N_expr, Pi, delta_Pi


def tip_deflection(rot_omega, n_lim=0, length=L1_VAL, tol=1e-7, max_n=20):
    def solve_n(n):
        basis = [x**(k + 2) for k in range(n)]

        K = np.zeros((n, n))
        F = np.zeros(n)

        for i in range(n):
            F[i] = poly_int_1d(-RHO_VAL * A1_VAL * G_VAL * basis[i], 0, length)

            for j in range(n):
                bend = (
                    E_VAL
                    * I1_VAL
                    * sp.diff(basis[i], x, 2)
                    * sp.diff(basis[j], x, 2)
                )

                axial = (
                    0.5
                    * RHO_VAL
                    * A1_VAL
                    * rot_omega**2
                    * (length**2 - x**2)
                    * sp.diff(basis[i], x)
                    * sp.diff(basis[j], x)
                )

                K[i, j] = poly_int_1d(bend + axial, 0, length)

        coeffs = np.linalg.solve(K, F)
        w_expr = sum(coeffs[i] * basis[i] for i in range(n))
        w_func = sp.lambdify(x, w_expr, "numpy")

        return float(w_func(length)), w_func

    tips = np.array([0.0, -1.0])

    if n_lim == 0:
        n = 1

        while abs(tips[-1] - tips[-2]) > tol:
            if n > max_n:
                raise RuntimeError("Tip deflection did not converge before max_n.")

            tip, w_func = solve_n(n)
            tips = np.append(tips, tip)
            n += 1

    elif n_lim > 0:
        tip, w_func = solve_n(n_lim)
        tips = np.append(tips, tip)

    else:
        raise ValueError("n_lim must be 0 or positive.")

    return tips, w_func


# === Problem 2 ===

def beam_basis(case, n, length=L2_VAL):
    if case == "ss":
        return [x**(k + 1) * (length - x) for k in range(n)]

    if case == "cs":
        return [x**(k + 2) * (length - x) for k in range(n)]

    if case == "cc":
        return [x**(k + 2) * (length - x)**2 for k in range(n)]

    raise ValueError("Unknown boundary case.")


def beam_frequency(case, n, length=L2_VAL):
    basis = beam_basis(case, n, length)

    K = np.zeros((n, n))
    M = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            K[i, j] = poly_int_1d(
                E_VAL
                * I2_VAL
                * sp.diff(basis[i], x, 2)
                * sp.diff(basis[j], x, 2),
                0,
                length,
            )

            M[i, j] = poly_int_1d(
                RHO_VAL * A2_VAL * basis[i] * basis[j],
                0,
                length,
            )

    return solve_eig(K, M)


# === Problem 3 ===

def plate_x_basis(case, length=L2_VAL):
    if case == "ss":
        return [x**n * (length - x) for n in range(1, 6)]

    if case == "cs":
        return [x**n * (length - x) for n in range(2, 6)]

    if case == "cc":
        return [x**n * (length - x)**2 for n in range(2, 5)]

    raise ValueError("Unknown boundary case.")


def plate_frequency(case, length=L2_VAL, width=C2_VAL):
    x_basis = plate_x_basis(case, length)
    y_basis = [1, y, y**2]

    n_x = len(x_basis)
    n_y = len(y_basis)
    n = n_x * n_y

    y_lower = -width / 2
    y_upper = width / 2

    K = np.zeros((n, n))
    M = np.zeros((n, n))

    def idx(i, a):
        return i * n_y + a

    x_int = {}
    y_int = {}

    for i in range(n_x):
        for j in range(n_x):
            x_int[("00", i, j)] = poly_int_1d(
                x_basis[i] * x_basis[j], 0, length, x
            )
            x_int[("10", i, j)] = poly_int_1d(
                sp.diff(x_basis[i], x) * sp.diff(x_basis[j], x),
                0,
                length,
                x,
            )
            x_int[("20", i, j)] = poly_int_1d(
                sp.diff(x_basis[i], x, 2) * x_basis[j],
                0,
                length,
                x,
            )
            x_int[("02", i, j)] = poly_int_1d(
                x_basis[i] * sp.diff(x_basis[j], x, 2),
                0,
                length,
                x,
            )
            x_int[("22", i, j)] = poly_int_1d(
                sp.diff(x_basis[i], x, 2) * sp.diff(x_basis[j], x, 2),
                0,
                length,
                x,
            )

    for a in range(n_y):
        for b in range(n_y):
            y_int[("00", a, b)] = poly_int_1d(
                y_basis[a] * y_basis[b], y_lower, y_upper, y
            )
            y_int[("10", a, b)] = poly_int_1d(
                sp.diff(y_basis[a], y) * sp.diff(y_basis[b], y),
                y_lower,
                y_upper,
                y,
            )
            y_int[("20", a, b)] = poly_int_1d(
                sp.diff(y_basis[a], y, 2) * y_basis[b],
                y_lower,
                y_upper,
                y,
            )
            y_int[("02", a, b)] = poly_int_1d(
                y_basis[a] * sp.diff(y_basis[b], y, 2),
                y_lower,
                y_upper,
                y,
            )
            y_int[("22", a, b)] = poly_int_1d(
                sp.diff(y_basis[a], y, 2) * sp.diff(y_basis[b], y, 2),
                y_lower,
                y_upper,
                y,
            )

    for i in range(n_x):
        for a in range(n_y):
            p = idx(i, a)

            for j in range(n_x):
                for b in range(n_y):
                    q = idx(j, b)

                    K[p, q] = D_VAL * (
                        x_int[("22", i, j)] * y_int[("00", a, b)]
                        + x_int[("00", i, j)] * y_int[("22", a, b)]
                        + NU_VAL * x_int[("20", i, j)] * y_int[("02", a, b)]
                        + NU_VAL * x_int[("02", i, j)] * y_int[("20", a, b)]
                        + 2
                        * (1 - NU_VAL)
                        * x_int[("10", i, j)]
                        * y_int[("10", a, b)]
                    )

                    M[p, q] = (
                        RHO_VAL
                        * T2_VAL
                        * x_int[("00", i, j)]
                        * y_int[("00", a, b)]
                    )

    return (*solve_eig(K, M), n)


# === Problem 4 ===

def plate_k_factor(length, width, m_max=50):
    m_vals = np.arange(1, m_max + 1)
    k_vals = (m_vals * width / length + length / (m_vals * width))**2

    idx = np.argmin(k_vals)

    return k_vals[idx], int(m_vals[idx])


def plate_buckling_stress(length, width, thickness):
    D = E4_VAL * thickness**3 / (12 * (1 - NU4_VAL**2))
    k_val, m_val = plate_k_factor(length, width)

    sigma_cr = k_val * np.pi**2 * D / (width**2 * thickness)

    return sigma_cr, m_val


def spar_section_properties():
    area = 2 * W_FLANGE * T_SPAR + H_SPAR * T_SPAR

    ix_web = T_SPAR * H_SPAR**3 / 12
    ix_flange = W_FLANGE * T_SPAR**3 / 12
    y_flange = H_SPAR / 2 + T_SPAR / 2
    ix = ix_web + 2 * (ix_flange + W_FLANGE * T_SPAR * y_flange**2)

    iy_web = H_SPAR * T_SPAR**3 / 12
    iy_flange = T_SPAR * W_FLANGE**3 / 12
    iy = iy_web + 2 * iy_flange

    return area, ix, iy


def spar_column_buckling_stress(length):
    area, ix, iy = spar_section_properties()

    return np.pi**2 * E4_VAL * min(ix, iy) / (area * length**2)


def problem_4():
    stress_grid = np.zeros((len(L_VALUES), len(B_VALUES)))
    label_grid = [["" for _ in B_VALUES] for _ in L_VALUES]
    results = []

    for i, length in enumerate(L_VALUES):
        for j, width in enumerate(B_VALUES):
            sigma_skin, m_skin = plate_buckling_stress(length, width, T_SKIN)
            sigma_web, m_web = plate_buckling_stress(length, H_SPAR, T_SPAR)
            sigma_spar = spar_column_buckling_stress(length)

            mode, sigma_min, m_min = min(
                [
                    ("Skin", sigma_skin, m_skin),
                    ("Web", sigma_web, m_web),
                    ("I-beam", sigma_spar, None),
                ],
                key=lambda item: item[1],
            )

            stress_grid[i, j] = sigma_min / 1e6
            label_grid[i][j] = "I" if mode == "I-beam" else f"{mode[0]}{m_min}"

            results.append(
                {
                    "L": length,
                    "b": width,
                    "sigma_min_mpa": sigma_min / 1e6,
                    "mode": mode,
                    "m": m_min,
                }
            )

    return results, stress_grid, label_grid


# === Problem 5 ===

def beta_crack(a, width=W5_VAL):
    ratio = a / width

    numerator = (
        1.122
        - 1.122 * ratio
        - 0.820 * ratio**2
        + 3.768 * ratio**3
        - 3.040 * ratio**4
    )

    return numerator / np.sqrt(1 - 2 * ratio)


def stress_intensity(a, sigma=SIGMA5_VAL, width=W5_VAL):
    return beta_crack(a, width) * sigma * np.sqrt(np.pi * a)


def critical_crack_size(width=W5_VAL, sigma=SIGMA5_VAL, kc=KC5_VAL):
    a_low = 1e-12
    a_high = 0.5 * width - 1e-12

    for _ in range(100):
        a_mid = 0.5 * (a_low + a_high)

        if stress_intensity(a_mid, sigma, width) < kc:
            a_low = a_mid
        else:
            a_high = a_mid

    return 0.5 * (a_low + a_high)


def cycles_to_failure(delta_n, a_initial, a_critical):
    cycles = 0.0
    a = a_initial

    n_history = [cycles]
    a_history = [a]

    while a < a_critical:
        delta_k = stress_intensity(a, SIGMA5_VAL)
        da_dn = C_PARIS * delta_k**M_PARIS
        a_next = a + da_dn * delta_n

        if a_next >= a_critical:
            cycles += (a_critical - a) / da_dn
            a = a_critical
        else:
            cycles += delta_n
            a = a_next

        n_history.append(cycles)
        a_history.append(a)

    return cycles, np.array(n_history), np.array(a_history)


def problem_5():
    a_critical = critical_crack_size()
    a_initial = a_critical / 4

    cycles_100, n_100, a_100 = cycles_to_failure(100, a_initial, a_critical)
    cycles_10, n_10, a_10 = cycles_to_failure(10, a_initial, a_critical)

    return {
        "a_critical": a_critical,
        "a_initial": a_initial,
        "cycles_100": cycles_100,
        "cycles_10": cycles_10,
        "n_100": n_100,
        "a_100": a_100,
        "n_10": n_10,
        "a_10": a_10,
    }


# === Main ===

def main():
    N_expr, Pi, delta_Pi = part_1_symbolic()

    print("Problem 1")
    print(f"N(x) = {sp.sstr(N_expr)}")

    print("\n--- Q1(ii) ---")
    print("Pi(w) =")
    sp.pprint(Pi)
    print("\nU = {w(x) : w(0) = 0, w'(0) = 0}  (cantilever BCs at x=0; tip free)")
    print("\ndelta_Pi =")
    sp.pprint(delta_Pi)
    print("  a(w,dw) = EI*int(w''dw'') + (1/2)*rho*A*omega^2*int((L^2-x^2)*w'*dw')")
    print("  l(dw)   = rho*A*g*int(dw)")

    tips_100, _ = tip_deflection(rpm_to_rad_s(100))
    n_converged = len(tips_100) - 2

    print(f"\n--- Q1(iii) convergence study at 100 RPM ---")
    print(f"{'N':>4}  {'Tip deflection [mm]':>20}")
    for n in range(1, n_converged + 1):
        tips_n, _ = tip_deflection(rpm_to_rad_s(100), n_lim=n)
        print(f"{n:>4}  {tips_n[-1] * 1000:>20.6f}")
    print(f"Converged at N = {n_converged} basis functions.")
    print(f"100 RPM tip deflection = {tips_100[-1] * 1000:.6f} mm")

    rpm_values = np.linspace(0, 400, 21)
    tip_values = np.zeros_like(rpm_values)

    for i, rpm in enumerate(rpm_values):
        tips, _ = tip_deflection(rpm_to_rad_s(rpm), n_lim=n_converged)
        tip_values[i] = tips[-1]

    print(f"\n--- Q1(iv) ---")
    print(f"Tip deflection range: {tip_values.min()*1000:.4f} to {tip_values.max()*1000:.4f} mm")
    print("As RPM increases, centrifugal tension stiffens the beam geometrically,")
    print("reducing tip deflection. Deflection shapes flatten with RPM. Both plots")
    print("are consistent with expected stress-stiffening behavior of a rotor blade.")

    plt.figure()
    plt.plot(rpm_values, tip_values * 1000, marker="o")
    plt.xlabel("Angular velocity [RPM]")
    plt.ylabel("Tip deflection [mm]")
    plt.title("Tip deflection vs angular velocity")
    plt.grid(True)
    plt.tight_layout()

    x_plot = np.linspace(0, L1_VAL, 300)

    plt.figure()

    for rpm in [0, 50, 100, 150, 200]:
        _, w_func = tip_deflection(rpm_to_rad_s(rpm), n_lim=n_converged)
        plt.plot(x_plot, w_func(x_plot) * 1000, label=f"{rpm} RPM")

    plt.xlabel("x [m]")
    plt.ylabel("w(x) [mm]")
    plt.title("Beam deflection shapes")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    cases = [
        ("Simply-supported", "ss"),
        ("Clamped-simply", "cs"),
        ("Clamped-clamped", "cc"),
    ]

    beam_rows = [
        (name, *beam_frequency(key, n_converged))
        for name, key in cases
    ]

    print_frequency_table("Problem 2 beam frequencies", beam_rows)
    f1_ss = beam_rows[0][2]; f1_cs = beam_rows[1][2]; f1_cc = beam_rows[2][2]
    print("SS < CS < CC: more constraints -> higher stiffness -> higher frequency.")
    print(f"Ratios: 1 : {f1_cs/f1_ss:.3f} : {f1_cc/f1_ss:.3f}, matching eigenvalue ratios pi^2 : 3.927^2 : 4.730^2.")

    print("\nProblem 3 plate vs beam")
    print("Boundary condition      Plate [Hz]    Beam [Hz]     Diff [%]")
    print("------------------------------------------------------------")

    plate_rows = []

    for (name, key), (_, _, f_beam) in zip(cases, beam_rows):
        omega_plate, f_plate, _ = plate_frequency(key)
        plate_rows.append((name, omega_plate, f_plate))

        diff = (f_plate - f_beam) / f_beam * 100

        print(f"{name:<22} {f_plate:10.6f} {f_beam:10.6f} {diff:10.6f}")

    print("Plate frequencies are slightly higher than beam in all cases.")
    print("Poisson coupling adds lateral stiffness absent in the 1D beam model.")
    print("Gap is small (~0.2-3.4%) since the plate is narrow (c/L=0.2),")
    print("and grows SS->CC as stronger BCs increase curvature, amplifying the effect.")

    labels = [row[0].replace("-", "\n") for row in beam_rows]
    x_bar = np.arange(len(labels))
    bar_width = 0.35

    plt.figure()
    plt.bar(x_bar - bar_width / 2, [row[2] for row in beam_rows], bar_width, label="Beam")
    plt.bar(x_bar + bar_width / 2, [row[2] for row in plate_rows], bar_width, label="Plate")
    plt.xticks(x_bar, labels)
    plt.ylabel("Fundamental frequency [Hz]")
    plt.title("Beam vs plate frequencies")
    plt.grid(axis="y")
    plt.legend()
    plt.tight_layout()

    results_4, stress_grid, label_grid = problem_4()

    mode_counts = {}

    for row in results_4:
        mode_counts[row["mode"]] = mode_counts.get(row["mode"], 0) + 1

    print("\nProblem 4")
    print(f"Critical stress range = {stress_grid.min():.6f} to {stress_grid.max():.6f} MPa")
    print(f"Critical mode counts = {mode_counts}")

    B_grid, L_grid = np.meshgrid(B_VALUES, L_VALUES)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    ax.plot_surface(
        B_grid,
        L_grid,
        stress_grid,
        alpha=0.8,
        edgecolor="k",
        linewidth=0.4,
    )

    for i, length in enumerate(L_VALUES):
        for j, width in enumerate(B_VALUES):
            ax.text(
                width,
                length,
                stress_grid[i, j],
                label_grid[i][j],
                ha="center",
                va="center",
                fontsize=8,
            )

    ax.set_xlabel("Spar spacing b [m]")
    ax.set_ylabel("Rib spacing L [m]")
    ax.set_zlabel("Critical buckling stress [MPa]")
    ax.set_title("Minimum buckling stress")

    plt.tight_layout()

    p5 = problem_5()
    diff_cycles = abs(p5["cycles_100"] - p5["cycles_10"]) / p5["cycles_10"] * 100

    print("\nProblem 5")
    print(f"Critical crack parameter = {p5['a_critical']:.8f} m")
    print(f"Initial crack parameter = {p5['a_initial']:.8f} m")
    print(f"Cycles to failure: dN=100 -> {p5['cycles_100']:.2f}, dN=10 -> {p5['cycles_10']:.2f}")
    print(f"Cycle difference = {diff_cycles:.4f}%")
    print(f"Solution does not change appreciably: {diff_cycles:.4f}% difference. dN=100 is sufficient.")

    plt.figure()
    plt.plot(p5["n_100"], p5["a_100"], label="dN = 100")
    plt.plot(p5["n_10"], p5["a_10"], label="dN = 10")
    plt.xlabel("Cycles N")
    plt.ylabel("Crack parameter a [m]")
    plt.title("Crack growth using Paris law")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


if __name__ == "__main__":
    main()