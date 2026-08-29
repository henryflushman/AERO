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
# ║   Course      :  AERO421 - Spacecraft Attitude Dynamics     ║
# ║                                     and Controls            ║
# ║   Assignment  :  Assignment 6                               ║
# ║   Date        :  May 15, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# From system
import numpy as np

# From directory:
from ADCS import (
    R3,
    Quaternion,
    is_valid_dcm,
    rotation_error_angle,
)

np.set_printoptions(precision=4, suppress=True)


# === Helper functions ===

# normalizing
def unit(v):
    v = np.asarray(v, dtype=float).reshape(3)
    n = np.linalg.norm(v)

    if n < 1e-12:
        raise ValueError("Cant normalize a zero vector")

    return v / n


def standardize_quaternion(q):
    """
    Normalize quaternion and choose the equivalent sign with q4 > 0.

    A quaternion has 4 elements:
        q = [q1, q2, q3, q4]

    where q4 is the scalar part.
    """
    q = np.asarray(q, dtype=float).reshape(4)
    q = q / np.linalg.norm(q)

    if q[3] < 0:
        q = -q

    return q


# printing
def print_vec(name, v):
    print(f"{name} =")
    print(np.asarray(v).reshape(-1, 1))
    print()


def print_mat(name, M):
    print(f"{name} =")
    print(np.asarray(M))
    print()


# triad
def triad_basis(v1, v2):
    """
    TRIAD basis construction.

    t1 = v1 / ||v1||
    t2 = (v1 x v2) / ||v1 x v2||
    t3 = t1 x t2

    The output matrix has the triad vectors as columns:
        C_bt = [t1_b, t2_b, t3_b]
    """
    t1 = unit(v1)
    t2 = unit(np.cross(v1, v2))
    t3 = np.cross(t1, t2)

    return np.column_stack((t1, t2, t3))


# attitude profile
def attitude_profile_matrix(body_vectors, ref_vectors, weights):
    """
    Attitude profile matrix.

    From the trace form of Wahba's problem:

        B = sum_k w_k s_bk s_ak^T

    Here:
        s_bk = measured vector in body frame
        s_ak = reference vector in inertial/ECI frame

    np.outer(b, r) computes b r^T.
    """
    B = np.zeros((3, 3))

    for w, b, r in zip(weights, body_vectors, ref_vectors):
        B += w * np.outer(unit(b), unit(r))

    return B


# davenport q
def davenport_q_method(body_vectors, ref_vectors, weights):
    """
    Davenport q-method.

    Build the 4x4 K matrix and solve:

        K q = lambda q

    The optimal quaternion is the eigenvector corresponding to the largest
    eigenvalue.
    """
    B = attitude_profile_matrix(body_vectors, ref_vectors, weights)

    S = B + B.T
    k22 = np.trace(B)

    k12 = np.array([
        B[1, 2] - B[2, 1],
        B[2, 0] - B[0, 2],
        B[0, 1] - B[1, 0],
    ])

    K = np.zeros((4, 4))
    K[:3, :3] = S - k22 * np.eye(3)
    K[:3, 3] = k12
    K[3, :3] = k12
    K[3, 3] = k22

    eigenvalues, eigenvectors = np.linalg.eigh(K)

    index_max = np.argmax(eigenvalues)
    lambda_max = eigenvalues[index_max]
    q_opt = eigenvectors[:, index_max]

    q_opt = standardize_quaternion(q_opt)

    return q_opt, lambda_max, B, S, k22, k12, K


def adjugate_3x3(A):
    """
    Classical adjoint / adjugate of a 3x3 matrix.

    This is used in QUEST.
    """
    A = np.asarray(A, dtype=float)

    a, b, c = A[0]
    d, e, f = A[1]
    g, h, i = A[2]

    return np.array([
        [e * i - f * h, c * h - b * i, b * f - c * e],
        [f * g - d * i, a * i - c * g, c * d - a * f],
        [d * h - e * g, b * g - a * h, a * e - b * d],
    ])


def quest_method(S, k12, k22, lambda_max):
    """
    QUEST reconstruction.

    q = 1 / sqrt(gamma^2 + x^T x) * [x; gamma]

    where:
        x = (alpha I + beta S + S^2) k12
        alpha = lambda^2 - k22^2 + tr(adj(S))
        beta = lambda - k22
        gamma = (lambda + k22) alpha - det(S)
    """
    I3 = np.eye(3)

    adjS = adjugate_3x3(S)
    tr_adjS = np.trace(adjS)
    detS = np.linalg.det(S)

    alpha = lambda_max**2 - k22**2 + tr_adjS
    beta = lambda_max - k22
    gamma = (lambda_max + k22) * alpha - detS

    # Eq. 25.17:
    # x = (alpha I + beta S + S^2) k12
    x = (alpha * I3 + beta * S + S @ S) @ k12

    q = np.append(x, gamma)
    q = standardize_quaternion(q)

    return q, alpha, beta, gamma, x, tr_adjS, detS


def dcm_from_q(q):
    return Quaternion(q).dcm


# === Analysis ===

def main():
    deg = np.pi / 180.0

    # === Given ===
    n_e_G = np.array([-1.0, 0.0, 0.0])
    n_s_G = np.array([0.0, 1.0, 0.0])
    rot_about_z = 45.0

    # === Part a ===
    C_bG_a = R3(rot_about_z * deg)

    # === Part b ===
    n_e_b = C_bG_a @ n_e_G
    n_s_b = C_bG_a @ n_s_G

    # === Part c, d, e ===
    C_bt = triad_basis(n_e_b, n_s_b)
    C_Gt = triad_basis(n_e_G, n_s_G)

    C_bG_triad_manual = C_bt @ C_Gt.T

    # === Part f ===
    body_vectors = [n_e_b, n_s_b]
    ref_vectors = [n_e_G, n_s_G]
    weights = [1.0, 1.0]

    q_dav, lambda_max, B, S, k22, k12, K = davenport_q_method(
        body_vectors,
        ref_vectors,
        weights,
    )

    C_bG_dav = dcm_from_q(q_dav)

    # quest reconstruction
    q_quest, alpha, beta, gamma, x, tr_adjS, detS = quest_method(
        S,
        k12,
        k22,
        lambda_max,
    )

    C_bG_quest = dcm_from_q(q_quest)

    # results
    print("=" * 70)
    print("PART (a): Known attitude")
    print("=" * 70)
    print_mat("C_bG from 45 deg rotation about z_G", C_bG_a)

    print("=" * 70)
    print("PART (b): Earth and sun vectors in body frame")
    print("=" * 70)
    print_vec("n_e_b", n_e_b)
    print_vec("n_s_b", n_s_b)

    print("=" * 70)
    print("PART (c): TRIAD basis vectors")
    print("=" * 70)
    print_mat("C_bt = [t1_b t2_b t3_b]", C_bt)
    print_mat("C_Gt = [t1_G t2_G t3_G]", C_Gt)

    print("=" * 70)
    print("PARTS (d)-(e): TRIAD attitude matrix")
    print("=" * 70)
    print_mat("C_bG = C_bt C_Gt.T", C_bG_triad_manual)

    print("=" * 70)
    print("PART (f): Davenport q-method")
    print("=" * 70)
    print_mat("B = sum w_k s_bk s_Gk.T", B)
    print_mat("S = B + B.T", S)
    print_vec("k12", k12)
    print(f"k22 = sigma = tr(B) = {k22:.12f}")
    print()
    print_mat("K matrix", K)
    print(f"lambda_max = {lambda_max:.12f}")
    print()
    print_vec("q_Davenport = [q1 q2 q3 q4]^T", q_dav)
    print_mat("C_bG from Davenport q-method", C_bG_dav)

    print("=" * 70)
    print("PART (f): QUEST reconstruction")
    print("=" * 70)
    print(f"tr(adj(S)) = {tr_adjS:.12f}")
    print(f"det(S)     = {detS:.12f}")
    print(f"alpha      = {alpha:.12f}")
    print(f"beta       = {beta:.12f}")
    print(f"gamma      = {gamma:.12f}")
    print()
    print_vec("x", x)
    print_vec("q_QUEST = [q1 q2 q3 q4]^T", q_quest)
    print_mat("C_bG from QUEST", C_bG_quest)

    print("=" * 70)
    print("Validation checks")
    print("=" * 70)
    print(f"is_valid_dcm(C_bG_a) = {is_valid_dcm(C_bG_a)}")

    print(
        "TRIAD vs part (a) error angle [deg]      = "
        f"{rotation_error_angle(C_bG_triad_manual, C_bG_a) / deg:.12e}"
    )

    print(
        "Davenport vs part (a) error angle [deg]  = "
        f"{rotation_error_angle(C_bG_dav, C_bG_a) / deg:.12e}"
    )

    print(
        "QUEST vs part (a) error angle [deg]      = "
        f"{rotation_error_angle(C_bG_quest, C_bG_a) / deg:.12e}"
    )


if __name__ == "__main__":
    main()