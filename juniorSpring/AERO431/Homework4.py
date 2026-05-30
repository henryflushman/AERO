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
# ║   Assignment  :  Homework 4                                 ║
# ║   Date        :  May 25, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# From system
import sys
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

# From directory
try:
    from libFRC_v1_1 import Laminate  # type: ignore
except ModuleNotFoundError:
    libFRC_path = Path(__file__).parent.parent.parent / "juniorWinter" / "AERO331"
    sys.path.insert(0, str(libFRC_path))
    from libFRC_v1_1 import Laminate  # type: ignore



# Global Function
def natural_frequencies_from_KM(K, M, nmodes=None):
    eigvals = eigh(K, M, eigvals_only=True)
    eighvals = np.real(eigvals)
    eighvals = eighvals[eigvals > 1e-7]
    omegas = np.sqrt(np.sort(eighvals))
    freqs_hz = omegas / (2.0 * np.pi)
    
    if nmodes is not None:
        freqs_hz = freqs_hz[:nmodes]
    
    return freqs_hz


# === Question 1 ===

def cantilever_bar_frequencies(E, rho, A, L, N=10, nmodes=3):
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    
    for i in range(1, N + 1):
        for j in range(1, N + 1):
            K[i - 1, j - 1] = E * A / L * (i * j) / (i + j - 1)
            M[i - 1, j - 1] = rho * A * L / (i + j + 1)
    
    return natural_frequencies_from_KM(K, M, nmodes)

def cantilever_beam_frequencies(E, rho, A, I, L, N=10, nmodes=3):
    K = np.zeros((N, N))
    M = np.zeros((N, N))
    
    for i in range(1, N + 1):
        p = i + 1
        
        for j in range(1, N + 1):
            q = j + 1
            
            K[i - 1, j - 1] = (
                E * I / L**3
                * (p * (p - 1) * q * (q - 1))
                / (p + q - 3)
            )
            
            M[i - 1, j - 1] = rho * A * L / (p + q + 1)
            
    return natural_frequencies_from_KM(K, M, nmodes)
    
    
# === Question 2 and 3 ===

def isotropic_plate_D(E, nu, h):
    D0 = E * h**3 / (12.0 * (1.0 - nu**2))
    
    D = D0 * np.array([
        [1.0, nu, 0.0],
        [nu, 1.0, 0.0],
        [0.0, 0.0, (1.0 - nu) / 2.0]
    ])
    
    return D


def make_sine_basis(mmax, nmax):
    return [(m, n) for m in range(1, mmax + 1) for n in range(1, nmax + 1)]

def simply_supported_plate_frequencies(
    D,
    rho_areal,
    a,
    b,
    mmax=4,
    nmax=4,
    nmodes=6,
    nq=None,
):
    basis = make_sine_basis(mmax, nmax)
    nbasis = len(basis)
    
    if nq is None:
        nq = max(80, 8 * max(mmax, nmax) + 10)
    
    gx, wx = np.polynomial.legendre.leggauss(nq)
    gy, wy = np.polynomial.legendre.leggauss(nq)
    
    x = 0.5 * a * (gx + 1.0)
    y = 0.5 * b * (gy + 1.0)
    
    wx = 0.5 * a * wx
    wy = 0.5 * b * wy
    
    X, Y = np.meshgrid(x, y, indexing='ij')
    W2 = np.outer(wx, wy)
    
    phi = []
    kappa = []
    
    for m, n in basis:
        sx = np.sin(m * np.pi * X / a)
        cx = np.cos(m * np.pi * X / a)
        
        sy = np.sin(n * np.pi * Y / b)
        cy = np.cos(n * np.pi * Y / b)
        
        ph = sx * sy
        
        w_xx = -((m * np.pi / a) ** 2) * ph
        w_yy = -((n * np.pi / b) ** 2) * ph
        two_w_xy = 2.0 * (m * np.pi / a) * (n * np.pi / b) * cx * cy
        
        phi.append(ph)
        kappa.append(np.stack([w_xx, w_yy, two_w_xy], axis=0))
    
    K = np.zeros((nbasis, nbasis))
    M = np.zeros((nbasis, nbasis))
    
    for i in range(nbasis):
        for j in range(i, nbasis):
            D_kappa_j = np.einsum("ab,bxy->axy", D, kappa[j])
            stiffness_density = np.einsum("axy,axy->xy", kappa[i], D_kappa_j)
            
            Kij = np.sum(W2 * stiffness_density)
            Mij = np.sum(W2 * rho_areal * phi[i] * phi[j])
            
            K[i, j] = Kij
            K[j, i] = Kij
            
            M[i, j] = Mij
            M[j, i] = Mij
            
        
    return natural_frequencies_from_KM(K, M, nmodes)

def symmetric_layup(theta_deg):
    return [theta_deg, 0.0, 90.0, 90.0, 0.0, theta_deg]

def laminate_D_from_libFRC(theta_deg, total_thickness):
    laminate = Laminate(theta=symmetric_layup(theta_deg), t=total_thickness)
    return laminate.D


# === Main Execution ===

if __name__ == "__main__":
    
    # Material properties
    E_al = 70e9
    rho_al = 2700.0
    nu_al = 0.30
    
    # Q1 geometry
    L = 2.0
    c = 0.10
    t = 0.01
    
    A = c * t
    
    Iyy = c * t**3 / 12.0
    Izz = t * c**3 / 12.0
    
    f_bar = cantilever_bar_frequencies(E_al, rho_al, A, L, N=10, nmodes=3)
    f_xz = cantilever_beam_frequencies(E_al, rho_al, A, Iyy, L, N=10, nmodes=3)
    f_xy = cantilever_beam_frequencies(E_al, rho_al, A, Izz, L, N=10, nmodes=3)
    
    print("Q1 first three frequencies [Hz]")
    print("   Axial bar:      ", f_bar)
    print("   Bending (xz):   ", f_xz)
    print("   Bending (xy):   ", f_xy)
    
    all_freqs = []
    
    for label, values in [
        ("bar axial", f_bar),
        ("beam xz", f_xz),
        ("beam xy", f_xy),
    ]:
        for mode_number, value in enumerate(values, start=1):
            all_freqs.append((value, label, mode_number))
            
    all_freqs.sort(key=lambda item: item[0])
    for value, label, mode_number in all_freqs:
        print(f"   {value:9.2f} Hz    {label}, mode {mode_number}")
        
    
    # Q1 OBSERVATIONS
    #   The lowest frequencies are the beam bending modes, then the axial bar modes.
    # Bending in the xz plane is easier than the xy plane because t is much smaller than c.    
        
    
    # Q2 isotropic aluminum plate
    a = 1.0
    b = 0.25
    h = 0.05
    
    D_al_plate = isotropic_plate_D(E_al, nu_al, h)
    rho_areal_al = rho_al * h
    
    f_q2 = simply_supported_plate_frequencies(
        D_al_plate,
        rho_areal_al,
        a,
        b,
        mmax=4,
        nmax=4,
        nmodes=1,
    )[0]
    
    print(f"\nQ2 aluminum plate fundamental freqeuncy: {f_q2:.2f} Hz")
    
    # Q3 composite plate
    rho_comp = 1600.0
    rho_areal_comp = rho_comp * h
    
    mmax = 5
    nmax = 5
    
    # i. theta = 30
    theta = 30.0
    
    D_comp_30 = laminate_D_from_libFRC(theta, h)
    
    f_q3_30 = simply_supported_plate_frequencies(
        D_comp_30,
        rho_areal_comp,
        a,
        b,
        mmax=mmax,
        nmax=nmax,
        nmodes=1,
    )[0]
    
    print(f"\nQ3(i) composite plate f1 at theta = {theta:.0f} deg: {f_q3_30:.2f} Hz")
    
    # ii. sweep theta from -90 to 90 in 15 degree increments
    thetas = np.arange(-90.0, 90.0 + 0.1, 15.0)
    f_theta = []
    
    for theta_i in thetas:
        D_comp = laminate_D_from_libFRC(theta_i, h)
        
        f1 = simply_supported_plate_frequencies(
            D_comp,
            rho_areal_comp,
            a,
            b,
            mmax=mmax,
            nmax=nmax,
            nmodes=1,
        )[0]
        
        f_theta.append(f1)
        
    f_theta = np.array(f_theta)
    
    print("\nQ3(ii) theta sweep")
    for theta_i, f_i in zip(thetas, f_theta):
        print(f"   theta = {theta_i:6.1f} deg: f1 = {f_i:10.2f} Hz")
        
    plt.figure()
    plt.plot(thetas, f_theta, marker='o')
    plt.xlabel(r"$\theta$ [deg]")
    plt.ylabel("fundamental frequency [Hz]")
    plt.title(r"Simply supported $[\theta, 0, 90]_s$ Carbon/Epoxy plate")
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    
    # Q3 OBSERVATIONS:
    #   The fundamental frequency is symmetric about theta=0. It is lowest at theta=0,
    # at about 1202.1 Hz, and increases as the magnitude of theta increases, reaching a
    # maximum at theta=90 or theta=-90 at about 2998 Hz. This changes the laminate structure