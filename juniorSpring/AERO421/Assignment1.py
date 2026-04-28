"""
Gravity-Gradient Torque Calculation
====================================
Problem:
    Spacecraft with principal moments of inertia Ix=100, Iy=120, Iz=80 kg·m².
    Attitude described by 3-2-1 (yaw-pitch-roll) Euler angles:
        φ (roll) = θ (pitch) = ψ (yaw) = π/4 rad
    Orbital position in ECI:  R_o = [0, 0, R0]ᵀ  with R0 = 7000 km.
    μ = 3.986 × 10⁵ km³/s²

Find: Gravity-gradient torque expressed in spacecraft body coordinates.

Theory:
    The gravity-gradient torque is:
        τ_gg = (3μ / R₀³) · ê_r_B × (I · ê_r_B)

    where ê_r_B is the unit nadir vector (Earth-centre direction) expressed
    in body frame:
        ê_r_B = C_BG · ê_r_G,   ê_r_G = -R̂_o (points from s/c toward Earth)

    C_BG is the passive 3-2-1 DCM mapping ECI → Body frame.
"""

import sys
import os
import numpy as np
import analysis.helpers.ADCS as adcs

# ═══════════════════════════════════════════════════════════════════════════════
# 1. Given Parameters
# ═══════════════════════════════════════════════════════════════════════════════

# Gravitational parameter [m³/s²]
mu_m = 3.986e14          # m³/s²

# Orbital position in ECI [m]
R0 = 7000.0e3            # m
R0_G = np.array([0.0, 0.0, R0])

# Euler angles (3-2-1 sequence):
phi   = np.pi / 4        # roll
theta = np.pi / 4        # pitch
psi   = np.pi / 4        # yaw

# Principal moments of inertia [kg·m²]
Ix, Iy, Iz = 100.0, 120.0, 80.0
J_b = np.diag([Ix, Iy, Iz])

# ═══════════════════════════════════════════════════════════════════════════════
# 2. Build DCM  C_BG  (ECI → Body),  3-2-1 sequence
# ═══════════════════════════════════════════════════════════════════════════════

C_bG = adcs.DCM([psi, theta, phi], sequence=(3, 2, 1)).matrix   # 3×3 ndarray

# ═══════════════════════════════════════════════════════════════════════════════
# 3. Convert Vectors to Body Frame
# ═══════════════════════════════════════════════════════════════════════════════

R0_b = C_bG @ R0_G

# ═══════════════════════════════════════════════════════════════════════════════
# 4. Gravity-Gradient Torque
# ═══════════════════════════════════════════════════════════════════════════════

JR0_b = J_b @ R0_b

scalar_prefactor = 3 * mu_m / (R0 ** 5)

T_g = scalar_prefactor * np.cross(R0_b, JR0_b)

# ═══════════════════════════════════════════════════════════════════════════════
# 5. Print Results
# ═══════════════════════════════════════════════════════════════════════════════

sep = "=" * 62

print(sep)
print("  Gravity-Gradient Torque  –  Results")
print(sep)

print("\n── Attitude (3-2-1 Euler angles) ──────────────────────────")
print(f"   Roll  φ = {np.degrees(phi):.2f}°  ({phi:.6f} rad)")
print(f"   Pitch θ = {np.degrees(theta):.2f}°  ({theta:.6f} rad)")
print(f"   Yaw   ψ = {np.degrees(psi):.2f}°  ({psi:.6f} rad)")

print("\n── Direction Cosine Matrix  C_BG (ECI → Body) ─────────────")
for row in C_bG:
    print(f"   [{row[0]:+.6f}   {row[1]:+.6f}   {row[2]:+.6f}]")

print("\n── Orbital Position ───────────────────────────────────────")
print(f"   In ECI frame  : {R0_G}")
print(f"   In Body frame : [{R0_b[0]:+.6f},  {R0_b[1]:+.6f},  {R0_b[2]:+.6f}]")

print("\n── Gravity-Gradient Prefactor ─────────────────────────────")
print(f"   3μ/R₀⁵ = {scalar_prefactor:.6e}  1/(m²·s²)")

print("\n── Gravity-Gradient Torque (Body Frame) ───────────────────")
print(f"   τ_gg = [{T_g[0]:+.6e},  {T_g[1]:+.6e},  {T_g[2]:+.6e}]  N·m")
print(f"   |τ_gg| = {np.linalg.norm(T_g):.6e}  N·m")

print(f"\n   τ_gg_x = {T_g[0]:+.4e}  N·m")
print(f"   τ_gg_y = {T_g[1]:+.4e}  N·m")
print(f"   τ_gg_z = {T_g[2]:+.4e}  N·m")

print("\n" + sep)


