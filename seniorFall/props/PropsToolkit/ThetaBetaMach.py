


from __future__ import annotations

import numpy as np



def theta_from_beta_mach(beta: float, mach: float, gamma: float, degrees: bool = True) -> float:
    if degrees:
        beta = np.radians(beta)
    theta = np.arctan(2.0 * 1.0/np.tan(beta) * (mach**2 * np.sin(beta)**2 - 1) / (mach**2 * (gamma + np.cos(2.0 * beta) + 2.0)))
    if degrees:
        theta = np.degrees(theta)
    return theta


def normal_mach1_from_mach1(mach: float, beta: float, degrees: bool = True) -> float:
    if degrees:
        beta = np.radians(beta)
    normal_mach = mach * np.sin(beta)
    return normal_mach


def normal_mach2_from_normal_mach1(normal_mach1: float, gamma: float) -> float:
    normal_mach2 = np.sqrt((1.0 + ((gamma - 1.0) / 2.0) * normal_mach1**2.0) / (gamma * normal_mach1**2.0 - (gamma - 1.0) / 2.0))
    return normal_mach2

def mach2_from_normal_mach2(normal_mach: float, beta: float, degrees: bool = True) -> float:
    if degrees:
        beta = np.radians(beta)
    mach = normal_mach / np.sin(beta)
    return mach


class ObliqueShock:
    def __init__(self, mach: float, beta: float, gamma: float, degrees: bool = True):
        self.mach = mach
        self.beta = beta
        self.gamma = gamma
        self.degrees = degrees

    def theta(self) -> float:
        return theta_from_beta_mach(self.beta, self.mach, self.gamma, self.degrees)

    def normal_mach(self) -> float:
        return normal_mach1_from_mach1(self.mach, self.beta, self.degrees)

    def normal_mach2(self) -> float:
        return normal_mach2_from_normal_mach1(self.normal_mach(), self.gamma)

    def mach2_from_normal2(self) -> float:
        return mach2_from_normal_mach2(self.normal_mach2(), self.beta, self.degrees)