""" 8.12 Orbits In Class Problem """

# imports
import numpy as np
import Orbits_Functions as of
import matplotlib.pyplot as plt

# Constants
mu_e = 398600
mu_s = 1.327e20
mu_j = 1.267e17
rearth = 6378
rs2j = -778.6e6
rs2e = 149.6e6
rjupiter = 71490


# Transfer Orbit
rp_t = rs2e
ra_t = rs2j
ecct = (ra_t-rp_t)/(ra_t+rp_t)
at = (ra_t+rp_t)/2
ht = np.sqrt(rp_t*mu_s*(1+ecct*np.cos(0)))
vp_t = ht/rp_t
va_t = ht/ra_t
v_e = np.sqrt(mu_s/rs2e)
v_j = np.sqrt(mu_s/rs2j)
v_inf_e = vp_t - v_e
v_inf_j = va_t - v_j
h_t = vp_t @ rp_t
print(h_t)

# Jupiter orbit
rp_j = 200000+rjupiter
vp_j = np.sqrt(v_inf_j**2 + (2*mu_j)/rp_j)
ecc_j = 1 + (rp_j*v_inf_j**2)/mu_j


# Earth Orbit


