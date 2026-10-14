"""Lunar-surface user transmit power needed to reach an Earth-Moon L1 spacecraft (Ka-band RF).

Screening-level link budget at the longest line-of-sight range. Prints the
transmitter power a user needs for a given dish diameter and data rate.

Every input is tagged [SOURCED] or [ASSUMED]. Assumed values still cite the
source that bounds or motivates them.

Link equation (all dB):
    EIRP_req + G/T_rx = Eb/N0_req + 10log10(Rb) + margin
                        + FSPL + L_pointing + L_circuit - 10log10(k)
    P_tx = EIRP_req - G_tx
"""

import numpy as np

# ---------------------------------------------------------------- constants
C_MPS = 299_792_458.0              # [SOURCED] exact SI speed of light
K_BOLTZMANN = 1.380649e-23         # [SOURCED] exact SI value (2019 redefinition), NIST CODATA
T_REF_K = 290.0                    # [SOURCED] standard noise-figure reference temperature (IEEE/IRE)

# ---------------------------------------------------------------- geometry
# Earth-Moon distance and L1-to-Earth distance are the same inputs as the optical
# budget: Condon et al. 2014 slide 3 (384,400 km) and Nayak (ed.) 2025 Table 7.1
# (356,510 km L1-to-Earth screening range).
EARTH_MOON_KM = 384_400.0          # [SOURCED] Condon et al., 2014
L1_TO_EARTH_KM = 356_510.0         # [SOURCED] Nayak (ed.), 2025, Table 7.1
L1_TO_MOON_KM = EARTH_MOON_KM - L1_TO_EARTH_KM    # derived: 27,890 km
MOON_RADIUS_KM = 1737.4            # [SOURCED] NASA NSSDC Moon Fact Sheet (volumetric mean radius)
#   https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
# Longest line of sight: a surface point on the limb as seen from L1 (zero elevation).
MAX_RANGE_M = np.sqrt(L1_TO_MOON_KM**2 - MOON_RADIUS_KM**2) * 1e3

# ---------------------------------------------------------------- link inputs
FREQUENCY_HZ = 27.25e9             # [SOURCED] middle of the 27.0-27.5 GHz lunar surface-to-orbit band
#   NASA, "LunaNet Overview": https://nescacademy.nasa.gov/downloadFile/1208
#   SFCG lunar PNT briefing: https://www.unoosa.org/documents/pdf/icg/2024/WG-B_Lunar_PNT_Jun24/LunarPNT_Jun24_02_01.pdf

REQUIRED_EB_N0_DB = 2.6           # [SOURCED] CCSDS AR4JA LDPC rate 2/3, shortest block (k=1024), BPSK/QPSK
#   Hamkins, JPL IPN Progress Report 42-184 (2011), Table 2: 1.80 dB at CWER ~1e-4 for k=1024.
#   The k=4096 and k=16384 rate-1/2 codes need less (1.25 and 0.95 dB), so 1.80 dB is the
#   highest of the three rate-1/2 block lengths. Simulated 200-iteration soft-decision decoder
#   on an AWGN channel with perfect synchronization (no implementation loss).
#   https://ipnpr.jpl.nasa.gov/progress_report/42-184/184D.pdf

POINTING_LOSS_DB = 2.0 + 2.0       # [ASSUMED by source] 2 dB off-boresight at each end
CIRCUIT_LOSS_DB = 1.0 + 2.0        # [ASSUMED by source] 1 dB transmit + 2 dB receive circuit/implementation
#   NASA LEGS performance measures (values labeled "assumed" there):
#   https://explorers.larc.nasa.gov/2023APPROBE/pdf_files/LEGS Performance Measures for SMD AOs 202503.pdf

LINK_MARGIN_DB = 3.0               # [ASSUMED] 3 dB rule-of-thumb margin policy
#   "The Role of Margin in Link Design and Optimization" (NASA S3VI knowledge base):
#   https://s3vi.ndc.nasa.gov/ssri-kb/static/resources/15-0025_A1b.pdf
#   JPL IPN Progress Report 42-222 (3 dB margin confirmed for 26 GHz near-Earth links):
#   https://ipnpr.jpl.nasa.gov/progress_report/42-222/42-222B.html

# ---------------------------------------------------------------- antennas
APERTURE_EFFICIENCY = 0.60         # [ASSUMED] within the typical 50-70% range for parabolic dishes
#   W1GHZ Antenna Handbook ch. 4 (55% reasonable, up to 70%):
#   https://mirror.unpad.ac.id/orari/orari-diklat/teknik/2.4ghz/w1ghz-antenna-handbook/chap4.pdf
#   radio-electronics.com, "Parabolic Reflector Antenna Gain" (50-70%):
#   https://www.radio-electronics.com/info/antennas/parabolic/parabolic-reflector-antenna-gain.php

L1_DISH_DIAMETER_M = 1.0           # [ASSUMED] design choice for the L1 receive dish (not from a source)

# L1 receiver system noise temperature: antenna temperature + LNA noise temperature.
ANTENNA_NOISE_K = 290.0            # [ASSUMED] conservative body-filled-beam figure, ITU-R small-sat link budget
#   https://itu.int/en/ITU-R/space/workshops/2015-prague-small-sat/Presentations/ITU-linkbudget.pdf
#   (that figure is for an Earth-viewing antenna; the lunar surface is assumed comparable, and the
#    sunlit surface is hotter, so this is optimistic by day)
RECEIVER_NOISE_FIGURE_DB = 2.6     # [ASSUMED] from a commercial Ka-band receiver (Stellant Systems, 28-31 GHz)
#   https://www.everythingrf.com/products/receivers/stellant-systems/131-1845-ka-band-communication-receiver

# ---------------------------------------------------------------- table axes
USER_DISH_DIAMETERS_M = (0.1, 0.3, 0.5, 1.0)
DATA_RATES_BPS = (1e6, 10e6, 100e6)


def db(x):
    return 10.0 * np.log10(x)


def dish_gain_db(diameter_m):
    wavelength = C_MPS / FREQUENCY_HZ
    return db(APERTURE_EFFICIENCY * (np.pi * diameter_m / wavelength) ** 2)


def free_space_path_loss_db(distance_m):
    wavelength = C_MPS / FREQUENCY_HZ
    return 20.0 * np.log10(4.0 * np.pi * distance_m / wavelength)


def system_noise_temperature_k():
    receiver_k = T_REF_K * (10.0 ** (RECEIVER_NOISE_FIGURE_DB / 10.0) - 1.0)
    return ANTENNA_NOISE_K + receiver_k


def l1_g_over_t_db():
    return dish_gain_db(L1_DISH_DIAMETER_M) - db(system_noise_temperature_k())


def required_eirp_dbw(data_rate_bps):
    """EIRP the user must radiate toward L1, independent of user dish size."""
    required_cn0_dbhz = REQUIRED_EB_N0_DB + db(data_rate_bps) + LINK_MARGIN_DB
    return (
        required_cn0_dbhz
        + free_space_path_loss_db(MAX_RANGE_M)
        + POINTING_LOSS_DB
        + CIRCUIT_LOSS_DB
        + db(K_BOLTZMANN)
        - l1_g_over_t_db()
    )


def required_tx_power_w(diameter_m, data_rate_bps):
    return 10.0 ** ((required_eirp_dbw(data_rate_bps) - dish_gain_db(diameter_m)) / 10.0)


def main():
    print(f"Range (limb of Moon seen from L1): {MAX_RANGE_M / 1e3:,.0f} km, f = {FREQUENCY_HZ / 1e9:.2f} GHz")
    print(f"L1 dish {L1_DISH_DIAMETER_M} m: G/T = {l1_g_over_t_db():.1f} dB/K "
          f"(Tsys = {system_noise_temperature_k():.0f} K)")
    print("Required EIRP (dBW): " + ", ".join(
        f"{r / 1e6:.0f} Mbps = {required_eirp_dbw(r):.1f}" for r in DATA_RATES_BPS))
    print()
    header = f"{'User dish':>10} {'Gain':>9} " + " ".join(f"{r / 1e6:>8.0f} Mbps" for r in DATA_RATES_BPS)
    print(header)
    for diameter in USER_DISH_DIAMETERS_M:
        powers = " ".join(f"{required_tx_power_w(diameter, r):>7.2f} W   " for r in DATA_RATES_BPS)
        print(f"{diameter:>8.1f} m {dish_gain_db(diameter):>6.1f} dBi {powers}")


if __name__ == "__main__":
    main()