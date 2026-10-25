"""Lunar link using the numerical assumptions in SuperDish.m.

Run: python superdish_lunar_link.py
Requires Python's standard library only.

DIRECTION: L1 satellite TRANSMITS; lunar user RECEIVES. This reproduces the
MATLAB equations, whose 'surface-to-L1' comment conflicts with their antenna
roles. It does not calculate the reverse user-to-satellite link.

Separate RF/electrical labels correct the original ambiguous power table.
The original invalid .332*bandwidth*EbN0 capacity claim is not reproduced.
"""

import math

# ==================== CONSTRAINT VALUES FOR THIS COMPARISON ====================
# Original assumed maximum satellite-to-lunar-user slant range [m]. The MATLAB
# calls this the Moon-limb range. Retained as provided, not verified against NRHO.
MAX_RANGE_M = 86_500 * 1e3

# Original lunar useful-data-rate case [bits/s], NOT the Earth link's 5 Gbps.
# For numerical fidelity, Eb/N0 is interpreted per useful bit at these rates.
# Add more rates to this list for a trade. No explicit waveform overhead is added.
DATA_RATES_BPS = [9e9]

# Receiving lunar user's physical dish diameter [m], original value.
USER_DISH_DIAMETER_M = 1.0

# ==================== SUBJECT TO CHANGE: ORIGINAL ASSUMPTIONS ==================
# Transmitting satellite dish diameters [m]. Original lunar case is 3 m.
# Add sizes to trade. No manufacturer or mesh technology is implied for this dish.
SAT_DISH_DIAMETERS_M = [5.0]

# Original lunar carrier frequency [Hz]. Keep 40 GHz; the selected 25.5-27 GHz
# band applies to your EARTH downlink, not this independent lunar comparison.
FREQUENCY_HZ = 38.75e9

# Original channel bandwidth [Hz]. Without modulation, code, rolloff, framing
# and modem limits, this does not establish a usable-data-rate ceiling.
BANDWIDTH_HZ = 3.5e9

# Original required energy per bit / noise spectral density [dB]. Its waveform
# and error-rate basis are unspecified; assume useful-bit energy for this copy.
# This is not symbol Es/N0 or signal/noise across the entire 2 GHz bandwidth.
REQUIRED_EB_N0_DB = 5.0

# Original combined satellite + lunar-user pointing allowance [dB]. Not derived
# from beamwidth, angular error or a statistical pointing-accuracy specification.
POINTING_LOSS_DB = 2.0 + 2.0

# Original combined circuit allowance [dB]. Its transmitter/receiver split is
# unspecified. Retained in the same lumped power equation as MATLAB; assigning
# receiver feed thermal noise later requires consistent noise reference planes.
CIRCUIT_LOSS_DB = 1.0 + 2.0 + 1.0

# Original extra signal reserve [dB], added above the specified Eb/N0 threshold.
LINK_MARGIN_DB = 3.0

# Original effective aperture efficiency [fraction], applied to BOTH antennas.
# Assumed net efficiency, not a verified supplier specification. Do not add again
# any surface, mesh, blockage or illumination losses already represented by it.
APERTURE_EFFICIENCY = 0.60

# Original receiving antenna noise [K] and receiver noise figure [dB]. The latter
# is converted using the 290 K reference, then added to antenna noise. This is
# the original simplified total, not a resolved lunar environment/feed model.
ANTENNA_NOISE_K = 290.0
RECEIVER_NOISE_FIGURE_DB = 2.6

# Original RF-output / amplifier-electrical-input efficiency [fraction]. MATLAB
# divided power by 0.4 inside its transmit-power function. Here print both powers.
# Electrical input excludes modem, converters, mechanisms and spacecraft loads.
AMPLIFIER_EFFICIENCY = 0.40

# No atmospheric, extra polarization, implementation or cochannel-interference
# loss is added: the original lunar equations contain none of these terms.

# ==================== CONSTANTS ====================
C_MPS = 299792458.0
K_BOLTZMANN = 1.380649e-23
T_REF_K = 290.0


def db(value):
    return 10 * math.log10(value)


def dish_gain_dbi(diameter_m):
    wavelength = C_MPS / FREQUENCY_HZ
    return db(APERTURE_EFFICIENCY * (math.pi * diameter_m / wavelength)**2)


def validate_inputs():
    for name in ("MAX_RANGE_M", "USER_DISH_DIAMETER_M", "FREQUENCY_HZ",
                 "BANDWIDTH_HZ", "ANTENNA_NOISE_K"):
        value = globals()[name]
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not DATA_RATES_BPS or not SAT_DISH_DIAMETERS_M:
        raise ValueError("Supply at least one data rate and satellite dish size")
    for value in list(DATA_RATES_BPS) + list(SAT_DISH_DIAMETERS_M):
        if not math.isfinite(value) or value <= 0:
            raise ValueError("Rates and dish sizes must be finite and positive")
    for name in ("APERTURE_EFFICIENCY", "AMPLIFIER_EFFICIENCY"):
        value = globals()[name]
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"{name} must be in (0, 1]")
    for name in ("POINTING_LOSS_DB", "CIRCUIT_LOSS_DB", "LINK_MARGIN_DB",
                 "RECEIVER_NOISE_FIGURE_DB"):
        value = globals()[name]
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if not math.isfinite(REQUIRED_EB_N0_DB):
        raise ValueError("Required Eb/N0 must be finite")


def calculate_budget():
    validate_inputs()
    wavelength = C_MPS / FREQUENCY_HZ
    receiver_noise = T_REF_K * (10**(RECEIVER_NOISE_FIGURE_DB/10) - 1)
    system_noise = ANTENNA_NOISE_K + receiver_noise
    user_gain = dish_gain_dbi(USER_DISH_DIAMETER_M)
    user_gt = user_gain - db(system_noise)
    fspl = 20 * math.log10(4 * math.pi * MAX_RANGE_M / wavelength)
    rows = []
    for diameter in SAT_DISH_DIAMETERS_M:
        gain = dish_gain_dbi(diameter)
        for rate in DATA_RATES_BPS:
            required_cn0 = REQUIRED_EB_N0_DB + db(rate) + LINK_MARGIN_DB
            # Same EIRP expression as MATLAB: includes ALL lumped circuit loss.
            # Equivalent EIRP = amplifier RF[dBW] + dish gain[dBi]. Conventional
            # antenna EIRP requires subtracting the separately identified TX
            # feed-loss portion, which is not specified in the original.
            eirp = (required_cn0 + fspl + POINTING_LOSS_DB + CIRCUIT_LOSS_DB
                    + db(K_BOLTZMANN) - user_gt)
            rf = 10**((eirp - gain)/10)
            rows.append({"diameter_m": diameter, "gain_dbi": gain,
                         "useful_rate_bps": rate, "required_cn0_dbhz": required_cn0,
                         "equivalent_eirp_dbw": eirp, "rf_output_w": rf,
                         "amplifier_electrical_w": rf / AMPLIFIER_EFFICIENCY})
    return {"wavelength_m": wavelength, "receiver_noise_k": receiver_noise,
            "system_noise_k": system_noise, "user_gain_dbi": user_gain,
            "user_gt_dbk": user_gt, "fspl_db": fspl, "rows": rows}


def main():
    b = calculate_budget()
    print("MOON -- L1 SATELLITE TRANSMITS TO LUNAR USER")
    print(f"Range: {MAX_RANGE_M/1000:,.0f} km, f = {FREQUENCY_HZ/1e9:.2f} GHz")
    print(f"User receiving dish {USER_DISH_DIAMETER_M:.2f} m: gain = {b['user_gain_dbi']:.2f} dBi")
    print(f"G/T = {b['user_gt_dbk']:.2f} dB/K (Tsys = {b['system_noise_k']:.2f} K)")
    print(f"Channel: {BANDWIDTH_HZ/1e9:g} GHz; required Eb/N0: {REQUIRED_EB_N0_DB:g} dB; reserve: {LINK_MARGIN_DB:g} dB")
    print(f"Amplifier RF-to-electrical efficiency: {AMPLIFIER_EFFICIENCY:.0%}")
    print(f"\n{'Sat dish [m]':>13}{'Gain [dBi]':>13}{'Useful [Gbps]':>15}{'Equiv EIRP [dBW]':>19}{'RF [W]':>13}{'Amp elec [W]':>16}")
    for row in b['rows']:
        print(f"{row['diameter_m']:>13.2f}{row['gain_dbi']:>13.2f}{row['useful_rate_bps']/1e9:>15.3f}"
              f"{row['equivalent_eirp_dbw']:>19.2f}{row['rf_output_w']:>13.2f}{row['amplifier_electrical_w']:>16.2f}")
    print("\nLOSS / NOISE AUDIT")
    for label, value in (("Free-space path loss", b['fspl_db']),
                         ("Pointing (combined)", POINTING_LOSS_DB),
                         ("Circuit (combined)", CIRCUIT_LOSS_DB),
                         ("Design reserve", LINK_MARGIN_DB)):
        print(f"{label}: {value:.3f} dB")
    print(f"Antenna / receiver noise: {ANTENNA_NOISE_K:.2f} / {b['receiver_noise_k']:.2f} K")
    print("\nBandwidth feasibility UNVERIFIED: waveform, coding, overhead and modem unspecified.")
    print("No capacity is inferred by substituting Eb/N0 for SNR.")
    print("Equivalent EIRP retains the original lumped circuit-loss convention.")
    print("Electrical input excludes modem, converters, mechanisms and other spacecraft loads.")
    print("This models satellite-to-user transmission; reverse-link receiver assumptions differ.")


if __name__ == '__main__':
    main()
