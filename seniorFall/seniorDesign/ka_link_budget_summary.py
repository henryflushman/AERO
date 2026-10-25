"""Earth downlink: SuperDish assumptions, mesh reflectors, 5 Gbps useful data.

Run: python superdish_mesh_5gbps.py
Standard library only. Earth link only.

Retains the numerical Earth-link assumptions in the uploaded SuperDish.m.
Changes spacecraft dish sizes to 12/18/22 m, useful rate to 5 Gbps,
and frequency to the center of the selected 25.5-27 GHz band.
Separates average RF output from amplifier electrical input and removes the
invalid capacity claim. Mesh performance is represented by assumed net aperture
efficiency; no specific reflector's flight performance is certified here.
"""

import math

# CONSTRAINT VALUES FOR THIS COMPARISON
# Useful customer bits/s. This is now 5 Gbps, rather than the original 2 Gbps.
# For numerical fidelity, the original Eb/N0 is treated as referenced to THESE
# useful bits. No additional coding/framing overhead is invented. A real modem
# threshold must use this same energy/bit definition or be converted to it.
USEFUL_DATA_RATE_BPS = 18e9
# Earth receiving dish diameter [m], unchanged from SuperDish.m.
EARTH_DISH_DIAMETER_M = 18.0
# Original assumed Earth-center-to-L1 distance [km] and Earth radius [km].
# The original cites Nayak (ed.), 2025, Table 7.1 for the first quantity; that
# citation has not been verified here. The tangent-range construction is kept.
# This is not a worst-case range propagated over an actual NRHO trajectory.
L1_TO_EARTH_KM = 384400.0
EARTH_RADIUS_KM = 6378.0
# SUBJECT TO CHANGE -- ORIGINAL ASSUMPTIONS RETAINED
# Spacecraft deployed mesh reflector diameters [m]. These replace the original
# 5 m dish. Mesh construction enters through the EFFECTIVE aperture efficiency;
# do not subtract an additional mesh/surface loss if already included in it.
# These are size candidates, not verified exact-frequency flight products.
SC_MESH_DISH_DIAMETERS_M = [12, 18, 22]
# Selected frequency band [Hz]. Use its midpoint for this single-carrier budget.
# These are the study band edges, not a confirmed licensed channel assignment.
BAND_LOW_HZ = 25.5e9
BAND_HIGH_HZ = 27.0e9
FREQUENCY_HZ = (BAND_LOW_HZ + BAND_HIGH_HZ) / 2  # 26.25 GHz
# Original channel bandwidth [Hz]. Without modulation, coding, rolloff and
# modem limits, bandwidth alone does NOT establish 5 Gbps waveform feasibility.
BANDWIDTH_HZ = BAND_HIGH_HZ - BAND_LOW_HZ  # 1.5 GHz upper envelope
# Original receiver threshold [dB], retained without adding implementation loss.
# Its modulation, code and error-rate basis are unspecified in SuperDish.m.
# This is Eb/N0, NOT Es/N0 or SNR. Here it is interpreted per useful payload bit.
REQUIRED_EB_N0_DB = 12.05
# Original combined spacecraft + ground pointing losses [dB]. Fixed losses
# across the dish sizes imply progressively tighter angular pointing accuracy.
# No beamwidth-based pointing calculation is introduced in this reproduction.
POINTING_LOSS_DB = 2.0 + 2.0
# Original combined circuit-loss allowance [dB]. Its TX/RX allocation is not
# specified, so it remains lumped in the required-power equation exactly as
# in the original. It must NOT also be added to the noise model without first
# establishing consistent receiver reference planes and thermal contributions.
CIRCUIT_LOSS_DB = 1.0 + 2.0
# Original atmospheric attenuation [dB]. No site/weather availability has been
# established for this value. It remains 10 dB rather than your Python's 2 dB.
ATM_LOSS_DB = 10.0
# Original additional signal reserve [dB], above the receiver threshold.
LINK_MARGIN_DB = 3.0
# Original effective aperture efficiency for BOTH antennas, unchanged at 0.60.
# For the mesh reflectors this is a provisional net efficiency, NOT a supplier
# guarantee. Feed taper, mesh reflection, surface errors and blockage affect it.
APERTURE_EFFICIENCY = 0.60
# Original antenna noise temperature [K] and receiver noise figure [dB].
# Their sum through the noise-figure conversion gives Tsys ~527.7 K.
# This preserves the original simplified model; receiving-feed thermal noise,
# sky/weather consistency and lunar beam coupling have not been resolved.
ANTENNA_NOISE_K = 290.0
RECEIVER_NOISE_FIGURE_DB = 2.6
# Original 40% RF-output / amplifier-electrical-input efficiency, previously
# hidden inside the power function. This is not aperture efficiency or backoff.
# Electrical input here excludes modem, converters, gimbals and spacecraft loads.
AMPLIFIER_EFFICIENCY = 0.40

# CONSTANTS
C_MPS = 299792458.0
K_BOLTZMANN = 1.380649e-23
T_REF_K = 290.0  # Noise-figure reference temperature in the original model.


def db(value):
    return 10 * math.log10(value)


def dish_gain_dbi(diameter_m):
    wavelength_m = C_MPS / FREQUENCY_HZ
    return db(APERTURE_EFFICIENCY * (math.pi * diameter_m / wavelength_m)**2)


def validate_inputs():
    for name in ("USEFUL_DATA_RATE_BPS", "EARTH_DISH_DIAMETER_M", "L1_TO_EARTH_KM",
                 "EARTH_RADIUS_KM", "FREQUENCY_HZ", "BANDWIDTH_HZ", "ANTENNA_NOISE_K"):
        value = globals()[name]
        if not math.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if not (math.isfinite(BAND_LOW_HZ) and math.isfinite(BAND_HIGH_HZ)
            and 0 < BAND_LOW_HZ < BAND_HIGH_HZ):
        raise ValueError("Band edges must be positive, finite and ordered")
    if (FREQUENCY_HZ - BANDWIDTH_HZ/2 < BAND_LOW_HZ
            or FREQUENCY_HZ + BANDWIDTH_HZ/2 > BAND_HIGH_HZ):
        raise ValueError("The channel extends outside the selected frequency band")
    if not SC_MESH_DISH_DIAMETERS_M:
        raise ValueError("Supply at least one mesh dish size")
    for diameter in SC_MESH_DISH_DIAMETERS_M:
        if not math.isfinite(diameter) or diameter <= 0:
            raise ValueError("Mesh dish diameters must be finite and positive")
    if L1_TO_EARTH_KM <= EARTH_RADIUS_KM:
        raise ValueError("Tangent range requires spacecraft outside Earth")
    for name in ("APERTURE_EFFICIENCY", "AMPLIFIER_EFFICIENCY"):
        value = globals()[name]
        if not math.isfinite(value) or not 0 < value <= 1:
            raise ValueError(f"{name} must be in (0, 1]")
    for name in ("POINTING_LOSS_DB", "CIRCUIT_LOSS_DB", "ATM_LOSS_DB",
                 "LINK_MARGIN_DB", "RECEIVER_NOISE_FIGURE_DB"):
        value = globals()[name]
        if not math.isfinite(value) or value < 0:
            raise ValueError(f"{name} must be finite and nonnegative")
    if not math.isfinite(REQUIRED_EB_N0_DB):
        raise ValueError("Required Eb/N0 must be finite")


def calculate_budget():
    """Return link inputs derived in the original convention and dish results."""
    validate_inputs()
    range_m = math.sqrt(L1_TO_EARTH_KM**2 - EARTH_RADIUS_KM**2) * 1000
    wavelength_m = C_MPS / FREQUENCY_HZ
    receiver_noise_k = T_REF_K * (10**(RECEIVER_NOISE_FIGURE_DB/10) - 1)
    system_noise_k = ANTENNA_NOISE_K + receiver_noise_k
    ground_gain = dish_gain_dbi(EARTH_DISH_DIAMETER_M)
    ground_gt = ground_gain - db(system_noise_k)
    fspl = 20 * math.log10(4 * math.pi * range_m / wavelength_m)
    required_cn0 = REQUIRED_EB_N0_DB + db(USEFUL_DATA_RATE_BPS) + LINK_MARGIN_DB

    # Preserve the original 'required EIRP' equation, including ALL circuit loss.
    # This is an equivalent amplifier-output-referenced quantity = P_RF[dBW]+Gtx.
    # Actual antenna EIRP needs the TX share of circuit/feed loss subtracted.
    # The original does not specify the split of TX and RX circuit losses.
    equivalent_eirp = (required_cn0 + fspl + POINTING_LOSS_DB + CIRCUIT_LOSS_DB
                       + ATM_LOSS_DB + db(K_BOLTZMANN) - ground_gt)
    rows = []
    for diameter in SC_MESH_DISH_DIAMETERS_M:
        gain = dish_gain_dbi(diameter)
        rf_w = 10**((equivalent_eirp - gain)/10)
        electrical_w = rf_w / AMPLIFIER_EFFICIENCY
        rows.append({"diameter_m": diameter, "gain_dbi": gain,
                     "rf_output_w": rf_w, "amplifier_electrical_w": electrical_w})
    return {"range_m": range_m, "wavelength_m": wavelength_m,
            "receiver_noise_k": receiver_noise_k, "system_noise_k": system_noise_k,
            "ground_gain_dbi": ground_gain, "ground_gt_dbk": ground_gt,
            "fspl_db": fspl, "required_cn0_dbhz": required_cn0,
            "equivalent_eirp_dbw": equivalent_eirp, "rows": rows}


def main():
    b = calculate_budget()
    print("EARTH -- MESH DISH TRADE, SUPERDISH ASSUMPTIONS WITH UPDATED BAND")
    print(f"Range (Earth tangent): {b['range_m']/1000:,.0f} km, f = {FREQUENCY_HZ/1e9:.2f} GHz")
    print(f"Earth dish {EARTH_DISH_DIAMETER_M:.1f} m: G/T = {b['ground_gt_dbk']:.2f} dB/K (Tsys = {b['system_noise_k']:.2f} K)")
    print(f"Useful data rate: {USEFUL_DATA_RATE_BPS/1e9:.2f} Gbps; channel bandwidth: {BANDWIDTH_HZ/1e9:.2f} GHz")
    print(f"Required Eb/N0: {REQUIRED_EB_N0_DB:.2f} dB per useful bit; reserve: {LINK_MARGIN_DB:.2f} dB")
    print(f"Required equivalent EIRP (original convention): {b['equivalent_eirp_dbw']:.2f} dBW")
    print(f"Amplifier RF-to-electrical efficiency: {AMPLIFIER_EFFICIENCY:.0%}")
    print(f"\n{'Mesh dish':>12}{'Gain [dBi]':>15}{'RF output [W]':>18}{'Amp electrical [W]':>22}")
    for row in b['rows']:
        print(f"{row['diameter_m']:>9.1f} m{row['gain_dbi']:>15.2f}{row['rf_output_w']:>18.2f}{row['amplifier_electrical_w']:>22.2f}")
    print("\nLOSS / NOISE AUDIT")
    for label, value in (("Free-space path loss", b['fspl_db']),
                         ("Pointing (combined)", POINTING_LOSS_DB),
                         ("Circuit (combined)", CIRCUIT_LOSS_DB),
                         ("Atmosphere", ATM_LOSS_DB), ("Design reserve", LINK_MARGIN_DB)):
        print(f"{label}: {value:.3f} dB")
    print(f"Antenna / receiver noise: {ANTENNA_NOISE_K:.2f} / {b['receiver_noise_k']:.2f} K")
    print("\nBandwidth feasibility: UNVERIFIED (waveform/overhead/modem unspecified).")
    print("No capacity is inferred by treating Eb/N0 as SNR.")
    print(f"Selected band: {BAND_LOW_HZ/1e9:g}-{BAND_HIGH_HZ/1e9:g} GHz; carrier at band center.")
    print("Mesh efficiency, fixed pointing losses, amplifier efficiency and threshold remain assumptions.")
    print("No extra polarization, component-feed, implementation or reuse-interference losses were added.")
    print("Electrical input excludes modem, converters, gimbals and other spacecraft loads.")


if __name__ == '__main__':
    main()
