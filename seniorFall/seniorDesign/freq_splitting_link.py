#!/usr/bin/env python3
"""Earth link budget matching the simplified link-budget explanation.

Run:
    python earth_link_budget_5gbps.py
    python earth_link_budget_5gbps.py --json
    python earth_link_budget_5gbps.py --output-json results.json

This is a preliminary sizing model, not a simulation of the custom modem.
The 30 dB downlink requirement is an estimate based on an ETSI receiver
benchmark. The extra 3 dB is design margin and is added only once.
Station comparisons reuse the supplied losses as provisional scenarios; they
are not verified atmospheric predictions at the current channel frequencies.
"""

import argparse
import copy
import json
import math
from pathlib import Path

C_M_S = 299_792_458.0
K_BOLTZMANN_J_K = 1.380649e-23
NOISE_FIGURE_REFERENCE_K = 290.0

SOURCES = {
    "lunanet": {
        "title": "NASA LunaNet Interoperability Specification, Version 5",
        "url": "https://www.nasa.gov/wp-content/uploads/2025/02/lunanet-interoperability-specification-v5-baseline.pdf",
        "location": "Appendix B, Tables B-35 to B-38",
        "use": "Candidate uplink/downlink bands; does not authorize this mission.",
    },
    "snr_benchmark": {
        "title": "ETSI TS 102 991 V1.3.1 (2016), DVB-C2 implementation guidelines",
        "url": "https://www.etsi.org/deliver/etsi_ts/102900_102999/102991/01.03.01_60/ts_102991v010301p.pdf",
        "location": "Table 26, printed page 142; Table 6, printed page 36",
        "use": "Table 26 gives 30.1 dB for long-block 1024-QAM 5/6 in the Gaussian-channel receiver case. We round to 30 dB as a preliminary estimate, not a measured custom-modem threshold. Table 6 gives 27.15 dB at BER 1e-6. DVB-C2 is a cable/OFDM benchmark, not this lunar single-carrier modem.",
    },
    "mesh_efficiency": {
        "title": "Chahat et al., CubeSat Deployable Ka-Band Mesh Reflector Antenna Development for Earth Science Missions",
        "url": "https://ieeexplore.ieee.org/document/7440821",
        "use": "Reports 52% aperture efficiency at 35.75 GHz; only a benchmark for the assumed 50%, not proof for an 18-22 m antenna.",
    },
    "noise_temperature": {
        "title": "Analog Devices, System Noise-Figure Analysis for Modern Radio Receivers",
        "url": "https://www.analog.com/en/resources/technical-articles/system-noisefigure-analysis-for-modern-radio-receivers.html",
        "use": "Noise factor, equivalent noise temperature and cascaded receiver analysis.",
    },
    "ground_lna_example": {
        "title": "Analog Devices HMC752",
        "url": "https://www.analog.com/en/products/hmc752.html",
        "use": "2.5 dB LNA benchmark at 24-28 GHz; does not establish the assumed 2.6 dB complete receiver.",
    },
    "satellite_lna_example": {
        "title": "Qorvo QPA2628",
        "url": "https://www.qorvo.com/products/p/QPA2628",
        "use": "1.6 dB typical LNA at 22-31.5 GHz; does not establish the assumed 3 dB complete satellite receiver or flight life.",
    },
    "qam_example": {
        "title": "Analog Devices, High Dynamic Range IF Receiver Simplifies Design of Next-Generation Microwave Point-to-Point Modems",
        "url": "https://www.analog.com/en/resources/technical-articles/high-dynamic-range-if-receiver-simplifies-design-of-next-generation-uw-point-to-point-modems.html",
        "use": "112 MHz 1024-QAM terrestrial implementation; supports the technique, not our full bandwidth or spacecraft design.",
    },
    "station_identity": {
        "title": "Official station identity references (not antenna access or loss data)",
        "urls": {
            "Goldstone_Madrid_Canberra": "https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301O.pdf",
            "Malargue": "https://www.esa.int/Newsroom/Press_Releases/From_la_Pampa_to_the_Universe_ESA_completes_deep-space_network",
            "Hartebees_Hartebeesthoek": "https://www.sarao.ac.za/about/hartrao/hartrao-research-programmes/hartrao-26m-radio-telescope-details/",
            "Usuda": "https://ddor.nav.isas.jaxa.jp/station_coordinates/",
        },
        "use": "Identifies the proposed sites only. These sources do not establish our assumed 18 m ground dish, station access, current Ka-band capability, or numerical atmospheric losses.",
    },
    "station_loss_numbers": {
        "title": "User-supplied ground-station comparison table",
        "url": None,
        "use": "Numerical source for the six station losses. Uplink losses were supplied at 22.85 GHz and downlink losses at 37.5 GHz, with reported P.618 5% exceedance and fallback behavior. Elevation and complete model inputs were not supplied, so these numbers cannot be independently reproduced from the pasted table.",
    },
    "itu_p618": {
        "title": "ITU-R P.618, Earth-space propagation prediction methods",
        "url": "https://www.itu.int/rec/R-REC-P.618/en",
        "use": "Supports the propagation-model method, not the exact pasted numerical station losses. No P.618 model is executed by this script.",
    },
    "itu_p676": {
        "title": "ITU-R P.676, Attenuation by atmospheric gases and related effects",
        "url": "https://www.itu.int/rec/R-REC-P.676/en",
        "use": "Supports frequency-dependent gaseous attenuation modeling. The supplied 37.5 GHz losses do not establish attenuation in our 25.5-27 GHz downlink band.",
    },
    "amplifier_effects": {
        "title": "CCSDS 130.11-G-2, SCCC Summary of Definition and Performance",
        "url": "https://ccsds.org/Pubs/130x11g2.pdf",
        "use": "Amplifier nonlinearity, filtering and backoff; no specific backoff is claimed for this custom waveform.",
    },
}

# CONSTRAINTS: mission targets and candidate bands. Change the targets here
# if the architecture changes. Bands still need coordination/authorization.
CONSTRAINTS = {
    "simultaneous_satellites": 2,
    "useful_downlink_rate_bps": 5e9,
    "useful_uplink_rate_bps": 1e9,
    "uplink_band_hz": [22.55e9, 23.15e9],  # LunaNet source above.
    "downlink_band_hz": [25.5e9, 27.0e9],
    "spacecraft_design_life_years": 15,      # User requirement, not qualification.
}

# DESIGN ASSUMPTIONS: editable preliminary values, not hardware guarantees.
DESIGN_ASSUMPTIONS = {
    "guard_fraction_of_slot": 0.05,       # 5% total, split between both edges.
    "ground_diameter_m": 18.0,
    "spacecraft_diameter_options_m": [18.0, 22.0],
    "ground_aperture_efficiency": 0.60,
    "spacecraft_aperture_efficiency": 0.50,
    "slant_range_m": 384_347e3,           # Retained input, not verified NRHO maximum.
    "pointing_loss_db": 4.0,             # Total, 2 dB allocated at each end.
    "downlink_tx_circuit_loss_db": 1.0,  # Assumed split of the original 3 dB.
    "downlink_rx_feed_loss_db": 2.0,
    "uplink_tx_circuit_loss_db": 1.0,
    "uplink_rx_feed_loss_db": 2.0,
    "earth_rx_feed_physical_temperature_k": 290.0,
    "satellite_rx_feed_physical_temperature_k": 290.0,
    "atmospheric_loss_db": 10.0,         # Adverse scenario, no availability assigned.
    "earth_antenna_noise_temperature_k": 290.0,
    "earth_receiver_noise_figure_db": 2.6,
    "satellite_antenna_noise_temperature_k": 290.0,
    "satellite_receiver_noise_figure_db": 3.0,
    "downlink_required_snr_db": 30.0,     # Rounded ETSI benchmark; see SOURCES.
    "uplink_required_snr_db": None,      # No selected/verified uplink waveform yet.
    "design_margin_db": 3.0,
    "satellite_amplifier_efficiency": 0.40,  # Must hold at actual operating backoff.
    "ground_amplifier_efficiency": 0.40,
    "modulation_order": 1024,
    "symbol_rate_baud": 640e6,
    "ldpc_code_rate": 5.0 / 6.0,
    "rolloff": 0.10,
    "overhead_fraction_after_coding": 0.05,  # Framing, pilots, customer packet overhead.
}

# STATION INPUTS: only the atmospheric-loss scenario and assumed ground dish
# differ from the main calculation. The 18 m dishes are design assumptions,
# not a claim about the installed or available antennas at these facilities.
# Numerical source: the user's pasted comparison, NOT an independent ITU data
# table. Applying its losses unchanged to our actual frequencies is provisional.
STATIONS = [
    {"name": "Goldstone", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 2.59, "downlink_atmospheric_loss_db": 2.07},
    {"name": "Malargue", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 2.34, "downlink_atmospheric_loss_db": 2.38},
    {"name": "Madrid", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 3.18, "downlink_atmospheric_loss_db": 3.49},
    {"name": "Hartebees", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 4.67, "downlink_atmospheric_loss_db": 5.21},
    {"name": "Usuda", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 6.34, "downlink_atmospheric_loss_db": 8.93},
    {"name": "Canberra", "ground_diameter_m": 18.0,
     "uplink_atmospheric_loss_db": 4.31, "downlink_atmospheric_loss_db": 5.39},
]
STATION_LOSS_PROVENANCE = {
    "source": "User-supplied comparison table",
    "source_uplink_frequency_hz": 22.85e9,
    "source_downlink_frequency_hz": 37.5e9,
    "reported_loss_percent_exceeded": 5.0,
    "status": "Provisional loss scenarios; not reproduced or validated at current frequencies",
    "assumption": "Supplied losses are held constant across the current uplink/downlink channels for comparison only. Rerun a site-specific propagation model before treating the results as expected station requirements.",
}

NOTES = [
    "30 dB is an estimated downlink receiver requirement; 33 dB is the target including the separate 3 dB design margin.",
    "SNR is C/(k*Tsys*B), with B equal to the 712.5 MHz allowed signal envelope, not the 704 MHz ideal waveform width.",
    "The custom waveform capacity is about 5.067 Gbps; required power and Eb/N0 use the 5.000 Gbps customer target.",
    "The original 3 dB circuit allowance is assigned as 1 dB transmit loss plus 2 dB receive feed loss in each direction. This split is a modeling assumption, not a measured RF chain. Noise is referenced to the LNA input after that receive feed; its thermal emission is included.",
    "The assumed 290 K antenna temperatures are inputs, not automatically the physical temperature of a cold sky or a spacecraft antenna.",
    "Atmospheric noise is represented only through the assumed antenna noise temperature; it is not added again. 10 dB atmospheric attenuation has no assigned outage probability.",
    "RF power is the required average modulated transmitter output before the lumped circuit-loss allowance. Required saturated amplifier rating is higher because of backoff; no backoff value is established here.",
    "Electrical power is RF power divided by assumed operating efficiency. It excludes modem, converter, thermal-control and other spacecraft loads.",
    "A 15-year life requires mission-specific radiation, reliability, thermal and end-of-life assessment; it cannot be demonstrated by this link budget.",
    "No practical uplink power is reported until a separate uplink receiver SNR requirement is supplied. The downlink 30 dB estimate is not copied to uplink.",
    "Station results use our present channel frequencies, dishes, 33 dB target and RF equations, substituting only the station atmospheric-loss scenario and its editable assumed ground dish. Powers from the old model are not copied.",
    "The six numerical losses are user-supplied scenarios originally at 22.85 GHz uplink and 37.5 GHz downlink. Their reported 5% exceedance cannot be verified without the omitted elevation and model inputs; the 37.5 GHz losses are not verified values for our 26 GHz downlink.",
    "Station noise temperatures and range remain fixed at the original budget inputs. This isolates the loss/dish comparison; atmospheric noise, geometry and elevation should be modeled per station for final predictions.",
]


def db10(value):
    """Convert a positive power ratio to dB."""
    return 10.0 * math.log10(value)


def equivalent_receiver_temperature(noise_figure_db):
    return NOISE_FIGURE_REFERENCE_K * (10.0 ** (noise_figure_db / 10.0) - 1.0)


def antenna_gain_dbi(diameter_m, frequency_hz, efficiency):
    wavelength_m = C_M_S / frequency_hz
    return db10(efficiency * (math.pi * diameter_m / wavelength_m) ** 2)


def free_space_loss_db(range_m, frequency_hz):
    return 20.0 * math.log10(4.0 * math.pi * range_m * frequency_hz / C_M_S)


def shannon_snr_db(rate_bps, bandwidth_hz):
    return db10(math.expm1(math.log(2.0) * rate_bps / bandwidth_hz))


def validate(constraints, design):
    count = constraints["simultaneous_satellites"]
    if not isinstance(count, int) or isinstance(count, bool) or count < 1:
        raise ValueError("simultaneous_satellites must be a positive integer")
    for key in ["useful_downlink_rate_bps", "useful_uplink_rate_bps", "spacecraft_design_life_years"]:
        if not math.isfinite(constraints[key]) or constraints[key] <= 0:
            raise ValueError(f"{key} must be positive and finite")
    for key in ["uplink_band_hz", "downlink_band_hz"]:
        lo, hi = constraints[key]
        if not (math.isfinite(lo) and math.isfinite(hi) and 0 < lo < hi):
            raise ValueError(f"{key} must contain increasing positive frequencies")
    for key in ["ground_aperture_efficiency", "spacecraft_aperture_efficiency", "satellite_amplifier_efficiency", "ground_amplifier_efficiency", "ldpc_code_rate"]:
        if not 0 < design[key] <= 1:
            raise ValueError(f"{key} must be in (0, 1]")
    for key in ["guard_fraction_of_slot", "overhead_fraction_after_coding"]:
        if not 0 <= design[key] < 1:
            raise ValueError(f"{key} must be in [0, 1)")
    for key in ["ground_diameter_m", "slant_range_m", "symbol_rate_baud"]:
        if not math.isfinite(design[key]) or design[key] <= 0:
            raise ValueError(f"{key} must be positive and finite")
    if not design["spacecraft_diameter_options_m"] or any(not math.isfinite(v) or v <= 0 for v in design["spacecraft_diameter_options_m"]):
        raise ValueError("spacecraft_diameter_options_m must contain positive values")
    for key in ["pointing_loss_db", "downlink_tx_circuit_loss_db", "downlink_rx_feed_loss_db", "uplink_tx_circuit_loss_db", "uplink_rx_feed_loss_db", "earth_rx_feed_physical_temperature_k", "satellite_rx_feed_physical_temperature_k", "atmospheric_loss_db", "design_margin_db", "rolloff", "earth_antenna_noise_temperature_k", "satellite_antenna_noise_temperature_k", "earth_receiver_noise_figure_db", "satellite_receiver_noise_figure_db"]:
        if not math.isfinite(design[key]) or design[key] < 0:
            raise ValueError(f"{key} must be finite and nonnegative")
    order = design["modulation_order"]
    if not isinstance(order, int) or order < 2 or order & (order - 1):
        raise ValueError("modulation_order must be an integer power of two")
    for key in ["downlink_required_snr_db", "uplink_required_snr_db"]:
        if design[key] is not None and not math.isfinite(design[key]):
            raise ValueError(f"{key} must be finite or None")
    if design["downlink_required_snr_db"] is None:
        raise ValueError("A downlink SNR requirement is needed")


def make_channels(constraints, design):
    channels = []
    count = constraints["simultaneous_satellites"]
    for direction, key in [("uplink", "uplink_band_hz"), ("downlink", "downlink_band_hz")]:
        band_lo, band_hi = constraints[key]
        slot = (band_hi - band_lo) / count
        guard_each_edge = slot * design["guard_fraction_of_slot"] / 2.0
        for index in range(count):
            lo = band_lo + index * slot
            hi = lo + slot
            channels.append({
                "direction": direction,
                "satellite_id": index + 1,
                "slot_low_hz": lo,
                "slot_high_hz": hi,
                "center_frequency_hz": (lo + hi) / 2.0,
                "slot_bandwidth_hz": slot,
                "envelope_low_hz": lo + guard_each_edge,
                "envelope_high_hz": hi - guard_each_edge,
                "available_bandwidth_hz": slot - 2.0 * guard_each_edge,
                "guard_each_edge_hz": guard_each_edge,
            })
    return channels


def station_comparison_rows(base_rows, stations, design, direction):
    """Apply station scenarios using the same receiver-power equation as baseline.

    The source frequency is recorded separately from the actual channel frequency.
    This is a what-if loss comparison, not an atmospheric-model calculation.
    """
    output = []
    down = direction == "downlink"
    efficiency = design["satellite_amplifier_efficiency"] if down else design["ground_amplifier_efficiency"]
    loss_key = "downlink_atmospheric_loss_db" if down else "uplink_atmospheric_loss_db"
    source_key = "source_downlink_frequency_hz" if down else "source_uplink_frequency_hz"
    for station in stations:
        for key in ["ground_diameter_m", "uplink_atmospheric_loss_db", "downlink_atmospheric_loss_db"]:
            value = station[key]
            if not math.isfinite(value) or value < 0 or (key == "ground_diameter_m" and value == 0):
                raise ValueError(f"Invalid {key} for station {station['name']}")
        for base in base_rows:
            row = copy.deepcopy(base)
            row["station_name"] = station["name"]
            row["ground_diameter_m"] = station["ground_diameter_m"]
            row["atmospheric_loss_db"] = station[loss_key]
            row["atmospheric_loss_source"] = STATION_LOSS_PROVENANCE["source"]
            row["atmospheric_loss_source_frequency_hz"] = STATION_LOSS_PROVENANCE[source_key]
            row["loss_scenario_status"] = STATION_LOSS_PROVENANCE["status"]
            row["ground_gain_dbi"] = antenna_gain_dbi(
                station["ground_diameter_m"], row["center_frequency_hz"],
                design["ground_aperture_efficiency"],
            )
            row["other_signal_losses_db"] = (
                design["pointing_loss_db"] + row["tx_circuit_loss_db"]
                + row["rx_feed_loss_db"] + station[loss_key]
            )
            if row["target_snr_db"] is not None:
                # Solve Pt from C = Pt + Gt + Gr - FSPL - other losses.
                row["required_carrier_dbw"] = row["noise_power_dbw"] + row["target_snr_db"]
                row["required_rf_power_dbw"] = (
                    row["required_carrier_dbw"] + row["free_space_loss_db"]
                    + row["other_signal_losses_db"] - row["spacecraft_gain_dbi"]
                    - row["ground_gain_dbi"]
                )
                row["required_rf_power_w"] = 10.0 ** (row["required_rf_power_dbw"] / 10.0)
                row["amplifier_electrical_power_w"] = row["required_rf_power_w"] / efficiency
            output.append(row)
    return output


def build_results(constraints=None, design=None, stations=None):
    """Return all document values as a JSON-compatible dictionary.

    Optional dictionaries replace the full defaults; edit the two config
    dictionaries above or pass modified copies to study different assumptions.
    """
    c = copy.deepcopy(CONSTRAINTS if constraints is None else constraints)
    a = copy.deepcopy(DESIGN_ASSUMPTIONS if design is None else design)
    station_inputs = copy.deepcopy(STATIONS if stations is None else stations)
    validate(c, a)
    channels = make_channels(c, a)
    te_earth = equivalent_receiver_temperature(a["earth_receiver_noise_figure_db"])
    te_sat = equivalent_receiver_temperature(a["satellite_receiver_noise_figure_db"])
    earth_feed_loss = 10.0 ** (a["downlink_rx_feed_loss_db"] / 10.0)
    satellite_feed_loss = 10.0 ** (a["uplink_rx_feed_loss_db"] / 10.0)
    # Available antenna noise is attenuated through the feed, while the feed
    # adds thermal noise. All temperatures below are at the LNA input.
    earth_antenna_at_lna = (a["earth_antenna_noise_temperature_k"] / earth_feed_loss
        + a["earth_rx_feed_physical_temperature_k"] * (1.0 - 1.0 / earth_feed_loss))
    satellite_antenna_at_lna = (a["satellite_antenna_noise_temperature_k"] / satellite_feed_loss
        + a["satellite_rx_feed_physical_temperature_k"] * (1.0 - 1.0 / satellite_feed_loss))
    t_earth = earth_antenna_at_lna + te_earth
    t_sat = satellite_antenna_at_lna + te_sat
    if min(t_earth, t_sat) <= 0:
        raise ValueError("System noise temperatures must be positive")
    bandwidth = next(x["available_bandwidth_hz"] for x in channels if x["direction"] == "downlink")
    bits_per_symbol = int(math.log2(a["modulation_order"]))
    coded_rate = a["symbol_rate_baud"] * bits_per_symbol
    information_rate = coded_rate * a["ldpc_code_rate"]
    useful_capacity = information_rate * (1.0 - a["overhead_fraction_after_coding"])
    occupied_bandwidth = a["symbol_rate_baud"] * (1.0 + a["rolloff"])
    waveform = {
        "modulation": f'{a["modulation_order"]}-QAM',
        "bits_per_symbol": bits_per_symbol,
        "symbol_rate_baud": a["symbol_rate_baud"],
        "ldpc_code_rate": a["ldpc_code_rate"],
        "rolloff": a["rolloff"],
        "overhead_fraction_after_coding": a["overhead_fraction_after_coding"],
        "coded_rate_bps": coded_rate,
        "information_rate_after_coding_bps": information_rate,
        "useful_capacity_bps": useful_capacity,
        "occupied_bandwidth_hz": occupied_bandwidth,
        "available_bandwidth_hz": bandwidth,
        "bandwidth_headroom_hz": bandwidth - occupied_bandwidth,
        "payload_headroom_bps": useful_capacity - c["useful_downlink_rate_bps"],
        "bandwidth_fits": occupied_bandwidth <= bandwidth,
        "payload_rate_fits": useful_capacity >= c["useful_downlink_rate_bps"],
    }
    downlinks = []
    uplinks = []
    for channel in channels:
        down = channel["direction"] == "downlink"
        tx_loss = a["downlink_tx_circuit_loss_db"] if down else a["uplink_tx_circuit_loss_db"]
        rx_loss = a["downlink_rx_feed_loss_db"] if down else a["uplink_rx_feed_loss_db"]
        loss = a["pointing_loss_db"] + tx_loss + rx_loss + a["atmospheric_loss_db"]
        t_sys = t_earth if down else t_sat
        rate = c["useful_downlink_rate_bps"] if down else c["useful_uplink_rate_bps"]
        b = channel["available_bandwidth_hz"]
        f = channel["center_frequency_hz"]
        noise = db10(K_BOLTZMANN_J_K * t_sys * b)
        fspl = free_space_loss_db(a["slant_range_m"], f)
        ground_gain = antenna_gain_dbi(a["ground_diameter_m"], f, a["ground_aperture_efficiency"])
        required_snr = a["downlink_required_snr_db"] if down else a["uplink_required_snr_db"]
        shannon = shannon_snr_db(rate, b)
        efficiency = a["satellite_amplifier_efficiency"] if down else a["ground_amplifier_efficiency"]
        for diameter in a["spacecraft_diameter_options_m"]:
            sat_gain = antenna_gain_dbi(diameter, f, a["spacecraft_aperture_efficiency"])
            target = None if required_snr is None else required_snr + a["design_margin_db"]
            rf_dbw = None if target is None else noise + target + fspl + loss - ground_gain - sat_gain
            rf_w = None if rf_dbw is None else 10.0 ** (rf_dbw / 10.0)
            row = {
                "satellite_id": channel["satellite_id"],
                "sat_diameter_m": diameter,
                "center_frequency_hz": f,
                "bandwidth_hz": b,
                "useful_rate_bps": rate,
                "spacecraft_gain_dbi": sat_gain,
                "ground_gain_dbi": ground_gain,
                "free_space_loss_db": fspl,
                "other_signal_losses_db": loss,
                "tx_circuit_loss_db": tx_loss,
                "rx_feed_loss_db": rx_loss,
                "noise_reference_plane": "LNA input, after receive feed",
                "system_temperature_k": t_sys,
                "noise_power_dbw": noise,
                "required_snr_db": required_snr,
                "design_margin_db": a["design_margin_db"],
                "target_snr_db": target,
                "required_ebn0_useful_db": None if required_snr is None else required_snr + db10(b / rate),
                "target_ebn0_useful_db": None if target is None else target + db10(b / rate),
                "shannon_snr_db": shannon,
                "shannon_ebn0_useful_db": shannon + db10(b / rate),
                "required_carrier_dbw": None if target is None else noise + target,
                "required_rf_power_dbw": rf_dbw,
                "required_rf_power_w": rf_w,
                "amplifier_electrical_power_w": None if rf_w is None else rf_w / efficiency,
            }
            if down:
                downlinks.append(row)
            else:
                uplinks.append(row)
    return {
        "config": {
            "constraints": c,
            "design_assumptions": a,
            "earth_receiver_equivalent_temperature_k": te_earth,
            "satellite_receiver_equivalent_temperature_k": te_sat,
            "earth_antenna_temperature_at_lna_k": earth_antenna_at_lna,
            "satellite_antenna_temperature_at_lna_k": satellite_antenna_at_lna,
            "earth_system_temperature_k": t_earth,
            "satellite_system_temperature_k": t_sat,
        },
        "waveform": waveform,
        "channels": channels,
        "station_inputs": station_inputs,
        "station_loss_provenance": copy.deepcopy(STATION_LOSS_PROVENANCE),
        "station_downlinks": station_comparison_rows(downlinks, station_inputs, a, "downlink"),
        "station_uplinks": station_comparison_rows(uplinks, station_inputs, a, "uplink"),
        "downlinks": downlinks,
        "uplink": {"status": "Threshold not selected" if a["uplink_required_snr_db"] is None else "User-supplied threshold; waveform still needs validation", "rows": uplinks},
        "notes": NOTES,
        "sources": SOURCES,
    }


def print_table(headers, rows):
    rows = [[str(value) for value in row] for row in rows]
    widths = [max(len(headers[i]), max((len(row[i]) for row in rows), default=0)) for i in range(len(headers))]
    print(" | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))))
    print("-+-".join("-" * width for width in widths))
    for row in rows:
        print(" | ".join(row[i].ljust(widths[i]) for i in range(len(headers))))


def print_results(results):
    """Print only the final power comparison for each spacecraft dish size."""
    reference = results["downlinks"][0]
    print(
        f'FINAL DOWNLINK POWER: '
        f'{reference["required_snr_db"]:g} dB SNR + '
        f'{reference["design_margin_db"]:g} dB margin = '
        f'{reference["target_snr_db"]:g} dB target\n'
    )
    print_table(
        ["Spacecraft dish m", "Satellite", "Frequency GHz", "RF output W", "Amplifier electrical W"],
        [
            [
                f'{row["sat_diameter_m"]:g}',
                row["satellite_id"],
                f'{row["center_frequency_hz"]/1e9:.3f}',
                f'{row["required_rf_power_w"]:.1f}',
                f'{row["amplifier_electrical_power_w"]:.1f}',
            ]
            for row in sorted(
                results["downlinks"],
                key=lambda row: (row["sat_diameter_m"], row["satellite_id"]),
            )
        ],
    )
    if results["uplink"]["rows"][0]["required_snr_db"] is not None:
        print("\nFINAL UPLINK POWER: user-supplied receiver threshold")
        print_table(
            ["Spacecraft dish m", "Satellite", "RF output W", "Amplifier electrical W"],
            [
                [f'{r["sat_diameter_m"]:g}', r["satellite_id"],
                 f'{r["required_rf_power_w"]:.1f}', f'{r["amplifier_electrical_power_w"]:.1f}']
                for r in sorted(results["uplink"]["rows"],
                                key=lambda row: (row["sat_diameter_m"], row["satellite_id"]))
            ],
        )
    print("\nSTATION DOWNLINK POWER: provisional supplied-loss scenarios")
    print_station_power_table(results["station_downlinks"], results["config"]["design_assumptions"])
    if results["uplink"]["rows"][0]["required_snr_db"] is not None:
        print("\nSTATION UPLINK POWER: user-supplied receiver threshold")
        print_station_power_table(results["station_uplinks"], results["config"]["design_assumptions"])
    print("\nStation losses came from a 22.85/37.5 GHz reference model; using them unchanged at our channel frequencies is provisional.")
    waveform = results["waveform"]
    if not waveform["bandwidth_fits"] or not waveform["payload_rate_fits"]:
        print("\nWARNING: The waveform does not meet the selected bandwidth/rate targets.")


def print_station_power_table(rows, design):
    """Group final RF and electrical powers by station and frequency channel."""
    diameters = design["spacecraft_diameter_options_m"]
    headers = ["Station", "Sat", "GHz", "Atm dB", "Ground m"]
    headers += [f"RF {diameter:g}m W" for diameter in diameters]
    headers += [f"DC {diameter:g}m W" for diameter in diameters]
    groups = {}
    for row in rows:
        key = (row["station_name"], row["satellite_id"])
        groups.setdefault(key, {})[row["sat_diameter_m"]] = row
    printable = []
    for (name, satellite_id), group in groups.items():
        first = next(iter(group.values()))
        values = [name, satellite_id, f'{first["center_frequency_hz"]/1e9:.3f}',
                  f'{first["atmospheric_loss_db"]:.2f}', f'{first["ground_diameter_m"]:g}']
        for field in ["required_rf_power_w", "amplifier_electrical_power_w"]:
            values += ["TBD" if group[d][field] is None else f'{group[d][field]:.1f}' for d in diameters]
        printable.append(values)
    print_table(headers, printable)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--json", action="store_true", help="Print complete results as JSON instead of text tables")
    parser.add_argument("--output-json", type=Path, help="Optionally save the complete results as JSON")
    parser.add_argument("--uplink-snr-db", type=float, help="Explicit assumed uplink receiver SNR threshold; adds the configured 3 dB margin")
    args = parser.parse_args()
    design = copy.deepcopy(DESIGN_ASSUMPTIONS)
    if args.uplink_snr_db is not None:
        design["uplink_required_snr_db"] = args.uplink_snr_db
    try:
        results = build_results(design=design)
    except (ValueError, OverflowError) as error:
        parser.error(str(error))
    if args.json:
        print(json.dumps(results, indent=2, allow_nan=False))
    else:
        print_results(results)
    if args.output_json:
        args.output_json.write_text(json.dumps(results, indent=2, allow_nan=False) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
