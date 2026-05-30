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
# ║   Course      :  AERO446 - Spacecraft Electrical and        ║
# ║                                     Electric Systems        ║
# ║   Assignment  :  Homework 5                                 ║
# ║   Date        :  May 15, 2026                               ║
# ╚═════════════════════════════════════════════════════════════╝


# === Imports ==================================================
import numpy as np
from dataclasses import dataclass

# ── Helper functions ──────────────────────
def section(title):
    width = 52
    print(f"\n{'═' * width}")
    print(f"  {title}")
    print(f"{'═' * width}")

def row(label, value, unit=""):
    print(f"  {label:<28} {value:>10.4f}  {unit}")
    
def linear_to_db(linear_value: float) -> float:
    """
    Converts a linear value to decibels (dB)
    """
    if type(linear_value) != float and type(linear_value) != int:
        raise TypeError("Linear value must be a number (int or float) to convert to dB.")
    if linear_value <= 0:
        raise ValueError("Linear value must be greater than zero to convert to dB.")
    return 10*np.log10(linear_value)

def db_to_linear(db_value: float) -> float:
    """
    Converts a decibel (dB) value to linear
    """
    if type(db_value) != float and type(db_value) != int:
        raise TypeError("dB value must be a number (int or float) to convert to linear.")
    return 10**(db_value/10)

@dataclass
class Amplifier:
    """
    Stores amplifier gain and noise figure in dB
    """
    name: str
    gain_db: float
    noise_figure_db: float
    
    @property
    def gain_linear(self) -> float:
        return db_to_linear(self.gain_db)
    @property
    def noise_factor_linear(self) -> float:
        return db_to_linear(self.noise_figure_db)

# ─────────────────────────────────────────

# ─────────────────────────────────────────
# Problem 1 — Gain from Frequency
# ─────────────────────────────────────────

# === Problem givens ===
f = 300  # MHz
eff_f = 0.60


def parabolic_antenna_gain_db(diameter_m: float, frequency_hz: float, efficiency: float) -> dict:
    """
    Parabolic antenna gain calculation:
    
    where,
        lambda = c / f
        A = pi * (D/2)^2
        A_eff = efficiency * *args
        G = 4*pi*A_eff / lambda^2
        G_db = 10*log10(G)
    """
    c = 3.00e8
    
    wavelength_m = c / frequency_hz
    physical_area_m2 = np.pi * (diameter_m / 2)**2
    effective_area_m2 = efficiency * physical_area_m2
    
    gain_linear = (4 * np.pi * effective_area_m2) / (wavelength_m**2)
    gain_db = linear_to_db(gain_linear)
    
    return {
        "wavelength_m": wavelength_m,
        "physical_area_m2": physical_area_m2,
        "effective_area_m2": effective_area_m2,
        "gain_linear": gain_linear,
        "gain_db": gain_db
    }
    
# ─────────────────────────────────────────
# Problem 2 — Short Answer Questions
# ─────────────────────────────────────────

def noise_figure_example(snr_in_linear: float, snr_out_linear: float) -> dict:
    """
    Noise factor:
        F = SNR_in / SNR_out
    Noise figure:
        NF_db = 10*log10(F)
    """
    noise_factor = snr_in_linear / snr_out_linear
    noise_figure_db = linear_to_db(noise_factor)
    
    return {
        "snr_in_linear": snr_in_linear,
        "snr_out_linear": snr_out_linear,
        "noise_factor": noise_factor,
        "noise_figure_db": noise_figure_db
    }
    
def db_power_percent_increase(db_value: float) -> float:
    """
    Convert dB power ratio to percent increase over original power.
    """
    linear_ratio = db_to_linear(db_value)
    return (linear_ratio - 1) * 100

# ─────────────────────────────────────────
# Problem 3 — Cascaded Receiver Analysis
# ─────────────────────────────────────────
        
def cascaded_reciever(order: str, amplifiers: dict[str, Amplifier]) -> dict:
    """
    Total gain:
        G_total = G1 * G2 * ... * Gn
        
    In dB, gains add:
        G_total_db = G1_db + G2_db + ... + Gn_db
    
    Total noise factor using Friis cascade equation:
        F_total = F1 + (F2 - 1)/G1 + (F3 - 1)/(G1*G2) + ... + (Fn - 1)/(G1*G2*...*Gn-1)
        
    Note:
        Gains and noise figures must be converted to linear units before using the Friis equation.
    """
    gain_total_linear = 1.0
    gain_total_db = 0.0
    
    for stage_name in order:
        amp = amplifiers[stage_name]
        gain_total_linear += amp.gain_linear
        gain_total_db += amp.gain_db
    
    first_amp = amplifiers[order[0]]
    noise_factor_total_linear = first_amp.noise_factor_linear
    
    gain_product_before_stage = 1.0
    
    for i in range(1, len(order)):
        previous_stage_name = order[i-1]
        current_stage_name = order[i]
        
        previous_amp = amplifiers[previous_stage_name]
        current_amp = amplifiers[current_stage_name]
        
        gain_product_before_stage *= previous_amp.gain_linear
        noise_factor_total_linear += (
            current_amp.noise_factor_linear - 1.0
        ) / gain_product_before_stage
        
    return {
        "order": order,
        "gain_total_linear": gain_total_linear,
        "gain_total_db": gain_total_db,
        "noise_factor_total_linear": noise_factor_total_linear,
        "noise_figure_total_db": linear_to_db(noise_factor_total_linear)
    }
    
def print_receiver_table(title: str, orders: list[str], amplifiers: dict[str, Amplifier]) -> list[dict]:
    """
    Table output for results
    """
    print("\n" + title)
    print("-" * len(title))
    print(f"{'Order':<8}{'G_total (dB)':>15}{'F_total':>15}{'NF_total (dB)':>18}")
    
    results = []
    
    for order in orders:
        result = cascaded_reciever(order, amplifiers)
        results.append(result)
        
        print(
            f"{result['order']:<8}"
            f"{result['gain_total_db']:>15.2f}"
            f"{result['noise_factor_total_linear']:>15.6f}"
            f"{result['noise_figure_total_db']:>18.4f}"
        )
        
    return results

# ─────────────────────────────────────────
# Problem 4 — Comparing A to A+
# ─────────────────────────────────────────

def compare_a_to_a_plus(orders: list[str], original_amplifiers: dict[str, Amplifier]) -> None:
    """
    Compare original A reciever to A+ reciever with improved noise figure in the first stage.
    """
    a_plus_amplifiers = original_amplifiers.copy()
    a_plus_amplifiers["A"] = Amplifier(name="A+", gain_db=10.0, noise_figure_db=0.8)
    
    original_results = {
        result["order"]: result
        for result in print_receiver_table("Original A results", orders, original_amplifiers)
    }
    
    plus_results = {
        result["order"]: result
        for result in print_receiver_table("A+ replacement results", orders, a_plus_amplifiers)
    }
    
    print("\nProblem 4 comparison")
    print("--------------------")
    print(f"{'Order':<8}{'Gain change (dB)':>18}{'NF change (dB)':>18}")
    
    for order in orders:
        gain_change = plus_results[order]["gain_total_db"] - original_results[order]["gain_total_db"]
        nf_change = plus_results[order]["noise_figure_total_db"] - original_results[order]["noise_figure_total_db"]
        
        print(f"{order:<8}{gain_change:>18.2f}{nf_change:>18.4f}")
        
    print("\nRecommendation:")
    print("A+ improves the noise figure most when it is earlier in the chain,")
    print("but it always reduces total gain by 20 dB.")
    print("I would not recommend A+ unless the design has a strict noise-figure requirement")
    print("and can tolerate or replace the lost 20 dB of gain")

# ─────────────────────────────────────────
# Main Script
# ─────────────────────────────────────────

def main() -> None:
    # Problem 1
    section("Problem 1 — Gain from Frequency")
    
    p1 = parabolic_antenna_gain_db(
        diameter_m=1.0,
        frequency_hz=300e6,
        efficiency=0.60
    )
    
    row("Wavelength", p1["wavelength_m"], "m")
    row("Physical Area", p1["physical_area_m2"], "m^2")
    row("Effective Area", p1["effective_area_m2"], "m^2")
    row("Gain (linear)", p1["gain_linear"])
    row("Gain (dB)", p1["gain_db"], "dB")
    
    # Problem 2
    section("Problem 2 — Short Answer Questions")
    
    normal_nf = noise_figure_example(snr_in_linear=100, snr_out_linear=80)
    negative_nf_example = noise_figure_example(snr_in_linear=100, snr_out_linear=125)
    
    row("Normal case NF (dB)", normal_nf["noise_figure_db"], "dB")
    row("Negative NF case NF (dB)", negative_nf_example["noise_figure_db"], "dB")
    
    percent_0p1_db = db_power_percent_increase(0.1)
    percent_1_db = db_power_percent_increase(1.0)
    ratio_of_percent_increases = percent_0p1_db / percent_1_db * 100
    
    row("Percent increase for 0.1 dB", percent_0p1_db, "%")
    row("Percent increase for 1 dB", percent_1_db, "%")
    row("Ratio of percent increases", ratio_of_percent_increases, "%")
    
    # Problem 3
    section("Problem 3 — Cascaded Receiver Analysis")
    
    amplifiers = {
        "A": Amplifier(name="A", gain_db=30.0, noise_figure_db=3.0),
        "B": Amplifier(name="B", gain_db=20, noise_figure_db=2),
        "C": Amplifier(name="C", gain_db=13, noise_figure_db=1.5),
    }

    orders = ["ABC", "ACB", "BAC", "BCA", "CAB"]

    print_receiver_table("Problem 3 receiver configurations", orders, amplifiers)

    # -----------------------------
    # Problem 4
    # -----------------------------
    compare_a_to_a_plus(orders, amplifiers)


if __name__ == "__main__":
    main()