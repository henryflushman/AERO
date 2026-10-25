#!/usr/bin/env python3
"""Calculate the tensile deliverables from the Team 13 Instron CSV exports.

INSTALL:  python -m pip install numpy matplotlib
RUN:      python tensile_lab_analysis.py --data-dir path/to/csv_files
          python tensile_lab_analysis.py --data-dir . --no-plots

Put both *_1.csv and *_2.csv exports in the data folder. The *_1 file contains
specimen-specific dimensions AND raw data; *_2 repeats the raw data. Each test
is analyzed once, using *_2 when available. No pandas or SciPy is required.

BEFORE USING RESULTS IN YOUR PRESENTATION:
1. Gauge lengths, widths, thicknesses, and Al 6061-T6 identification are filled
   from Carbon.xlsx, BONES sheet. Inches are converted using 25.4 mm/in.
   Carbon "L gage" is used as the original jaw-to-jaw length. Measurements are
   embedded below; Carbon.xlsx is not needed when running the script.
2. Inspect the *_review.png figures, set elastic_time_s=(start, end), and set
   sensor_valid_until_s to the last time the gauge/extensometer is reliable.
   Moduli without a chosen interval are PRELIMINARY. Offset yield requires
   BOTH settings; a numerical crossing alone cannot prove the gauge worked.
3. Set fracture_time_s to the last sample immediately BEFORE fracture. Without
   it, any sudden-drop estimate is explicitly a preliminary candidate. No drop
   means fracture and failure elongation are left blank.
4. If debonding is confirmed, set debond_time_s to its onset. Debond strain is
   the last strain value BEFORE that time. The automatic jump suggestion is
   only a review aid: specimen failure, slip, saturation, and electrical faults
   can resemble debonding. It is never accepted automatically.
5. Confirm material, sensor, and direction. Every supplied file says "Bonded
   Strain Gage". There is no explicitly identified extensometer channel. If a
   mislabeled test actually used one, change sensor to "extensometer" here.
   Extensometer strains are already strains; do NOT divide them by length.

UNITS / METHODS:
  A0 in mm^2; F in kN; sigma = 1000*F/A0 in MPa (engineering stress).
  epsilon_Instron = (extension - chosen zero)/original_length in mm/mm.
  Exported gauge strain is percent: epsilon = exported_value / 100.
  sigma = E*epsilon + b fitted by ordinary least squares, with an intercept.
  Offset line: sigma = E*(epsilon - 0.002) + b. Yield is the first continuous
  positive-to-negative intersection after the selected elastic interval,
  before UTS and before the sensor cutoff; interpolation is between samples.
  For a brittle composite, an offset intersection is an offset strength;
  it does not establish that the material exhibits metal-like plastic yield.
  Poisson ratio uses separate tests at a COMMON elastic stress interval:
       nu = -(d epsilon_transverse / d sigma)/(d epsilon_axial / d sigma).
  This assumes comparable material/layup and test conditions. Never pair times
  from separate tests or use a transverse gauge to compute Young's modulus.
  Elongation is an Instron estimate at fracture, not measured post-fracture
  permanent elongation; machine compliance and jaw slip can bias it.

OUTPUTS: specimen_summary.csv, poisson_summary.csv, replicate_summary.csv,
  reference_comparison.csv, analysis_notes.txt, example_calculations.txt,
  per-test cleaned CSVs, and PNG figures. Blank cells mean unavailable, not 0.
  Cleaned CSVs retain load/extension through the fracture candidate and blank
  sensor readings outside the usable sensor interval. Review plots retain raw
  data so that excluded portions remain visible for checking.

Scope: tensile calculations only. Workbook impact/fatigue data are outside
this script's scope. Theoretical composite modulus additionally requires actual
fiber/matrix properties, fiber fraction, and laminate architecture.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import numpy as np


# ======================== EDIT YOUR EXPERIMENT INPUTS ========================
# Directions inferred from Long/Par and Perp labels; confirm with lab notes.
SPECIMENS = {
    1: dict(material="aluminum", sensor="strain_gauge", direction="transverse"),
    2: dict(material="aluminum", sensor="strain_gauge", direction="axial"),
    3: dict(material="carbon", sensor="strain_gauge", direction="axial"),
    4: dict(material="carbon", sensor="strain_gauge", direction="axial"),
    5: dict(material="carbon", sensor="strain_gauge", direction="transverse"),
    6: dict(material="carbon", sensor="strain_gauge", direction="transverse"),
    7: dict(material="carbon", sensor="strain_gauge", direction="transverse"),
    8: dict(material="carbon", sensor="strain_gauge", direction="axial"),
}

# Measured defaults from Carbon.xlsx!BONES. Values are INCHES as recorded.
# Use gauge length (column C), NOT total specimen length (column D).
# Match carbon names in CSV labels, NOT workbook column G testing order:
# CSVs 3/4 are Lewis, 5/6 Austen, 7 Chesterton, 8 Tolkin/Tolkien.
# Multiple records with the same name are not automatically independent tests.
INCH_TO_MM = 25.4
GEOMETRY_IN = {
    "aluminum_axial": dict(length=2.765, width=0.4545, thickness=0.0585, row=3),
    "aluminum_transverse": dict(length=2.640, width=0.4975, thickness=0.0620, row=4),
    "Lewis": dict(length=4.117, width=0.7130, thickness=0.0670, row=7),
    "Chesterton": dict(length=4.399, width=0.6180, thickness=0.0630, row=10),
    "Austen": dict(length=4.1295, width=0.7090, thickness=0.0665, row=13),
    "Tolkien": dict(length=4.163, width=0.6805, thickness=0.0660, row=16),
    "Madilin": dict(length=4.6522, width=0.6715, thickness=0.0665, row=18),
}
# Madilin has no matching uploaded raw-data file, so is not analyzed.
SPECIMEN_GEOMETRY = {
    1: "aluminum_transverse", 2: "aluminum_axial", 3: "Lewis", 4: "Lewis",
    5: "Austen", 6: "Austen", 7: "Chesterton", 8: "Tolkien",
}
for specimen_number, geometry_name in SPECIMEN_GEOMETRY.items():
    measurement = GEOMETRY_IN[geometry_name]
    SPECIMENS[specimen_number].update(
        original_length_mm=measurement["length"] * INCH_TO_MM,
        width_mm=measurement["width"] * INCH_TO_MM,
        thickness_mm=measurement["thickness"] * INCH_TO_MM,
        area_mm2=measurement["width"] * measurement["thickness"] * INCH_TO_MM**2,
        geometry_source=f"Carbon.xlsx!BONES C{measurement['row']}:F{measurement['row']} ({geometry_name}); inches converted to mm",
    )
    if SPECIMENS[specimen_number]["material"] == "aluminum":
        SPECIMENS[specimen_number]["alloy"] = "6061-T6"

DEFAULTS = dict(
    original_length_mm=None,
    area_mm2=None,                  # Override Results Table 1 area if measured.
    width_mm=None,
    thickness_mm=None,
    geometry_source=None,
    elastic_time_s=None,            # Example (10.0, 20.0); inspect plots first.
    sensor_valid_until_s=None,      # Last reliable reading (inclusive).
    debond_time_s=None,             # Confirmed onset (exclusive sensor cutoff).
    fracture_time_s=None,           # Last intact point, before force collapse.
    start_time_s=None,              # Manual end of slack; otherwise 1% peak F.
    instron_zero_time_s=None,       # None preserves the instrument's zero.
    strain_zero_time_s=None,        # None preserves the instrument's zero.
    strain_sign=1.0,                # Set -1 only if wiring polarity is confirmed.
    strain_column=None,             # None selects Bonded Strain Gage by header.
    strain_unit=None,               # None uses CSV units; "%" or "mm/mm".
    min_fit_r2=0.98,
    force_sign=1.0,                 # Tensile force must be positive.
    include=True,
    replicate_group=None,          # Only group independent, comparable tests.
    alloy=None,                     # "1100-H14", "2024-T3", or "6061-T6".
)
# Example overrides (AFTER the geometry defaults above):
# SPECIMENS[2].update(elastic_time_s=(..., ...),
#                     sensor_valid_until_s=..., fracture_time_s=...)
# SPECIMENS[8]["sensor"] = "extensometer"  # ONLY if lab notes confirm this.
# SPECIMENS[3]["include"] = False          # ONLY if confirmed aborted/retest.

# Explicit pairs: (axial specimen number, transverse specimen number).
# Carbon is not paired automatically: confirm matching layup and usable data.
POISSON_PAIRS = [(2, 1)]
# POISSON_PAIRS += [(4, 7)]

# Values supplied in the assignment, not independently verified literature.
COURSE_REFERENCES = {
    "1100-H14": dict(E_GPa=69.0, yield_MPa=117.0, UTS_MPa=138.0, elongation_pct=20.0),
    "2024-T3": dict(E_GPa=73.0, yield_MPa=345.0, UTS_MPa=483.0, elongation_pct=18.0),
    "6061-T6": dict(E_GPa=67.0, yield_MPa=240.0, UTS_MPa=290.0, elongation_pct=8.0),
}
# Optional longitudinal rule of mixtures. Only valid for a suitable aligned
# unidirectional composite under axial loading, not a general woven laminate.
UD_COMPOSITE = None
# UD_COMPOSITE = dict(fiber_E_GPa=..., matrix_E_GPa=..., fiber_volume_fraction=...)
# ===========================================================================


def number(text):
    try:
        value = float(text)
        return value if np.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def read_export(path):
    """Read either raw-only or report-plus-raw exports; preserve CSV units."""
    with path.open(encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f))
    meta = {}
    for i, row in enumerate(rows):
        if "Specimen label" in row and "Area" in row:
            meta = dict(zip(row, rows[i + 2]))
            meta["_units"] = dict(zip(row, rows[i + 1]))
            break
    h = next((i for i, row in enumerate(rows)
              if row and row[0].strip() == "Time" and "Force" in row), None)
    if h is None:
        raise ValueError(f"No Time/Force raw-data header in {path.name}")
    names, units = rows[h], rows[h + 1]
    values = []
    for row in rows[h + 2:]:
        if len(row) != len(names) or number(row[0]) is None:
            continue
        values.append([np.nan if number(v) is None else float(v) for v in row])
    if len(values) < 5:
        raise ValueError(f"Too few raw-data rows in {path.name}")
    a = np.asarray(values)
    return {name: a[:, i] for i, name in enumerate(names)}, dict(zip(names, units)), meta


def unit_text(value):
    return value.strip().strip("()").strip().lower()


def strain_factor(unit):
    u = unit_text(unit)
    if u in ("%", "percent"):
        return 0.01
    if u in ("mm/mm", "in/in", "strain", "1", ""):
        if not u:
            raise ValueError("Missing strain units; set strain_unit explicitly")
        return 1.0
    if u in ("microstrain", "ue", "µe", "με"):
        return 1e-6
    raise ValueError(f"Unsupported strain units: {unit}")


def regress(x, y):
    """OLS y=m*x+b; report R² and n so the chosen interval is auditable."""
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    if len(x) < 5 or np.ptp(x) <= 1e-12 or np.ptp(y) <= 1e-12:
        raise ValueError("Insufficient nonconstant data for linear regression")
    m, b = np.polyfit(x, y, 1)
    residual = y - (m * x + b)
    r2 = 1.0 - np.sum(residual**2) / np.sum((y - np.mean(y))**2)
    return float(m), float(b), float(r2), len(x)


def yield_intersection(strain, stress, slope, intercept, first, last):
    """Interpolate a continuous offset intersection on the loading branch."""
    residual = stress - (slope * (strain - 0.002) + intercept)
    for j in range(max(1, first + 1), last + 1):
        if not (np.isfinite(residual[j - 1]) and np.isfinite(residual[j])):
            continue
        if strain[j] <= strain[j - 1]:
            continue
        if residual[j - 1] > 0 and residual[j] <= 0:
            f = residual[j - 1] / (residual[j - 1] - residual[j])
            return (float(strain[j - 1] + f * (strain[j] - strain[j - 1])),
                    float(stress[j - 1] + f * (stress[j] - stress[j - 1])))
    return None


def sudden_drop_candidate(t, force, peak):
    """Suggest onset of >=15% peak-force loss within 0.10 s after peak.

    This is deliberately labeled a candidate: unloading/slip can look similar.
    If the record ends without a drop, do not claim the last point is failure.
    """
    later = np.minimum(np.searchsorted(t, t + 0.10), len(t) - 1)
    eligible = (np.arange(len(t)) >= peak) & (force > 0.3 * force[peak])
    candidates = np.flatnonzero(eligible & (force - force[later] >= 0.15 * force[peak]))
    if not len(candidates):
        return None
    i = int(candidates[0])
    j = int(later[i])
    # Refine to the step preceding a fast collapse, after the candidate onset.
    dt = np.diff(t[i:j + 1])
    speed = -np.diff(force[i:j + 1]) / dt
    fast = np.flatnonzero(speed > 0.10 * np.max(speed))
    return i + int(fast[0]) if len(fast) else i


def at_time(t, value, name):
    if not np.isfinite(value) or value < t[0] or value > t[-1]:
        raise ValueError(f"{name}={value} outside recorded time interval")
    return int(np.clip(np.searchsorted(t, value, side="right") - 1, 0, len(t) - 1))


def save_csv(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields or ["status"])
        w.writeheader()
        for row in rows:
            w.writerow({k: ("" if v is None or isinstance(v, (float, np.floating))
                            and not np.isfinite(v) else v) for k, v in row.items()})


def analyze(specimen, config, folder):
    c = {**DEFAULTS, **config}
    if c["sensor"] not in ("strain_gauge", "extensometer", "none"):
        raise ValueError("sensor must be strain_gauge, extensometer, or none")
    if c["direction"] not in ("axial", "transverse"):
        raise ValueError("direction must be axial or transverse")
    stem = f"team13_tues_pm_{specimen}"
    paths = [folder / f"{stem}_{suffix}.csv" for suffix in (1, 2)]
    parsed = [read_export(p) if p.exists() else None for p in paths]
    if not any(parsed):
        raise FileNotFoundError(f"No exports found for {stem}")
    raw, units, meta2 = parsed[1] or parsed[0]
    meta = parsed[0][2] if parsed[0] else meta2
    notes = []
    if all(parsed):
        same = all(k in parsed[0][0] and parsed[0][0][k].shape == v.shape
                   and np.allclose(parsed[0][0][k], v, equal_nan=True)
                   for k, v in raw.items())
        notes.append("Paired raw exports match; analyzed once." if same else
                     "Paired raw exports differ; used _2; inspect both exports.")
    if unit_text(units["Time"]) != "s" or unit_text(units["Displacement"]) != "mm":
        raise ValueError("Expected time in seconds and displacement in mm")
    if unit_text(units["Force"]) not in ("kn", "n"):
        raise ValueError("Expected force in kN or N")
    col = c["strain_column"]
    if col is None:
        col = next((k for k in raw if "bonded strain" in k.lower()), None)
        if col is None:
            col = next((k for k in raw if "extensometer" in k.lower()), None)
    t, extension = raw["Time"].copy(), raw["Displacement"].copy()
    force = raw["Force"].copy() * c["force_sign"]
    if unit_text(units["Force"]) == "n":
        force /= 1000
    strain = np.full_like(t, np.nan)
    if c["sensor"] != "none" and col is not None:
        if col not in raw:
            raise ValueError(f"Strain column {col!r} does not exist")
        strain = raw[col].copy() * strain_factor(c["strain_unit"] or units[col]) * c["strain_sign"]
    else:
        notes.append("No usable configured strain channel.")
    good = np.isfinite(t) & np.isfinite(extension) & np.isfinite(force)
    removed = len(t) - int(good.sum())
    t, extension, force, strain = (a[good] for a in (t, extension, force, strain))
    if removed:
        notes.append(f"Removed {removed} rows with missing time/load/displacement.")
    if np.any(np.diff(t) <= 0):
        raise ValueError("Raw time must strictly increase; inspect duplicate/reversed rows")
    if c["original_length_mm"] is not None and c["original_length_mm"] <= 0:
        raise ValueError("original_length_mm must be positive")
    area = c["area_mm2"] if c["area_mm2"] is not None else number(meta.get("Area"))
    if area is None or area <= 0:
        raise ValueError("Enter positive area_mm2; no usable specimen-specific area")
    if c["area_mm2"] is None and unit_text(meta.get("_units", {}).get("Area", "")) != "mm^2":
        raise ValueError("Expected Results Table 1 area in mm^2")
    if c["geometry_source"]:
        notes.append(f"Geometry source: {c['geometry_source']}.")
        notes.append("Stress uses measured width*thickness area from workbook defaults; area_mm2 can be overridden.")
    else:
        notes.append("Stress uses Results Table 1 area, not generic inch header dimensions."
                     if c["area_mm2"] is None else "Stress uses user-entered area_mm2.")
    for key, arr in (("instron_zero_time_s", extension), ("strain_zero_time_s", strain)):
        if c[key] is not None:
            arr -= arr[at_time(t, c[key], key)]
    peak = int(np.argmax(force))
    if force[peak] <= 0:
        raise ValueError("No positive tensile load; check force_sign")
    start = (at_time(t, c["start_time_s"], "start_time_s") if c["start_time_s"] is not None
             else int(np.flatnonzero(force >= 0.01 * force[peak])[0]))
    if start >= peak:
        raise ValueError("start_time_s must precede maximum force")
    if c["fracture_time_s"] is not None:
        fracture = at_time(t, c["fracture_time_s"], "fracture_time_s")
        fracture_status = "manually selected pre-fracture point"
        if fracture < peak:
            raise ValueError("Fracture precedes recorded UTS; check test segmentation")
    else:
        fracture = sudden_drop_candidate(t, force, peak)
        fracture_status = "PRELIMINARY sudden-drop candidate" if fracture is not None else "unavailable: no clear fracture drop"
    end = fracture if fracture is not None else len(t) - 1
    stress = force * 1000 / area
    instron_strain = extension / c["original_length_mm"] if c["original_length_mm"] else np.full_like(t, np.nan)
    valid_end = end
    if c["sensor_valid_until_s"] is not None:
        valid_end = min(valid_end, at_time(t, c["sensor_valid_until_s"], "sensor_valid_until_s"))
    if c["debond_time_s"] is not None:
        at_time(t, c["debond_time_s"], "debond_time_s")  # Validate first.
        valid_end = min(valid_end, int(np.searchsorted(t, c["debond_time_s"])) - 1)
    if valid_end < start:
        raise ValueError("Sensor validity ends before usable loading begins")
    reviewed_cutoff = c["sensor_valid_until_s"] is not None or c["debond_time_s"] is not None
    fit_mask = (np.arange(len(t)) >= start) & (np.arange(len(t)) <= min(peak, valid_end))
    if c["elastic_time_s"] is not None:
        lo, hi = c["elastic_time_s"]
        if lo >= hi:
            raise ValueError("elastic_time_s must be (earlier time, later time)")
        at_time(t, lo, "elastic fit start")
        at_time(t, hi, "elastic fit end")
        if lo < t[start] or hi > t[min(peak, valid_end)]:
            raise ValueError("Elastic interval must be after slack and before peak/sensor cutoff")
        fit_mask &= (t >= lo) & (t <= hi)
    else:
        # Initial reproducible candidate only; not an automatic elastic-region claim.
        low = np.flatnonzero(force[:peak + 1] >= 0.10 * force[peak])[0]
        high = np.flatnonzero(force[:peak + 1] >= 0.30 * force[peak])[0]
        fit_mask &= (np.arange(len(t)) >= low) & (np.arange(len(t)) <= high)
        notes.append("Elastic interval PRELIMINARY: first 10-30% of peak load; review plot.")
    fit_mask &= np.isfinite(strain)
    fit, compliance = None, None
    try:
        compliance = regress(stress[fit_mask], strain[fit_mask])
        if np.unique(strain[fit_mask]).size < 5:
            raise ValueError("Strain channel is flat or has too few resolved strain levels")
        if c["direction"] == "axial":
            fit = regress(strain[fit_mask], stress[fit_mask])
            if fit[0] <= 0:
                fit = None
                raise ValueError("Nonpositive axial modulus; check wiring/channel and interval")
        if compliance[2] < c["min_fit_r2"]:
            notes.append(f"Low elastic fit R²={compliance[2]:.4f}; inspect toe/nonlinearity/noise.")
    except ValueError as e:
        fit, compliance = None, None
        notes.append(str(e))
    reviewed_fit = c["elastic_time_s"] is not None
    acceptable = (fit is not None and fit[2] >= c["min_fit_r2"])
    yield_point = None
    if fit is not None:
        if acceptable and reviewed_fit and reviewed_cutoff:
            fit_last = int(np.flatnonzero(fit_mask)[-1])
            yield_point = yield_intersection(strain, stress, fit[0], fit[1], fit_last,
                                             min(peak, valid_end))
            if yield_point is None:
                notes.append("No 0.2% offset intersection in valid pre-UTS sensor data; no yield reported.")
        else:
            notes.append("Yield withheld: select elastic_time_s and a reliable sensor cutoff; require acceptable R².")
    if not reviewed_cutoff:
        notes.append("Sensor cutoff unreviewed: cleaned/analysis sensor curve limited to fit interval.")
    if c["original_length_mm"] is None:
        notes.append("Missing original_length_mm: Instron strains and elongation are blank.")
    if fracture is None:
        notes.append("No clear fracture: failure stress/strain/elongation withheld; record may be incomplete or aborted.")
    elif c["fracture_time_s"] is None:
        notes.append("Fracture-derived quantities are PRELIMINARY until fracture_time_s is selected.")
    table_uts = number(meta.get("Tensile stress at Tensile strength"))
    if table_uts and abs(stress[peak] / table_uts - 1) > 0.01:
        notes.append("Raw-data UTS differs from exported summary UTS by >1%; raw data used. Check segmentation/area.")
    jump = None
    if np.isfinite(strain).sum() > 10 and compliance is not None:
        steps = np.abs(np.diff(strain))
        local = steps[fit_mask[:-1] & fit_mask[1:]]
        typical = np.nanmedian(local) if len(local) else 0
        # 0.0005 mm/mm is a screening threshold, not a debond criterion.
        threshold = max(0.0005, 50 * typical)
        ids = np.flatnonzero((steps > threshold) & (np.arange(len(steps)) >= start)
                             & (np.arange(len(steps)) < min(peak, end)))
        if len(ids):
            jump = int(ids[0] + 1)
            notes.append(f"Sensor jump near t={t[jump]:.6g} s; possible fault/slip/debonding, not a confirmed event.")
    debond = None
    if c["sensor"] == "strain_gauge" and c["debond_time_s"] is not None:
        debond_index = int(np.searchsorted(t, c["debond_time_s"])) - 1
        if debond_index >= 0:
            debond = strain[debond_index]
    sensor_end = valid_end if reviewed_cutoff else (int(np.flatnonzero(fit_mask)[-1]) if np.any(fit_mask) else start - 1)
    result = dict(
        specimen=specimen, label=meta.get("Specimen label", stem), material=c["material"],
        sensor=c["sensor"], direction=c["direction"], alloy=c["alloy"], area_mm2=area,
        width_mm=c["width_mm"], thickness_mm=c["thickness_mm"], geometry_source=c["geometry_source"],
        original_length_mm=c["original_length_mm"],
        ultimate_load_kN=float(force[peak]), ultimate_extension_mm=float(extension[peak]),
        UTS_MPa=float(stress[peak]), instron_strain_at_UTS=float(instron_strain[peak]),
        fracture_stress_MPa=float(stress[fracture]) if fracture is not None else None,
        instron_strain_at_fracture=float(instron_strain[fracture]) if fracture is not None else None,
        elongation_at_failure_pct=float(100 * instron_strain[fracture]) if fracture is not None else None,
        fracture_time_s=float(t[fracture]) if fracture is not None else None,
        fracture_status=fracture_status,
        E_GPa=fit[0] / 1000 if fit else None,
        elastic_intercept_MPa=fit[1] if fit else None,
        elastic_R2=compliance[2] if compliance else None,
        elastic_fit_points=compliance[3] if compliance else None,
        elastic_time_start_s=float(t[fit_mask][0]) if np.any(fit_mask) else None,
        elastic_time_end_s=float(t[fit_mask][-1]) if np.any(fit_mask) else None,
        elastic_status=(("reviewed interval" if reviewed_fit else "PRELIMINARY interval")
                        + ("; LOW R2" if compliance[2] < c["min_fit_r2"] else "")) if compliance else "unavailable",
        yield_MPa=yield_point[1] if yield_point else None,
        yield_sensor_strain=yield_point[0] if yield_point else None,
        debonding_strain=debond,
        debonding_time_s=c["debond_time_s"] if c["sensor"] == "strain_gauge" else None,
        sensor_jump_candidate_time_s=float(t[jump]) if jump is not None else None,
        sensor_cutoff_reviewed=reviewed_cutoff,
        replicate_group=c["replicate_group"], notes=" | ".join(notes),
    )
    return dict(result=result, config=c, t=t, extension=extension, force=force,
                stress=stress, strain=strain, instron_strain=instron_strain,
                start=start, peak=peak, end=end, fracture=fracture, fit_mask=fit_mask,
                fit=fit, compliance=compliance, sensor_end=sensor_end, yield_point=yield_point)


def poisson_analysis(a, b):
    """Compare strain compliances at matched stress, never matched test time."""
    ca, cb = a["config"], b["config"]
    if ca["material"] != cb["material"] or ca["direction"] != "axial" or cb["direction"] != "transverse":
        raise ValueError("Pair must be axial/transverse tests of matching material")
    if ca["sensor"] != "strain_gauge" or cb["sensor"] != "strain_gauge":
        raise ValueError("Requested Poisson calculation requires two strain-gauge tests")
    if a["compliance"] is None or b["compliance"] is None:
        raise ValueError("One strain channel has no usable elastic regression")
    # Restrict BOTH independently fitted intervals to their common stress range.
    lo = max(np.min(a["stress"][a["fit_mask"]]), np.min(b["stress"][b["fit_mask"]]))
    hi = min(np.max(a["stress"][a["fit_mask"]]), np.max(b["stress"][b["fit_mask"]]))
    if hi <= lo:
        raise ValueError("Selected elastic stress intervals do not overlap")
    fitted = []
    for data in (a, b):
        mask = data["fit_mask"] & (data["stress"] >= lo) & (data["stress"] <= hi)
        fitted.append(regress(data["stress"][mask], data["strain"][mask]))
    sa, sb = fitted[0][0], fitted[1][0]
    if sa <= 0 or sb >= 0:
        raise ValueError("Expected positive axial and negative transverse elastic compliance; inspect sign/resolution")
    if any(fit[2] < cfg["min_fit_r2"] for fit, cfg in zip(fitted, (ca, cb))):
        raise ValueError("Common-stress fit R² too low; review intervals")
    reviewed = all(d["config"]["elastic_time_s"] is not None and d["result"]["sensor_cutoff_reviewed"] for d in (a, b))
    return dict(poisson_ratio=-sb / sa, common_stress_min_MPa=lo, common_stress_max_MPa=hi,
                axial_compliance_per_MPa=sa, transverse_compliance_per_MPa=sb,
                axial_R2=fitted[0][2], transverse_R2=fitted[1][2],
                status="reviewed intervals; separate-specimen estimate" if reviewed else "PRELIMINARY separate-specimen estimate")


def make_plots(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r, t = data["result"], data["t"]
    stem = f"specimen_{r['specimen']}"
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), constrained_layout=True)
    for ax, y, label in zip(axes, (data["force"], data["strain"], data["extension"]),
                             ("Force (kN)", "Sensor strain (mm/mm)", "Instron extension (mm)")):
        ax.plot(t, y, lw=1)
        ax.set(ylabel=label, xlabel="Time (s)")
        if np.any(data["fit_mask"]):
            ts = t[data["fit_mask"]]
            ax.axvspan(ts[0], ts[-1], color="green", alpha=0.15, label="Elastic-fit interval")
        for time, color, label in ((r["fracture_time_s"], "red", "Fracture selection/candidate"),
                                   (data["config"]["debond_time_s"], "purple", "Confirmed debond onset"),
                                   (r["sensor_jump_candidate_time_s"], "orange", "Sensor jump to review"),
                                   (data["config"]["sensor_valid_until_s"], "black", "Chosen sensor cutoff")):
            if time is not None:
                ax.axvline(time, color=color, ls="--", label=label)
        ax.grid(alpha=0.25)
    axes[0].legend(fontsize=8)
    fig.suptitle(f"{r['label']} — RAW review data (not presentation curves)")
    fig.savefig(out / f"{stem}_review.png", dpi=160)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), constrained_layout=True)
    clean = slice(data["start"], data["end"] + 1)
    axes[0].plot(data["extension"][clean], data["force"][clean])
    axes[0].set(xlabel="Instron extension (mm)", ylabel="Load (kN)")
    length = r["original_length_mm"]
    x = data["instron_strain"] if length else data["extension"]
    axes[1].plot(x[clean], data["stress"][clean])
    axes[1].set(xlabel="Instron strain (mm/mm)" if length else "Instron extension (mm); original length missing",
                ylabel="Engineering stress (MPa)")
    for ax, xx, yy in ((axes[0], data["extension"], data["force"]), (axes[1], x, data["stress"])):
        for i, marker, label in ((data["peak"], "o", "UTS / maximum load"),
                                  (data["fracture"], "X", r["fracture_status"])):
            if i is not None:
                ax.scatter(xx[i], yy[i], marker=marker, s=45, label=label, zorder=4)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=7)
    fig.suptitle(r["label"])
    fig.savefig(out / f"{stem}_load_and_instron.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    sensor_end = data["sensor_end"]
    sensor_slice = slice(data["start"], sensor_end + 1)
    ax.plot(data["strain"][sensor_slice], data["stress"][sensor_slice], lw=1.3, label="Usable sensor interval")
    fit = data["fit"]
    if fit is not None:
        m, b, r2, n = fit
        xf = data["strain"][data["fit_mask"]]
        ax.plot(xf, m * xf + b, color="black", ls="--",
                label=f"E = {m/1000:.3g} GPa; R² = {r2:.4f}; n = {n}")
        hi = max(float(np.nanmax(data["strain"][sensor_slice])), float(np.max(xf)) + 0.002)
        xo = np.linspace(-b / m + 0.002, hi, 100)
        ax.plot(xo, m * (xo - 0.002) + b, ls=":", label="0.2% strain offset")
        if data["yield_point"] is not None:
            xe, ys = data["yield_point"]
            ax.scatter(xe, ys, color="red", zorder=5, label=f"Offset yield = {ys:.3g} MPa")
        ax.set_ylim(bottom=min(0, float(np.nanmin(data["stress"][sensor_slice]))),
                    top=max(float(np.nanmax(data["stress"][sensor_slice])) * 1.2, 1))
    elif data["compliance"] is not None:
        k, b, r2, n = data["compliance"]
        ss = data["stress"][data["fit_mask"]]
        ax.plot(k * ss + b, ss, "k--", label=f"Transverse compliance = {k:.3g} /MPa; R²={r2:.4f}")
    for i, label, marker in ((data["peak"], "UTS", "o"), (data["fracture"], "Fracture", "X")):
        if i is not None and i <= sensor_end and np.isfinite(data["strain"][i]):
            ax.scatter(data["strain"][i], data["stress"][i], s=45, marker=marker, label=label)
    interval_note = "\nUTS/fracture outside usable sensor interval; see Instron plot." if data["peak"] > sensor_end else ""
    ax.set(xlabel=f"{r['sensor'].replace('_', ' ')} strain (mm/mm)", ylabel="Engineering stress (MPa)")
    ax.set_title(f"{r['label']}\n{r['elastic_status']}{interval_note}", fontsize=10)
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=8, loc="best")
    ax.grid(alpha=0.25)
    fig.savefig(out / f"{stem}_sensor_stress_strain.png", dpi=170)
    plt.close(fig)


def export_cleaned(data, out):
    rows = []
    for i in range(data["start"], data["end"] + 1):
        sensor_ok = i <= data["sensor_end"]
        rows.append(dict(time_s=data["t"][i], extension_mm=data["extension"][i],
                         force_kN=data["force"][i], stress_MPa=data["stress"][i],
                         instron_strain=data["instron_strain"][i],
                         sensor_strain=data["strain"][i] if sensor_ok else None,
                         in_elastic_fit=bool(data["fit_mask"][i]), sensor_in_usable_interval=sensor_ok))
    save_csv(out / f"specimen_{data['result']['specimen']}_cleaned.csv", rows)


def plot_poisson(a, b, result, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    lo, hi = result["common_stress_min_MPa"], result["common_stress_max_MPa"]
    for data, label in ((a, "Axial"), (b, "Transverse")):
        mask = data["fit_mask"] & (data["stress"] >= lo) & (data["stress"] <= hi)
        x, y = data["stress"][mask], data["strain"][mask]
        m, intercept, _, _ = regress(x, y)
        line, = ax.plot(x, y, lw=1, label=f"{label} test {data['result']['specimen']}")
        ax.plot(x, m * x + intercept, "--", color=line.get_color())
    sa, sb = result["axial_compliance_per_MPa"], result["transverse_compliance_per_MPa"]
    ax.set(xlabel="Engineering stress (MPa)", ylabel="Sensor strain (mm/mm)",
           title=f"Poisson ratio: separate tests at common stress\nnu = -({sb:.4g}) / ({sa:.4g}) = {result['poisson_ratio']:.4f}\n{result['status']}")
    ax.legend()
    ax.grid(alpha=0.25)
    fig.savefig(out / f"poisson_{a['result']['specimen']}_{b['result']['specimen']}.png", dpi=170)
    plt.close(fig)


def aggregate_and_compare(all_data, out):
    groups, comparisons = [], []
    # Do not silently combine aborted tests, retests, different layups or sensors.
    keys = sorted({d["config"]["replicate_group"] for d in all_data.values()
                   if d["config"]["replicate_group"] is not None})
    for group in keys:
        selected = [d["result"] for d in all_data.values() if d["config"]["replicate_group"] == group]
        for prop in ("UTS_MPa", "E_GPa", "yield_MPa", "elongation_at_failure_pct"):
            vals = [r.get(prop) for r in selected if r.get(prop) is not None and np.isfinite(r[prop])]
            groups.append(dict(group=group, property=prop, n=len(vals),
                               mean=float(np.mean(vals)) if vals else None,
                               sample_std=float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
                               status="Review individual fit/fracture statuses and test independence"))
    save_csv(out / "replicate_summary.csv", groups,
             ["group", "property", "n", "mean", "sample_std", "status"])
    for data in all_data.values():
        alloy = data["config"]["alloy"]
        if alloy is None:
            continue
        if alloy not in COURSE_REFERENCES:
            raise ValueError(f"Unknown alloy {alloy!r}; update COURSE_REFERENCES")
        for prop, refprop in (("E_GPa", "E_GPa"), ("yield_MPa", "yield_MPa"),
                              ("UTS_MPa", "UTS_MPa"), ("elongation_at_failure_pct", "elongation_pct")):
            measured, reference = data["result"].get(prop), COURSE_REFERENCES[alloy][refprop]
            comparisons.append(dict(specimen=data["result"]["specimen"], alloy=alloy, property=prop,
                                    measured=measured, reference=reference,
                                    percent_difference=100 * (measured - reference) / reference if measured is not None else None,
                                    source="Assignment nominal values; check applicability and result statuses"))
    save_csv(out / "reference_comparison.csv", comparisons,
             ["specimen", "alloy", "property", "measured", "reference", "percent_difference", "source"])


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    here = Path(__file__).resolve().parent
    # Include the folder layout shown in the user's VS Code project. Defaults
    # are relative to THIS SCRIPT; an explicit --data-dir follows terminal cwd.
    candidates = [here, here / "data" / "tensile_data", here / "upload"]
    default_folder = next((p for p in candidates if any(p.glob("team13_tues_pm_*_*.csv"))), here / "data" / "tensile_data")
    parser.add_argument("--data-dir", type=Path, default=default_folder)
    parser.add_argument("--output-dir", type=Path, default=here / "tensile_results")
    parser.add_argument("--no-plots", action="store_true")
    args = parser.parse_args()
    if not args.data_dir.is_dir():
        parser.error(f"Data directory does not exist: {args.data_dir}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    all_data, errors, calculations = {}, [], []
    print("Specimen  UTS (MPa)   E (GPa)    Elastic fit status")
    for specimen, config in SPECIMENS.items():
        if not config.get("include", True):
            continue
        try:
            d = analyze(specimen, config, args.data_dir)
            export_cleaned(d, args.output_dir)
            if not args.no_plots:
                make_plots(d, args.output_dir)
            all_data[specimen] = d
            r = d["result"]
            estr = f"{r['E_GPa']:.3f}" if r["E_GPa"] is not None else "--"
            print(f"{specimen:8d}  {r['UTS_MPa']:9.3f}  {estr:>9}    {r['elastic_status']}")
            calculations.extend([
                f"\nSPECIMEN {specimen}: {r['label']}",
                f"Geometry source: {r['geometry_source']}",
                f"A0 = width*thickness = {r['width_mm']}*{r['thickness_mm']} = {r['area_mm2']:.6g} mm^2; L0 = {r['original_length_mm']} mm",
                f"A0 = {r['area_mm2']:.6g} mm^2; Fmax = {r['ultimate_load_kN']:.6g} kN",
                f"UTS = 1000*{r['ultimate_load_kN']:.6g}/{r['area_mm2']:.6g} = {r['UTS_MPa']:.6g} MPa",
                f"Extension at Fmax = {r['ultimate_extension_mm']:.6g} mm",
                f"epsilon_Instron at UTS = extension/L0 = {r['instron_strain_at_UTS']}",
                f"Elastic fit: {r['elastic_time_start_s']} to {r['elastic_time_end_s']} s; {r['elastic_status']}",
                f"E = slope/1000 = {r['E_GPa']} GPa; intercept = {r['elastic_intercept_MPa']} MPa; R^2 = {r['elastic_R2']}",
                f"Offset line sigma = E_MPa*(epsilon - 0.002) + intercept; intersection = ({r['yield_sensor_strain']}, {r['yield_MPa']} MPa)",
                f"Fracture: {r['fracture_status']}; time = {r['fracture_time_s']} s; stress = {r['fracture_stress_MPa']} MPa",
                f"Instron elongation at fracture = 100*extension/L0 = {r['elongation_at_failure_pct']} %",
                f"Confirmed last pre-debond sensor strain = {r['debonding_strain']}",
            ])
        except (ValueError, FileNotFoundError, KeyError) as e:
            errors.append(f"Specimen {specimen}: {e}")
            print(f"{specimen:8d}  SKIPPED: {e}")
    if not all_data:
        raise SystemExit("No specimens analyzed. " + "; ".join(errors))
    poisson_rows = []
    for axial, transverse in POISSON_PAIRS:
        row = dict(axial_specimen=axial, transverse_specimen=transverse)
        try:
            if axial not in all_data or transverse not in all_data:
                raise ValueError("A paired specimen was excluded or could not be analyzed")
            row.update(poisson_analysis(all_data[axial], all_data[transverse]))
            if not args.no_plots:
                plot_poisson(all_data[axial], all_data[transverse], row, args.output_dir)
            calculations.append(f"\nPoisson pair ({axial}, {transverse}): nu = -{row['transverse_compliance_per_MPa']:.6g}/{row['axial_compliance_per_MPa']:.6g} = {row['poisson_ratio']:.6g}; {row['status']}")
        except ValueError as e:
            row.update(poisson_ratio=None, status=f"unavailable: {e}")
        poisson_rows.append(row)
    save_csv(args.output_dir / "specimen_summary.csv", [d["result"] for d in all_data.values()])
    save_csv(args.output_dir / "poisson_summary.csv", poisson_rows)
    aggregate_and_compare(all_data, args.output_dir)
    notes = ["TENSILE ANALYSIS REVIEW NOTES", __doc__,
             "No extensometer method is inferred from a Bonded Strain Gage header.",
             "Course references contain no Poisson ratio; add a traceable reference for your actual material/layup.",
             "Carbon pairings and replicate groups require lab-note confirmation; specimen 3 may be an aborted/retest record.",
             "Geometry defaults: Carbon.xlsx!BONES, matched by name (not testing order). Carbon L gage is interpreted as jaw-to-jaw length.",
             "Madilin geometry (BONES row 18) is recorded in GEOMETRY_IN, but no matching raw test is available.",
             "Aluminum is identified as 6061-T6 in BONES!B2; comparison uses the assignment's nominal values."]
    for d in all_data.values():
        notes.extend([f"\nSPECIMEN {d['result']['specimen']}: {d['result']['label']}", d["result"]["notes"]])
    notes.extend(["\nERRORS / SKIPPED SPECIMENS", *errors])
    if UD_COMPOSITE is not None:
        p = UD_COMPOSITE
        vf = p["fiber_volume_fraction"]
        if not 0 <= vf <= 1 or p["fiber_E_GPa"] <= 0 or p["matrix_E_GPa"] <= 0:
            raise ValueError("Invalid UD composite constituent properties")
        theoretical = vf * p["fiber_E_GPa"] + (1 - vf) * p["matrix_E_GPa"]
        calculations.append(f"\nALIGNED UD COMPOSITE ONLY: E1 = Vf*Ef + (1-Vf)*Em = {theoretical:.6g} GPa")
    (args.output_dir / "analysis_notes.txt").write_text("\n".join(notes), encoding="utf-8")
    (args.output_dir / "example_calculations.txt").write_text("\n".join(calculations), encoding="utf-8")
    print(f"\nOutputs: {args.output_dir.resolve()}")
    print("Workbook geometry defaults loaded. Review elastic/sensor/fracture settings before reporting results.")
    if errors:
        raise SystemExit("Some specimens were skipped; see analysis_notes.txt")


if __name__ == "__main__":
    main()