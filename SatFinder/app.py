import os
import sys
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import streamlit as st

# Allow importing the sibling PolySpace repo without installing it.
POLYSPACE_SRC = Path(r"C:\VScode\PolySpace\src")
if POLYSPACE_SRC.exists() and str(POLYSPACE_SRC) not in sys.path:
    sys.path.insert(0, str(POLYSPACE_SRC))

try:
    from polyspace.SATCAT.satcat import SatCat
    from polyspace.SpaceObject.SpaceObject import SpaceObject
except Exception as exc:  # pragma: no cover - surfaced in app UI
    SatCat = None
    SpaceObject = None
    IMPORT_ERROR = exc
else:
    IMPORT_ERROR = None


def normalize_key(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip().lower().replace(" ", "_")


def normalize_record(record: dict) -> dict:
    out: dict[str, Any] = {}
    for key, val in (record or {}).items():
        nk = normalize_key(key)
        out[nk] = val
    return out


def coalesce_keys(record: dict, *keys: str) -> Any:
    for key in keys:
        if key in record and record[key] not in (None, "", "nan"):
            return record[key]
    return None


def normalize_intldes(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip().upper()


def build_query_window(epoch_range_days: int):
    end = datetime.utcnow()
    start = end - timedelta(days=max(epoch_range_days, 1))
    return start, end, f"{start:%Y-%m-%d}--{end:%Y-%m-%d}"


def chunk_list(items, chunk_size: int = 100):
    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def extract_identifier_values(df: pd.DataFrame, candidate_columns):
    values = []
    for col in candidate_columns:
        if col not in df.columns:
            continue
        values.extend(
            str(v).strip().upper()
            for v in df[col].dropna().tolist()
            if str(v).strip()
        )
    unique = []
    seen = set()
    for value in values:
        if value and value not in seen:
            seen.add(value)
            unique.append(value)
    return unique


def get_first_matching_column(df: pd.DataFrame, candidates):
    for name in candidates:
        if name in df.columns:
            return name
    return None


def fetch_space_track_rows(username: str, password: str, epoch_range_days: int):
    if not (username and password):
        raise ValueError("Space-Track username and password are required.")

    if SatCat is None:
        raise ImportError(f"PolySpace could not be imported: {IMPORT_ERROR}")

    satcat_client = SatCat(username, password, sessionTimeout=1800)
    try:
        start, end, epoch_window = build_query_window(epoch_range_days)
        gp_records = satcat_client._ExecuteQuery("gp", {"EPOCH": epoch_window})
    finally:
        satcat_client.Close()

    gp_rows = [normalize_record(r) for r in gp_records]
    gp_df = pd.DataFrame(gp_rows)

    satcat_records: list[dict] = []
    if not gp_df.empty:
        satcat_client = SatCat(username, password, sessionTimeout=1800)
        try:
            gp_identifiers = extract_identifier_values(
                gp_df,
                ["OBJECT_ID", "INTLDES", "INTERNATIONAL_DESIGNATOR", "INTDES", "INTLDES"],
            )

            norad_col = get_first_matching_column(gp_df, ["NORAD_CAT_ID", "NORAD_ID", "CATALOG_NUMBER", "NORADID"])
            gp_norads = []
            if norad_col is not None:
                gp_norads = [
                    str(int(float(v)))
                    for v in gp_df[norad_col].dropna().tolist()
                    if str(v).strip() and str(v).strip() not in {"nan", "None", ""}
                ]

            for identifier_list in chunk_list(gp_identifiers, 100):
                if not identifier_list:
                    continue
                try:
                    satcat_records.extend(satcat_client._ExecuteQuery("satcat", {"OBJECT_ID": identifier_list}))
                except ConnectionError:
                    continue

            for norad_list in chunk_list(gp_norads, 100):
                if not norad_list:
                    continue
                try:
                    satcat_records.extend(satcat_client._ExecuteQuery("satcat", {"NORAD_CAT_ID": norad_list}))
                except ConnectionError:
                    continue
        finally:
            satcat_client.Close()

    satcat_rows = [normalize_record(r) for r in satcat_records]
    sat_df = pd.DataFrame(satcat_rows)

    if gp_df.empty and sat_df.empty:
        return pd.DataFrame(), pd.DataFrame(), start, end

    return gp_df, sat_df, start, end


def canonicalize_df(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df

    out = df.copy()
    for col in list(out.columns):
        if col in {"index"}:
            continue
        out[col] = out[col].replace({"nan": np.nan, "None": np.nan, "null": np.nan})

    # Common SATCAT/GP naming normalization.
    rename_map = {
        "object_id": "INTLDES",
        "intl_des": "INTLDES",
        "international_designator": "INTLDES",
        "norad_cat_id": "NORAD_CAT_ID",
        "norad_id": "NORAD_CAT_ID",
        "object_name": "OBJECT_NAME",
        "satname": "OBJECT_NAME",
        "epoch": "EPOCH",
        "apogee": "APOGEE",
        "perigee": "PERIGEE",
        "semi_major_axis": "SEMIMAJOR_AXIS",
        "sma": "SEMIMAJOR_AXIS",
        "inclination": "INCLINATION",
        "inc": "INCLINATION",
        "eccentricity": "ECCENTRICITY",
        "ecc": "ECCENTRICITY",
        "object_type": "OBJECT_TYPE",
        "country_code": "COUNTRY_CODE",
    }

    out.columns = [rename_map.get(str(c).strip().lower(), str(c)) for c in out.columns]
    out["INTLDES"] = out.get("INTLDES", pd.Series(index=out.index, dtype="object")).map(normalize_intldes)
    if "NORAD_CAT_ID" in out.columns:
        out["NORAD_CAT_ID"] = pd.to_numeric(out["NORAD_CAT_ID"], errors="coerce")
    return out


def build_joined_dataset(gp_df: pd.DataFrame, sat_df: pd.DataFrame) -> pd.DataFrame:
    if gp_df.empty and sat_df.empty:
        return pd.DataFrame()

    gp = canonicalize_df(gp_df)
    sat = canonicalize_df(sat_df)

    if gp.empty:
        gp = pd.DataFrame(columns=["INTLDES", "OBJECT_NAME", "OBJECT_TYPE", "EPOCH"])
    if sat.empty:
        sat = pd.DataFrame(columns=["INTLDES", "OBJECT_NAME", "OBJECT_TYPE", "EPOCH"])

    gp = gp.copy()
    sat = sat.copy()

    gp["_match_key"] = gp.get("INTLDES", pd.Series([""] * len(gp), index=gp.index))
    sat["_match_key"] = sat.get("INTLDES", pd.Series([""] * len(sat), index=sat.index))

    gp["_match_key"] = gp["_match_key"].map(normalize_intldes)
    sat["_match_key"] = sat["_match_key"].map(normalize_intldes)

    joined = gp.merge(
        sat,
        left_on="_match_key",
        right_on="_match_key",
        how="outer",
        suffixes=("_gp", "_sat"),
    )

    if joined.empty:
        return joined

    joined["INTLDES"] = joined["_match_key"]
    for col in ["_match_key"]:
        if col in joined.columns:
            joined = joined.drop(columns=[col])

    # If the match is not on INTLDES, allow a fallback by NORAD_CAT_ID.
    if "NORAD_CAT_ID_gp" in joined.columns and "NORAD_CAT_ID_sat" in joined.columns:
        fallback_mask = joined["INTLDES"].isna() | (joined["INTLDES"] == "")
        if fallback_mask.any():
            joined.loc[fallback_mask, "INTLDES"] = joined.loc[
                fallback_mask, "NORAD_CAT_ID_gp"
            ].combine_first(joined.loc[fallback_mask, "NORAD_CAT_ID_sat"]).astype(str)

    return joined


@st.cache_data(show_spinner=False)
def load_dataset(username: str, password: str, epoch_range_days: int):
    gp_df, sat_df, start, end = fetch_space_track_rows(username, password, epoch_range_days)
    dataset = build_joined_dataset(gp_df, sat_df)
    return dataset, start, end


def infer_dtype(series: pd.Series):
    if series.empty:
        return "string"

    cleaned = series.dropna()
    if cleaned.empty:
        return "string"

    numeric = pd.to_numeric(cleaned, errors="coerce")
    if numeric.notna().sum() > 0 and numeric.notna().sum() == len(cleaned):
        if np.allclose(numeric, numeric.round()):
            return "int"
        return "float"

    bool_values = {"true", "false", "yes", "no", "y", "n", "1", "0"}
    lowered = cleaned.astype(str).str.strip().str.lower()
    if lowered.map(lambda x: x in bool_values).all():
        return "bool"

    try:
        pd.to_datetime(cleaned.astype(str), errors="raise")
        return "datetime"
    except Exception:
        return "string"


def filter_dataframe(df: pd.DataFrame, filters: dict[str, Any]) -> pd.DataFrame:
    out = df.copy()
    if out.empty:
        return out

    for column, config in filters.items():
        if not config or config.get("enabled") is False:
            continue
        dtype = infer_dtype(out[column]) if column in out.columns else None
        if dtype is None:
            continue

        if dtype in {"int", "float"}:
            op = config.get("op", "within")
            value = config.get("value")
            lower = config.get("lower")
            upper = config.get("upper")

            col_values = pd.to_numeric(out[column], errors="coerce")
            if op == ">":
                if value is not None:
                    out = out[col_values > float(value)]
            elif op == "<":
                if value is not None:
                    out = out[col_values < float(value)]
            elif op == "==":
                if value is not None:
                    out = out[col_values == float(value)]
            elif op == "within":
                if lower is not None: out = out[col_values >= float(lower)]
                if upper is not None: out = out[col_values <= float(upper)]

        elif dtype == "bool":
            choice = config.get("value", "all")
            if choice == "true":
                out = out[out[column].astype(str).str.lower().isin({"true", "1", "yes", "y"})]
            elif choice == "false":
                out = out[out[column].astype(str).str.lower().isin({"false", "0", "no", "n"})]

        elif dtype == "string":
            text = str(config.get("value", "")).strip()
            if text:
                out = out[out[column].astype(str).str.contains(text, case=False, na=False)]

        elif dtype == "datetime":
            start_value = config.get("start")
            end_value = config.get("end")
            if start_value:
                out = out[pd.to_datetime(out[column], errors="coerce") >= pd.to_datetime(start_value)]
            if end_value:
                out = out[pd.to_datetime(out[column], errors="coerce") <= pd.to_datetime(end_value)]

    return out


def get_display_columns(df: pd.DataFrame):
    preferred = [
        "INTLDES",
        "NORAD_CAT_ID",
        "OBJECT_NAME",
        "OBJECT_TYPE",
        "APOGEE",
        "PERIGEE",
        "SEMIMAJOR_AXIS",
        "INCLINATION",
        "ECCENTRICITY",
        "PERIOD",
        "EPOCH",
    ]
    existing = [c for c in preferred if c in df.columns]
    extras = [c for c in df.columns if c not in existing]
    return existing + extras


def build_sidebar_filters(df: pd.DataFrame):
    filters: dict[str, Any] = {}
    for column in sorted(df.columns):
        if column.startswith("_"):
            continue
        if column in {"INTLDES", "OBJECT_NAME", "OBJECT_TYPE"}:
            dtype = "string"
        else:
            dtype = infer_dtype(df[column])

        key = f"filter_{column}"
        if dtype in {"float", "int"}:
            numeric_values = pd.to_numeric(df[column], errors="coerce")
            valid = numeric_values.dropna()
            if valid.empty:
                continue
            lower = float(valid.min())
            upper = float(valid.max())

            if np.isclose(lower, upper):
                value = float(lower)
                filters[column] = {"enabled": True, "op": "==", "value": value}
                continue

            op = st.sidebar.selectbox(f"{column} operator", ["within", ">", "<", "=="], key=f"{key}_op")
            if op == "within":
                low, high = st.sidebar.slider(
                    f"{column} range",
                    min_value=float(lower),
                    max_value=float(upper),
                    value=(float(lower), float(upper)),
                    step=max((upper - lower) / 100, 1e-6),
                    key=f"{key}_range",
                )
                filters[column] = {"enabled": True, "op": "within", "lower": low, "upper": high}
            else:
                value = st.sidebar.number_input(
                    f"{column} value",
                    min_value=float(lower),
                    max_value=float(upper),
                    value=float(valid.median()),
                    step=max((upper - lower) / 100, 1e-6),
                    key=f"{key}_value",
                )
                filters[column] = {"enabled": True, "op": op, "value": value}

        elif dtype == "bool":
            choice = st.sidebar.selectbox(f"{column} value", ["all", "true", "false"], key=key)
            filters[column] = {"enabled": choice != "all", "value": choice}

        elif dtype == "datetime":
            start_val = st.sidebar.date_input(f"{column} start", value=pd.to_datetime(df[column].dropna().min()).date(), key=f"{key}_start")
            end_val = st.sidebar.date_input(f"{column} end", value=pd.to_datetime(df[column].dropna().max()).date(), key=f"{key}_end")
            filters[column] = {"enabled": True, "start": start_val.isoformat(), "end": end_val.isoformat()}

        else:
            value = st.sidebar.text_input(f"{column} contains", "", key=key)
            filters[column] = {"enabled": bool(value.strip()), "value": value}

    return filters


def main():
    st.set_page_config(page_title="SatFinder", page_icon="🛰️", layout="wide")
    st.title("SatFinder")
    st.caption("Fetch GP + SATCAT orbital data, match by INTLDES, then filter and inspect satellites.")

    with st.sidebar:
        st.header("Connection")
        username = st.text_input("Space-Track username", value=os.getenv("SPACETRACK_USERNAME", ""))
        password = st.text_input("Space-Track password", type="password", value=os.getenv("SPACETRACK_PASSWORD", ""))
        epoch_range_days = st.number_input("EPOCH range (days)", min_value=1, max_value=3650, value=30, step=1)
        refresh = st.button("Refresh data")

        if IMPORT_ERROR is not None:
            st.error(f"PolySpace import failed: {IMPORT_ERROR}")

    if refresh or "dataset" not in st.session_state:
        if not (username and password):
            st.warning("Provide credentials in the sidebar or set SPACETRACK_USERNAME and SPACETRACK_PASSWORD before refreshing.")
            return

        with st.spinner("Fetching GP and SATCAT records..."):
            dataset, start, end = load_dataset(username, password, int(epoch_range_days))
        st.session_state["dataset"] = dataset
        st.session_state["start"] = start
        st.session_state["end"] = end

    dataset = st.session_state.get("dataset", pd.DataFrame())
    if dataset.empty:
        st.info("No data returned for the current query. Try widening the EPOCH range or checking credentials.")
        return

    st.subheader("Data summary")
    st.write(f"Records loaded: {len(dataset)}")
    st.write(f"Window: {st.session_state.get('start')} to {st.session_state.get('end')}")

    filters = build_sidebar_filters(dataset)
    filtered = filter_dataframe(dataset, filters)

    st.sidebar.markdown("---")
    display_columns = st.sidebar.multiselect(
        "Visible columns",
        options=get_display_columns(filtered),
        default=get_display_columns(filtered)[:10],
    )

    st.subheader("Filtered results")
    if filtered.empty:
        st.warning("No satellites match the current filter set.")
        return

    st.dataframe(filtered[display_columns], use_container_width=True)
    csv_bytes = filtered[display_columns].to_csv(index=False).encode("utf-8")
    st.download_button("Download filtered CSV", data=csv_bytes, file_name="satfinder_filtered.csv", mime="text/csv")


if __name__ == "__main__":
    main()
