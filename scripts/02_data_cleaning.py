"""
02_data_cleaning.py
Wrangle AIHW PBS Excel workbooks and ABS linkage files into tidy panel datasets.

Inputs (from data/raw/):
    AIHW PBS supplementary data tables:
      - AIHW-HWE-098-PBS-ATC2-prescriptions-historical-data.xlsx (May 2002 – Jan 2015)
      - AIHW-HWE-098-PBS-ATC2-prescriptions-monthly-data.xlsx   (Jan 2015 – Dec 2025)
    ABS linkage files:
      - SEIFA_2021_LGA.xlsx
      - RA_2021_AUST.xlsx          (SA1 → Remoteness Area)
      - LGA_2021_AUST_allocation.xlsx  (MB → LGA)
      - MB_2021_AUST.xlsx          (MB → SA1, for crosswalk)
      - ERP_LGA_2001_2024.xlsx

Outputs (to data/processed/):
    - pbs_national_monthly.csv   : month × ATC2 × age_group × metric
    - pbs_state_monthly.csv      : month × state × ATC2 × age_group × metric
    - pbs_lga_monthly.csv        : month × LGA × ATC2 × age_group × metric
    - lga_characteristics.csv    : LGA-level SEIFA quintile + remoteness
    - erp_lga.csv                : annual ERP by LGA (2001–2024)

Usage:
    source .venv/bin/activate
    python scripts/02_data_cleaning.py
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=UserWarning, module="openpyxl")

# -- Paths --------------------------------------------------------------------
PROJ_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = PROJ_DIR / "data" / "raw"
PROC_DIR = PROJ_DIR / "data" / "processed"
PROC_DIR.mkdir(parents=True, exist_ok=True)

INTERVENTION_DATE = pd.Timestamp("2023-01-01")


# =============================================================================
# SECTION 1: SEIFA
# =============================================================================

def load_seifa():
    """Load ABS SEIFA 2021 IRSD at LGA level from Table 1."""
    path = RAW_DIR / "SEIFA_2021_LGA.xlsx"
    if not path.exists():
        print("WARNING: SEIFA_2021_LGA.xlsx not found.")
        return pd.DataFrame()

    print("\nLoading SEIFA 2021...")

    # Table 1: rows 0–4 metadata, row 5 headers, row 6+ data.
    # Fixed layout: col 0=LGA code, 1=LGA name, 2=IRSD Score, 3=IRSD Decile,
    # 4=IRSAD Score, 5=IRSAD Decile, 6=IER Score, 7=IER Decile,
    # 8=IEO Score, 9=IEO Decile, 10=Usual Resident Population
    raw = pd.read_excel(path, sheet_name="Table 1", header=None, skiprows=6)
    raw.columns = [
        "lga_code", "lga_name",
        "irsd_score", "irsd_decile", "irsad_score", "irsad_decile",
        "ier_score", "ier_decile", "ieo_score", "ieo_decile",
        "usual_resident_pop_2021",
    ]

    raw["lga_code"] = (raw["lga_code"].astype(str).str.strip()
                       .str.replace(r"\.0$", "", regex=True))
    raw = raw[raw["lga_code"].str.match(r"^\d+$")].copy()

    for col in ["irsd_score", "irsd_decile", "usual_resident_pop_2021"]:
        raw[col] = pd.to_numeric(raw[col], errors="coerce")

    raw = raw.dropna(subset=["irsd_score"])
    raw["seifa_quintile"] = np.ceil(raw["irsd_decile"] / 2).astype("Int64")

    seifa = raw[["lga_code", "lga_name", "irsd_score", "irsd_decile",
                 "seifa_quintile", "usual_resident_pop_2021"]].copy()

    print(f"  Loaded {len(seifa)} LGAs")
    print(f"  SEIFA quintile distribution:")
    print(seifa["seifa_quintile"].value_counts().sort_index().to_string())
    return seifa


# =============================================================================
# SECTION 2: Remoteness (MB → SA1 → RA, MB → LGA)
# =============================================================================

def load_remoteness():
    """Build LGA-to-Remoteness crosswalk via Mesh Block intermediary.

    Chain: MB_2021_AUST.xlsx (MB→SA1) + RA_2021_AUST.xlsx (SA1→RA)
           + LGA_2021_AUST_allocation.xlsx (MB→LGA)
    Then assign each LGA its modal Remoteness Area.
    """
    ra_path = RAW_DIR / "RA_2021_AUST.xlsx"
    lga_path = RAW_DIR / "LGA_2021_AUST_allocation.xlsx"
    mb_path = RAW_DIR / "MB_2021_AUST.xlsx"

    for p, name in [(ra_path, "RA_2021_AUST.xlsx"),
                    (lga_path, "LGA_2021_AUST_allocation.xlsx"),
                    (mb_path, "MB_2021_AUST.xlsx")]:
        if not p.exists():
            print(f"WARNING: {name} not found. Cannot build remoteness crosswalk.")
            return pd.DataFrame()

    print("\nBuilding LGA-to-Remoteness crosswalk...")

    # 1. MB → SA1 (from main structure allocation)
    mb = pd.read_excel(mb_path, sheet_name="MB_2021_AUST",
                       usecols=["MB_CODE_2021", "SA1_CODE_2021"],
                       dtype=str)
    mb["MB_CODE_2021"] = mb["MB_CODE_2021"].str.strip()
    mb["SA1_CODE_2021"] = mb["SA1_CODE_2021"].str.strip()
    print(f"  MB→SA1: {len(mb):,} Mesh Blocks")

    # 2. SA1 → RA
    ra = pd.read_excel(ra_path, sheet_name="SA1_RA_2021_AUST",
                       usecols=["SA1_CODE_2021", "RA_CODE_2021", "RA_NAME_2021"],
                       dtype=str)
    ra["SA1_CODE_2021"] = ra["SA1_CODE_2021"].str.strip()
    print(f"  SA1→RA: {len(ra):,} SA1s")

    # 3. MB → LGA
    lga = pd.read_excel(lga_path, sheet_name="LGA_2021_AUST",
                        usecols=["MB_CODE_2021", "LGA_CODE_2021"],
                        dtype=str)
    lga["MB_CODE_2021"] = lga["MB_CODE_2021"].str.strip()
    lga["LGA_CODE_2021"] = lga["LGA_CODE_2021"].str.strip()
    print(f"  MB→LGA: {len(lga):,} Mesh Blocks → {lga['LGA_CODE_2021'].nunique()} LGAs")

    # 4. Join: MB → (SA1, LGA) → RA
    mb_sa1_lga = mb.merge(lga, on="MB_CODE_2021", how="inner")
    mb_sa1_lga_ra = mb_sa1_lga.merge(ra, on="SA1_CODE_2021", how="inner")

    # 5. Modal remoteness per LGA
    lga_ra = (mb_sa1_lga_ra
              .groupby("LGA_CODE_2021")["RA_NAME_2021"]
              .agg(lambda x: x.mode().iloc[0] if len(x.mode()) > 0 else np.nan)
              .reset_index())
    lga_ra.columns = ["lga_code", "remoteness_area"]

    ra_order = {
        "Major Cities of Australia": 1,
        "Inner Regional Australia": 2,
        "Outer Regional Australia": 3,
        "Remote Australia": 4,
        "Very Remote Australia": 5,
    }
    lga_ra["remoteness_code"] = lga_ra["remoteness_area"].map(ra_order)

    # Filter out non-LGA codes and non-standard remoteness categories
    lga_ra = lga_ra[lga_ra["lga_code"].str.match(r"^\d+$")]
    lga_ra = lga_ra[lga_ra["remoteness_code"].notna()]

    print(f"  Crosswalk: {len(lga_ra)} LGAs")
    print(f"  Distribution:")
    print(lga_ra["remoteness_area"].value_counts().to_string())
    return lga_ra


# =============================================================================
# SECTION 3: ERP
# =============================================================================

def load_erp():
    """Load ABS ERP by LGA, 2001–2024 (32180DS0004, Table 1)."""
    path = RAW_DIR / "ERP_LGA_2001_2024.xlsx"
    if not path.exists():
        print("WARNING: ERP_LGA_2001_2024.xlsx not found.")
        return pd.DataFrame()

    print("\nLoading ERP by LGA...")

    # Multi-row header: row 5 = years (2001–2024), row 6 = labels, data from row 7
    raw = pd.read_excel(path, sheet_name="Table 1", header=None)
    years = raw.iloc[5].tolist()

    col_names = ["lga_code", "lga_name"]
    for i, yr in enumerate(years):
        if i < 2:
            continue
        try:
            col_names.append(str(int(yr)))
        except (ValueError, TypeError):
            col_names.append(f"col_{i}")

    data = raw.iloc[7:].copy()
    data.columns = col_names
    data["lga_code"] = (data["lga_code"].astype(str).str.strip()
                        .str.replace(r"\.0$", "", regex=True))
    data = data[data["lga_code"].str.match(r"^\d+$")]

    year_cols = [c for c in data.columns if re.match(r"^20\d{2}$", c)]
    erp_long = data.melt(id_vars=["lga_code", "lga_name"],
                         value_vars=year_cols,
                         var_name="year", value_name="erp")
    erp_long["year"] = erp_long["year"].astype(int)
    erp_long["erp"] = pd.to_numeric(erp_long["erp"], errors="coerce")
    erp_long = erp_long.dropna(subset=["erp"])

    erp_out = erp_long[["lga_code", "year", "erp"]].copy()
    print(f"  Loaded ERP for {erp_out['lga_code'].nunique()} LGAs, "
          f"years {erp_out['year'].min()}–{erp_out['year'].max()}")
    return erp_out


# =============================================================================
# SECTION 4: LGA characteristics
# =============================================================================

def build_lga_characteristics(seifa, remoteness):
    """Merge SEIFA and remoteness into a single LGA lookup."""
    print("\nMerging LGA characteristics...")
    if seifa.empty and remoteness.empty:
        return pd.DataFrame()
    if seifa.empty:
        return remoteness
    if remoteness.empty:
        return seifa

    lga_chars = seifa.merge(remoteness, on="lga_code", how="outer")
    print(f"  {len(lga_chars)} LGAs total")
    print(f"  SEIFA coverage: {lga_chars['seifa_quintile'].notna().sum()}")
    print(f"  Remoteness coverage: {lga_chars['remoteness_area'].notna().sum()}")
    return lga_chars


# =============================================================================
# SECTION 5: AIHW PBS data parser
# =============================================================================
# AIHW ATC2 workbook structure (verified):
#   Table 1: National — cols: Type of script, Age group, Value, MAY-2002...
#   Table 2: State   — cols: State, Type of script, Age group, Value, MAY-2002...
#   Table 3: LGA 0–18  — cols: State, LGA code, LGA name, Type of script, Value, MAY-2002...
#   Table 4: LGA 19–64 — cols: same
#   Table 5: LGA 65+   — cols: same
#   Table 6: LGA Total  — cols: same
# "Value" has two levels: "Benefits per ERP", "Scripts per ERP"
# Months as columns (wide format). Melt to long.
# =============================================================================

AGE_GROUP_BY_TABLE = {
    "Table 3": "0-18",
    "Table 4": "19-64",
    "Table 5": "65+",
    "Table 6": "Total",
}

MONTH_PATTERN = re.compile(
    r"^(JAN|FEB|MAR|APR|MAY|JUN|JUL|AUG|SEP|OCT|NOV|DEC)-(\d{4})$"
)


def parse_pbs_wide_sheet(filepath, sheet_name):
    """Parse one sheet of an AIHW PBS workbook from wide to long format."""
    raw = pd.read_excel(filepath, sheet_name=sheet_name, header=None,
                        dtype=str)
    if raw.shape[0] < 3:
        return pd.DataFrame()

    # Row 0 = title, Row 1 = headers, Row 2+ = data
    headers = raw.iloc[1].tolist()
    data = raw.iloc[2:].copy()
    data.columns = headers

    # Drop empty rows and footer rows
    data = data.dropna(how="all")
    first_col = headers[0]
    data = data[~data[first_col].astype(str).str.match(
        r"(?i)^(source|note|\.\.|nan)?$", na=True)]

    # Identify month columns vs ID columns
    month_cols = [c for c in headers if isinstance(c, str) and MONTH_PATTERN.match(c)]
    id_cols = [c for c in headers if c not in month_cols]

    if not month_cols:
        print(f"    No month columns found in {sheet_name}")
        return pd.DataFrame()

    # Rename "Value" ID column to "metric" BEFORE melting to avoid collision
    # with the melted numeric values column
    rename_before = {c: "metric" for c in id_cols if str(c).strip().lower() == "value"}
    data = data.rename(columns=rename_before)
    id_cols = ["metric" if str(c).strip().lower() == "value" else c for c in id_cols]

    # Melt wide → long
    long = data.melt(id_vars=id_cols, value_vars=month_cols,
                     var_name="month_str", value_name="rate")

    # Parse date
    long["date"] = pd.to_datetime(long["month_str"], format="%b-%Y")
    long["rate"] = pd.to_numeric(long["rate"], errors="coerce")
    long = long.dropna(subset=["rate"])
    long = long.drop(columns=["month_str"])

    return long


def _rename_pbs_cols(df, extra_renames=None):
    """Standardise AIHW PBS column names."""
    rename = {}
    for c in df.columns:
        cl = str(c).lower()
        if "type" in cl and "script" in cl:
            rename[c] = "atc_name"
        elif "age" in cl and "group" in cl:
            rename[c] = "age_group"
        elif "state" in cl:
            rename[c] = "state"
        elif "lga" in cl and "code" in cl:
            rename[c] = "lga_code"
        elif "lga" in cl and "name" in cl:
            rename[c] = "lga_name"
    if extra_renames:
        rename.update(extra_renames)
    return df.rename(columns=rename)


def parse_pbs_national(filepath):
    """Parse Table 1 (national level) from a PBS ATC workbook."""
    print(f"  Parsing Table 1 (national)...")
    df = parse_pbs_wide_sheet(filepath, "Table 1")
    if df.empty:
        return df
    df = _rename_pbs_cols(df)
    print(f"    {len(df):,} rows")
    return df


def parse_pbs_state(filepath):
    """Parse Table 2 (state level) from a PBS ATC workbook."""
    print(f"  Parsing Table 2 (state)...")
    df = parse_pbs_wide_sheet(filepath, "Table 2")
    if df.empty:
        return df
    df = _rename_pbs_cols(df)
    print(f"    {len(df):,} rows")
    return df


def parse_pbs_lga(filepath):
    """Parse Tables 3–6 (LGA level, one per age group)."""
    all_dfs = []
    for table, age_group in AGE_GROUP_BY_TABLE.items():
        print(f"  Parsing {table} (LGA, age {age_group})...")
        df = parse_pbs_wide_sheet(filepath, table)
        if df.empty:
            continue

        df = _rename_pbs_cols(df)
        df["age_group"] = age_group
        all_dfs.append(df)
        print(f"    {len(df):,} rows")

    if all_dfs:
        return pd.concat(all_dfs, ignore_index=True)
    return pd.DataFrame()


def load_all_pbs():
    """Load and combine ATC2 historical + monthly PBS data."""
    # Find ATC2 files (prioritise these for the analysis)
    atc2_files = sorted(RAW_DIR.glob("*ATC2*prescriptions*.xlsx"))
    if not atc2_files:
        # Fall back to any PBS files
        atc2_files = sorted(RAW_DIR.glob("*PBS*ATC*prescriptions*.xlsx"))

    if not atc2_files:
        print("WARNING: No PBS ATC2 data files found in data/raw/")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    print(f"\nFound {len(atc2_files)} ATC2 PBS file(s):")
    for f in atc2_files:
        print(f"  {f.name} ({f.stat().st_size / 1e6:.0f} MB)")

    all_national, all_state, all_lga = [], [], []

    for fpath in atc2_files:
        print(f"\n--- {fpath.name} ---")

        nat = parse_pbs_national(fpath)
        if not nat.empty:
            all_national.append(nat)

        st = parse_pbs_state(fpath)
        if not st.empty:
            all_state.append(st)

        lga = parse_pbs_lga(fpath)
        if not lga.empty:
            all_lga.append(lga)

    national = pd.concat(all_national, ignore_index=True) if all_national else pd.DataFrame()
    state = pd.concat(all_state, ignore_index=True) if all_state else pd.DataFrame()
    lga = pd.concat(all_lga, ignore_index=True) if all_lga else pd.DataFrame()

    # De-duplicate overlapping dates (Jan 2015 may appear in both files)
    for df in [national, state, lga]:
        if not df.empty and "date" in df.columns:
            group_cols = [c for c in df.columns if c not in ["rate", "date"]]
            df.drop_duplicates(subset=group_cols + ["date"], keep="last", inplace=True)

    return national, state, lga


# =============================================================================
# SECTION 6: Pivot metrics (Scripts per ERP / Benefits per ERP) to columns
# =============================================================================

def pivot_metrics(df):
    """Pivot 'metric' column so Scripts per ERP and Benefits per ERP become
    separate columns instead of rows. Halves row count, easier to analyse."""
    if df.empty or "metric" not in df.columns:
        return df

    id_cols = [c for c in df.columns if c not in ["metric", "rate"]]
    pivoted = df.pivot_table(index=id_cols, columns="metric",
                             values="rate", aggfunc="first").reset_index()
    pivoted.columns.name = None

    # Clean column names
    rename = {}
    for c in pivoted.columns:
        cl = str(c).lower()
        if "script" in cl and "erp" in cl:
            rename[c] = "scripts_per_erp"
        elif "benefit" in cl and "erp" in cl:
            rename[c] = "benefits_per_erp"
    pivoted = pivoted.rename(columns=rename)

    return pivoted


# =============================================================================
# SECTION 7: Add ITS analysis variables
# =============================================================================

def add_ts_variables(df):
    """Add time series variables for ITS modelling."""
    if "date" not in df.columns:
        return df

    df = df.sort_values("date").copy()
    min_date = df["date"].min()

    df["time_index"] = ((df["date"] - min_date).dt.days / 30.44).round().astype(int)
    df["post_reform"] = (df["date"] >= INTERVENTION_DATE).astype(int)
    reform_months = ((df["date"] - INTERVENTION_DATE).dt.days / 30.44).round()
    df["time_post_reform"] = np.where(
        df["post_reform"] == 1, reform_months, 0
    ).astype(int)

    df["year"] = df["date"].dt.year
    df["month_num"] = df["date"].dt.month
    df["quarter"] = df["date"].dt.quarter
    df["covid_period"] = (
        (df["date"] >= "2020-03-01") & (df["date"] <= "2021-12-31")
    ).astype(int)

    return df


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("PBS Copayment Reform — Data Cleaning Pipeline")
    print("=" * 60)

    # ---- ABS linkage files ----
    seifa = load_seifa()
    remoteness = load_remoteness()
    erp = load_erp()
    lga_chars = build_lga_characteristics(seifa, remoteness)

    if not lga_chars.empty:
        out = PROC_DIR / "lga_characteristics.csv"
        lga_chars.to_csv(out, index=False)
        print(f"\nSaved: {out.name} ({len(lga_chars)} LGAs)")

    if not erp.empty:
        out = PROC_DIR / "erp_lga.csv"
        erp.to_csv(out, index=False)
        print(f"Saved: {out.name} ({len(erp):,} rows)")

    # ---- PBS data ----
    national, state, lga = load_all_pbs()

    if national.empty and state.empty and lga.empty:
        print("\nNo PBS data parsed. Download AIHW files to data/raw/ first.")
        return

    # Pivot metrics to columns and add time series variables
    for name, df in [("national", national), ("state", state), ("lga", lga)]:
        if df.empty:
            continue

        df = pivot_metrics(df)
        df = add_ts_variables(df)

        # Merge LGA data with characteristics
        if name == "lga" and not lga_chars.empty and "lga_code" in df.columns:
            df["lga_code"] = (df["lga_code"].astype(str).str.strip()
                              .str.replace(r"\.0$", "", regex=True))
            df = df.merge(
                lga_chars[["lga_code", "seifa_quintile", "remoteness_area",
                           "remoteness_code"]],
                on="lga_code", how="left"
            )

        out = PROC_DIR / f"pbs_{name}_monthly.csv"
        df.to_csv(out, index=False)
        print(f"\nSaved: {out.name}")
        print(f"  Rows: {len(df):,}")
        print(f"  Date range: {df['date'].min().strftime('%b %Y')} – "
              f"{df['date'].max().strftime('%b %Y')}")
        if "atc_name" in df.columns:
            print(f"  ATC classes: {df['atc_name'].nunique()}")
        if "lga_code" in df.columns:
            print(f"  LGAs: {df['lga_code'].nunique()}")

    # ---- Summary ----
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE")
    print("=" * 60)
    print(f"\nOutputs in: {PROC_DIR}")
    for f in sorted(PROC_DIR.glob("*.csv")):
        print(f"  {f.name} ({f.stat().st_size / 1e6:.1f} MB)")
    print(f"\nIntervention: {INTERVENTION_DATE.strftime('%B %Y')}")
    print("\nIMPORTANT: AIHW data provides total prescriptions only (no")
    print("general/concessional split). For CausalImpact, consider using")
    print("age 65+ (mostly concessional) vs 19–64 (mostly general) as proxy.")
    print("\nNext step: Run 03_its_analysis.R")


if __name__ == "__main__":
    main()
