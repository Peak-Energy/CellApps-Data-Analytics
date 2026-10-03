"""Data in and out: raw-file locations, the Neware CSV loader, rest grouping, datapack tables.

Neware format facts this relies on (checked on every file):
- UTF-8 with BOM. The 15 °C files carry auxiliary columns whose temperature headers are duplicated,
  so columns are read by position.
- `Time` is the step-relative time printed as hh:mm:ss and floored to the second; pulses are logged
  every 0.1 s, so ten rows share each printed second. Rows sharing a second are spaced evenly in it.
- `Capacity(Ah)` is an unsigned per-step charge counter. Current is negative on discharge in the file;
  here it is flipped so discharge is positive.
- Some files contain unlogged pauses (the step clock stops, the wall clock does not); the pause is
  added back into the time axis.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- locations (read-only inputs)
ROOT = Path("/Users/dj/Documents/Projects/Battery Modeling Data/Veken/Veken-B-POR")
DATAPACK = Path("/Users/dj/cell-performance/cells/Veken170Ah-HT-B")
CACHE = Path.home() / ".cache" / "veken-b-ecm-reparam"          # parsed files; outside the repo

GITT = {15: ["CV_0010_CS2D0052_Pr_FastGITT_15C_.csv", "CV_0010_CS2D9049_Pr_FastGITT_15C_.csv"],
        25: ["CV_0010_CS2D7013_FastGITTCover20-2_25C_.csv", "CV_0010_CS2D9805_FastGITTCover20-2_25C_.csv"],
        45: ["CV_0010_CS2D2520_Pr_FastGITT_45C_.csv"]}
HPPC = {15: ["CV_0010_CS2D0052_Pr_HPPC_15C_.csv", "CV_0010_CS2D9049_Pr_HPPC_15C_.csv"],
        25: ["CV_0010_CS2D7013_Pr_HPPC_25C_.csv", "CV_0010_CS2D9805_Pr_HPPC_25C_.csv"],
        45: ["CV_0010_CS2D2520_Pr_HPPC_45C_.csv", "CV_0010_CS2D8533_Pr_HPPC_45C_ 1.csv"]}
RPT = {15: "CV_0010_CS2D7485_Pr_RPT-2_15C_.csv", 25: "CV_0010_CS2D7485_Pr_RPT_25C_.csv",
       45: "CV_0010_CS2D7485_Pr_RPT_45C_.csv", 60: "CV_0010_CS2D7485_Pr_RPT_60C.csv"}


def gitt_path(name): return ROOT / "GITT" / "Raw data" / name
def hppc_path(name): return ROOT / "HPPC" / "Raw data" / name
def rpt_path(name): return ROOT / "RPT" / name


def all_files():
    out = [("GITT", T, gitt_path(n)) for T, ns in GITT.items() for n in ns]
    out += [("HPPC", T, hppc_path(n)) for T, ns in HPPC.items() for n in ns]
    return out + [("RPT", T, rpt_path(n)) for T, n in RPT.items()]


def cell_id(path) -> str:
    m = re.search(r"CS2D\d{4}", Path(path).name)
    if not m:
        raise ValueError(f"no CS2Dxxxx token in {path}")
    return m.group(0)


# ---------------------------------------------------------------- Neware loader
CORE = {0: "datapoint", 1: "cycle", 2: "step", 3: "step_type", 4: "time_str",
        6: "current_raw", 7: "voltage", 8: "cap_ah", 14: "date"}
CORE_HEADER = ["DataPoint", "Cycle Index", "Step Index", "Step Type", "Time", "Total Time", "Current(A)",
               "Voltage(V)", "Capacity(Ah)", "Chg. Cap.(Ah)", "DChg. Cap.(Ah)", "Energy(Wh)", "Chg. Energy(Wh)",
               "DChg. Energy(Wh)", "Date", "Power(W)"]
AUX_HEADER = ["V1(V)", "V2(V)", "T1(℃)", "T2(℃)", "T1(℃)", "T2(℃)", "H2(%)", "H22(%)", "Aux. ΔV(V)", "Aux. ΔT(℃)"]
AUX_T = {18: "aux_T1a", 19: "aux_T2a", 20: "aux_T1b", 21: "aux_T2b"}   # T1a, T2a = cell probes; T1b = T2b = ambient
STEP_KIND = {"rest": "rest", "cc chg": "chg", "cc dchg": "dchg", "cp chg": "chg", "cp dchg": "dchg",
             "cv chg": "chg", "cccv chg": "chg"}
PAUSE_MIN_S = 30.0          # wall-clock jumps beyond the step clock larger than this are pauses
LOADER_VERSION = 4


def read_raw(path) -> pd.DataFrame:
    with open(path, encoding="utf-8-sig") as f:
        header = f.readline().rstrip("\r\n").split(",")
    if header[:16] != CORE_HEADER:
        raise ValueError(f"unexpected core header in {path}: {header[:16]}")
    cols = dict(CORE)
    if len(header) > 16:
        if header[16:26] != AUX_HEADER:
            raise ValueError(f"unexpected aux header in {path}: {header[16:]}")
        cols.update(AUX_T)
    idx = sorted(cols)
    df = pd.read_csv(path, encoding="utf-8-sig", header=None, skiprows=1, usecols=idx,
                     dtype={3: str, 4: str, 14: str}, engine="c")
    df.columns = [cols[i] for i in idx]
    return df


def reconstruct_time(seg, sec, kind, cap_ah, current):
    """Step-relative time with sub-second resolution: rows sharing a printed second are spaced
    evenly inside it. Also returns, for current-carrying rows, the difference to the time implied by
    the capacity counter (a check only; the counter's 1e-4 Ah resolution is too coarse at 8.5 A)."""
    n = len(sec)
    key = seg.astype(np.int64) * 10_000_000 + sec
    new_grp = np.r_[True, key[1:] != key[:-1]]
    grp = np.cumsum(new_grp) - 1
    grp_start = np.flatnonzero(new_grp)
    grp_size = np.diff(np.r_[grp_start, n])
    t = sec + (np.arange(n) - grp_start[grp]) / grp_size[grp]

    new_seg = np.r_[True, seg[1:] != seg[:-1]]
    dq = np.diff(cap_ah, prepend=0.0)
    dq[new_seg] = cap_ah[new_seg]
    i_abs = np.abs(current)
    i_prev = np.r_[i_abs[0], i_abs[:-1]]
    i_prev[new_seg] = i_abs[new_seg]
    i_mid = 0.5 * (i_abs + i_prev)
    with np.errstate(divide="ignore", invalid="ignore"):
        dt_cap = np.where(i_mid > 0, dq * 3600.0 / i_mid, 0.0)
    cs = np.cumsum(dt_cap)
    t_cap = cs - (cs - dt_cap)[np.flatnonzero(new_seg)][np.cumsum(new_seg) - 1]
    return t, np.where(kind != "rest", t_cap - t, np.nan)


def load(path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (rows, steps). Rows: t [s], current [A, discharge +], voltage, q_ah (cumulative
    discharge-positive charge), seg (step id), kind (rest/chg/dchg), mode (CC/CP/REST)."""
    df = read_raw(path)
    st = df["step_type"].str.strip().str.lower()
    unknown = set(st.unique()) - set(STEP_KIND)
    if unknown:
        raise ValueError(f"unknown step types {unknown} in {path}")
    df["kind"] = st.map(STEP_KIND)
    df["mode"] = df["step_type"].str.strip().str.split(" ").str[0].str.upper()
    df["current"] = -df.pop("current_raw")
    key = df["cycle"].to_numpy(np.int64) * 100_000 + df["step"].to_numpy(np.int64)
    df["seg"] = np.cumsum(np.r_[True, key[1:] != key[:-1]]) - 1
    seg = df["seg"].to_numpy()
    parts = df.pop("time_str").str.split(":", expand=True).astype(np.int64).to_numpy()
    sec = parts[:, 0] * 3600 + parts[:, 1] * 60 + parts[:, 2]
    df["t_rel"], df["t_cap_dev"] = reconstruct_time(seg, sec, df["kind"].to_numpy(),
                                                    df["cap_ah"].to_numpy(float), df["current"].to_numpy(float))
    # absolute time: steps chained end to start, plus unlogged pauses taken from the wall clock
    seg_end = df.groupby("seg", sort=True)["t_rel"].max().to_numpy()
    t_chain = np.r_[0.0, np.cumsum(seg_end)[:-1]][seg] + df["t_rel"].to_numpy()
    df["date"] = pd.to_datetime(df["date"], format="%Y-%m-%d %H:%M:%S")
    date_s = (df["date"] - df["date"].iloc[0]).dt.total_seconds().to_numpy()
    excess = np.diff(date_s) - np.diff(t_chain)
    df["pause_s"] = np.r_[0.0, np.where(excess > PAUSE_MIN_S, excess, 0.0)]
    df["t"] = t_chain + np.cumsum(df["pause_s"].to_numpy())
    # cumulative charge from the per-step counter, discharge positive
    sign = np.select([df["kind"] == "dchg", df["kind"] == "chg"], [1.0, -1.0], 0.0)
    seg_q = df.groupby("seg", sort=True)["cap_ah"].last().to_numpy() * pd.Series(sign).groupby(seg).first().to_numpy()
    df["q_ah"] = np.r_[0.0, np.cumsum(seg_q)[:-1]][seg] + sign * df["cap_ah"].to_numpy()

    g = df.groupby("seg", sort=True)
    steps = pd.DataFrame({
        "step_type": g["step_type"].first(), "kind": g["kind"].first(), "mode": g["mode"].first(),
        "n": g.size(), "t0": g["t"].first(), "dur_s": g["t_rel"].max(),
        "dt_median_s": g["t"].apply(lambda s: float(np.median(np.diff(s.to_numpy()))) if len(s) > 1 else np.nan),
        "I_mean": g["current"].mean(), "V0": g["voltage"].first(), "V_end": g["voltage"].last(),
        "cap_ah": g["cap_ah"].last(), "q_end": g["q_ah"].last(), "pause_s": g["pause_s"].sum(),
    })
    return df, steps


def load_cached(path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """load() with a pickle cache in CACHE, keyed on file name, size, mtime and loader version."""
    p = Path(path)
    st = p.stat()
    tag = f"{p.stem.replace(' ', '_')}_{st.st_size}_{int(st.st_mtime)}_v{LOADER_VERSION}"
    CACHE.mkdir(parents=True, exist_ok=True)
    f_s, f_g = CACHE / f"{tag}.rows.pkl", CACHE / f"{tag}.steps.pkl"
    if f_s.exists() and f_g.exists():
        return pd.read_pickle(f_s), pd.read_pickle(f_g)
    df, steps = load(p)
    df.to_pickle(f_s)
    steps.to_pickle(f_g)
    return df, steps


# ---------------------------------------------------------------- rests
def merged_rests(df: pd.DataFrame, steps: pd.DataFrame) -> pd.DataFrame:
    """One row per rest (consecutive Rest steps merged: Neware splits some rests at cycle
    boundaries), with the current step that preceded it."""
    is_rest = (steps.kind == "rest").to_numpy()
    sid = steps.index.to_numpy()
    rows, i = [], 0
    while i < len(sid):
        if not is_rest[i]:
            i += 1
            continue
        j = i
        while j + 1 < len(sid) and is_rest[j + 1]:
            j += 1
        prev = steps.iloc[i - 1] if i > 0 else None
        first, last = sid[i], sid[j]
        rows.append(dict(seg_first=first, seg_last=last, prev_seg=sid[i - 1] if i > 0 else -1,
                         prev_kind=prev.kind if prev is not None else "none",
                         prev_I=prev.I_mean if prev is not None else 0.0,
                         prev_dur=prev.dur_s if prev is not None else 0.0,
                         prev_cap=prev.cap_ah if prev is not None else 0.0,
                         prev_Vend=prev.V_end if prev is not None else np.nan,
                         t0=steps.loc[first, "t0"], dur_s=steps.loc[first:last, "dur_s"].sum(),
                         q_ah=steps.loc[last, "q_end"], V_end=steps.loc[last, "V_end"]))
        i = j + 1
    return pd.DataFrame(rows)


def rest_trace(df, r):
    """(t since rest start [s], V, rows) for a merged rest."""
    d = df[(df.seg >= r.seg_first) & (df.seg <= r.seg_last)]
    return d.t.to_numpy() - d.t.iloc[0], d.voltage.to_numpy(), d


# ---------------------------------------------------------------- datapack tables
OCV_COLS = ["Temperature", "SOC", "OCV-Charge", "OCV-Discharge", "Average"]
R_COLS = ["Temperature", "SOC", "Mean (mOhms)", "Error (mOhms)"]
C_COLS = ["Temperature", "SOC", "Mean (Farad)", "Error (Farad)"]
ENTROPY_COLS = ["SOC", "Charge", "Discharge"]
R_SOC_GRID = [0, 0.025, 0.05, 0.075, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6,
              0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.925, 0.95, 0.975, 1]
R_T_GRID = [-30, -20, -10, 0, 10, 15, 25, 45]


def read_table(path) -> pd.DataFrame:
    return pd.read_csv(path, encoding="utf-8-sig")


def write_table(df: pd.DataFrame, path, kind: str) -> None:
    """Write in the datapack format: LF, no BOM, sorted by Temperature then SOC, SOC as a fraction,
    Temperature as an integer. kind: 'ocv' | 'r' | 'c' | 'entropy'."""
    cols = {"ocv": OCV_COLS, "r": R_COLS, "c": C_COLS, "entropy": ENTROPY_COLS}[kind]
    if list(df.columns) != cols:
        raise ValueError(f"{kind}: columns {list(df.columns)} != {cols}")
    if df.isna().any().any():
        raise ValueError(f"{kind}: NaN in table")
    df = df.sort_values(["SOC"] if kind == "entropy" else ["Temperature", "SOC"], kind="stable")
    fmt = {"ocv": lambda x: f"{float(x):.6f}".rstrip("0").rstrip(".") or "0",
           "r": lambda x: f"{float(x):.9g}", "c": lambda x: f"{float(x):.9g}",
           "entropy": lambda x: f"{float(x):.3E}"}[kind]
    lines = [",".join(cols)]
    for row in df.itertuples(index=False):
        out = []
        for c, v in zip(cols, row):
            out.append(str(int(round(v))) if c == "Temperature" else f"{round(float(v), 4):g}" if c == "SOC" else fmt(v))
        lines.append(",".join(out))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
