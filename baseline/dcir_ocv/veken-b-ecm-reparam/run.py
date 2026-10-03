"""Rebuild everything:  python3 run.py

1. New OCV from GITT (+ HPPC rests).
2. 1RC from HPPC, fitted twice: with the new OCV and with the team's current OCV subtracted.
3. Two drop-in cell folders (same files and cell.yaml format as cell-performance/cells/Veken170Ah-HT-B):
     outputs/<run_date>/new_ocv/Veken170Ah-HT-B/   new ocv.csv + r0/r1/c1 fitted on it
     outputs/<run_date>/team_ocv/Veken170Ah-HT-B/  team ocv.csv (unchanged) + r0/r1/c1 fitted on it
   Every other file is copied unchanged from the team datapack.
4. Whole-file simulation of every RPT and HPPC record with the team datapack and both sets; figures;
   provenance.json. Prints the summary tables quoted in README.md.
"""
import filecmp
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from reparam import ecm, figures, load, ocv, validate

HERE = Path(__file__).resolve().parent
cfg = yaml.safe_load(open(HERE / "config.yaml"))
DATE = cfg["run_date"]
OUT = HERE / "outputs" / DATE
if OUT.exists():
    shutil.rmtree(OUT)
FIG = OUT / "figures"
FIG.mkdir(parents=True)
Q = {int(k): v for k, v in cfg["soc"]["capacity_ah"].items()}
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 200)
team = {k: load.read_table(load.DATAPACK / f"{k}.csv") for k in ["ocv", "r0", "r1", "c1", "entropy"]}

# ---------------------------------------------------------------- 1-2. OCV and the two 1RC fits
print("1/4 new OCV ...", flush=True)
ocv_new, ent_candidate, ocv_info = ocv.build(cfg)
print("2/4 1RC fits ...", flush=True)
sets = {}
for name, ocv_tab, leave_out in [("new OCV set", ocv_new, []),
                                 ("team OCV set", team["ocv"], cfg["hppc"]["team_ocv_leave_out"])]:
    fits, tau1, axes = ecm.fit_hppc(cfg, ocv_tab, leave_out)
    full, per_cell, ea, fallback = ecm.tables(fits)
    r0, r1, c1 = ecm.to_datapack(full)
    sets[name] = dict(ocv=ocv_tab, r0=r0, r1=r1, c1=c1, fits=fits, tau1=tau1, axes=axes, per_cell=per_cell,
                      ea=ea, fallback=fallback)

# ---------------------------------------------------------------- 3. cell folders
OFFSET_AH = {}
for T, names in load.GITT.items():
    for name in names:
        df, steps = load.load_cached(load.gitt_path(name))
        d = float(ocv.gitt_branches(df, steps, "cp")["anchor100"].q_ah - ocv.gitt_branches(df, steps, "c20")["anchor100"].q_ah)
        OFFSET_AH.setdefault(T, []).append(d)
OFFSET_AH = {T: float(np.mean(v)) for T, v in OFFSET_AH.items()}

ORIGIN_HPPC = ("-B Veken 170Ah-HT HPPC CS2D0052, CS2D9049 (15C), CS2D7013, CS2D9805 (25C), CS2D2520, CS2D8533 (45C); "
               "reparameterized in CellApps-Data-Analytics/baseline/dcir_ocv/veken-b-ecm-reparam")


def block(key, description, origin):
    lines = [f"  {key}:", f"    file: {key}.csv", "    description: >"]
    words, line = description.split(), "     "
    for w in words:
        if len(line) + len(w) + 1 > 100:
            lines.append(line); line = "     "
        line += " " + w
    lines += [line, "    source:", f"      origin: {origin}", f'      retrieved: "{DATE}"']
    return "\n".join(lines) + "\n"


def replace_block(text, key, new):
    pat = re.compile(rf"  {key}:\n    file: {key}\.csv\n.*?retrieved: \"[^\"]*\"\n", re.S)
    if not pat.search(text):
        raise ValueError(f"{key} block not found in cell.yaml")
    return pat.sub(lambda m: new, text, count=1)


def rc_blocks(text, set_name, tau1):
    on = "the new OCV (ocv.csv in this folder)" if set_name == "new OCV set" else "the team ocv.csv (unchanged, in this folder)"
    text = replace_block(text, "r0", block("r0",
        "Ohmic resistance: voltage change at the first complete sample (0.1 s) of the HPPC discharge pulses, median over "
        "pulse currents. 15/25/45 C measured (Mean = median of 2 cells, Error = std of 2 cells); -30 to 10 C per-SOC "
        f"Arrhenius through the measured rows. Fitted with {on} subtracted.", ORIGIN_HPPC))
    text = replace_block(text, "r1", block("r1",
        f"1RC resistance of the minutes-scale response, from HPPC SOC moves + 2 h rests with one time constant for all "
        f"SOC and temperatures (tau1 = {tau1:.0f} s), R0 held. Rows as for r0. Fitted with {on} subtracted."
        + ("" if set_name == "new OCV set" else " The 45 C move ending at 20 % SOC is left out (the cells' rested voltage "
           "at 25 % sits 9-10 mV above this OCV); 20 % is interpolated from its neighbours."), ORIGIN_HPPC))
    text = replace_block(text, "c1", block("c1", f"C1 = tau1 / R1 with tau1 = {tau1:.0f} s. Fitted with {on} subtracted.",
                                           ORIGIN_HPPC))
    return text


yaml_src = (load.DATAPACK / "cell.yaml").read_text()
status_old = ("# -B DATAPACK STATUS (2026-09-09): nominal values are MEASURED -B data\n"
              "# (mass / specific heat / capacity); OCV (RPT, 2026-08-04) and R0 / R1 / C1\n"
              "# (CAPE-215 HPPC, 2026-05-20) are -B files. ")
assert status_old in yaml_src, "cell.yaml status comment changed upstream"
for set_name, folder in [("new OCV set", "new_ocv"), ("team OCV set", "team_ocv")]:
    S = sets[set_name]
    dest = OUT / folder / "Veken170Ah-HT-B"
    shutil.copytree(load.DATAPACK, dest, ignore=shutil.ignore_patterns("Claude outputs", ".DS_Store"))
    for key, kind in [("r0", "r"), ("r1", "r"), ("c1", "c")]:
        load.write_table(S[key], dest / f"{key}.csv", kind)
    text = rc_blocks(yaml_src, set_name, S["tau1"])
    if set_name == "new OCV set":
        load.write_table(S["ocv"], dest / "ocv.csv", "ocv")
        text = replace_block(text, "ocv", block("ocv",
            "OCV (charge, discharge, average) vs SOC at 15, 25 and 45 C in 1 % SOC steps: average of the GITT rest "
            "samples drifting less than 1 mV/h; gaps filled with HPPC 2 h rest voltages of the same cells (15/45 C) or "
            "the 15/45 C curve shape (25 C). SOC = 1 - dQ/Q(T); 100 % = rested after a slow (C/20) charge to 3.45 V; "
            "0 % = rested after a normal-rate discharge to 1.5 V. A normal-rate (~0.24C) charge to 3.45 V ends at SOC "
            f"{1 - OFFSET_AH[15] / Q[15]:.3f} / {1 - OFFSET_AH[25] / Q[25]:.3f} / {1 - OFFSET_AH[45] / Q[45]:.3f} at "
            "15 / 25 / 45 C on this axis.",
            "-B Veken 170Ah-HT GITT CS2D0052, CS2D9049 (15C), CS2D7013, CS2D9805 (25C), CS2D2520 (45C) + HPPC 2 h rests "
            "of the same cells; reparameterized in CellApps-Data-Analytics/baseline/dcir_ocv/veken-b-ecm-reparam"))
        status = (f"# -B DATAPACK STATUS ({DATE}): nominal values are MEASURED -B data\n"
                  f"# (mass / specific heat / capacity); OCV (GITT reparam, {DATE}) and R0 / R1 / C1\n"
                  f"# (HPPC reparam fitted on that OCV, {DATE}) are -B files. ")
    else:
        status = (f"# -B DATAPACK STATUS ({DATE}): nominal values are MEASURED -B data\n"
                  f"# (mass / specific heat / capacity); OCV (unchanged, 2026-08-04) and R0 / R1 / C1\n"
                  f"# (HPPC reparam fitted on that OCV, {DATE}) are -B files. ")
    text = text.replace(status_old, status)
    yaml.safe_load(text)                                       # still valid YAML
    (dest / "cell.yaml").write_text(text)
    changed = {f.name for f in dest.iterdir() if not filecmp.cmp(f, load.DATAPACK / f.name, shallow=False)}
    expected = {"r0.csv", "r1.csv", "c1.csv", "cell.yaml"} | ({"ocv.csv"} if set_name == "new OCV set" else set())
    assert changed == expected, f"{folder}: unexpected changed files {changed ^ expected}"

# ---------------------------------------------------------------- 4. validation, figures, provenance
print("3/4 whole-file simulations ...", flush=True)
team_rc = ecm.from_datapack(team["r0"], team["r1"], team["c1"])
sim = {"team datapack": (team["ocv"], team_rc)}
sim.update({n: (S["ocv"], ecm.from_datapack(S["r0"], S["r1"], S["c1"])) for n, S in sets.items()})
rows, rpt_traces, hppc_traces, hppc_file, jumps = [], {}, {}, {}, []
for T, name in load.RPT.items():
    starts = {"new OCV set": {"SOC 1 at full charge": 1.0, "SOC at a normal-rate charge": 1.0 - OFFSET_AH.get(T, OFFSET_AH[45]) / Q[T]}}
    r, tr, j = validate.run_file("RPT", T, load.rpt_path(name), sim, Q[T], starts)
    rows += r; rpt_traces[T] = tr; jumps.append(j.assign(test="RPT", T=T, file=name))
for T, names in load.HPPC.items():
    for name in names:
        label = f"HPPC {T} °C, cell {load.cell_id(name)}"
        r, tr, j = validate.run_file("HPPC", T, load.hppc_path(name), sim, Q[T], {})
        rows += r; hppc_traces[label] = tr; hppc_file[label] = name
        jumps.append(j.assign(test="HPPC", T=T, file=name))
val = pd.DataFrame(rows)
jumps = pd.concat(jumps, ignore_index=True)

print("4/4 figures and provenance ...", flush=True)
order = ["team datapack", "new OCV set", "team OCV set"]
figures.ocv(ocv_new, team["ocv"], FIG / "ocv_new_vs_team.png")
figures.entropy(ent_candidate, team["entropy"], FIG / "entropy_candidate_vs_team.png")
for key, col, lab in [("r0", "Mean (mOhms)", "R0 (mΩ)"), ("r1", "Mean (mOhms)", "R1 (mΩ)"), ("c1", "Mean (Farad)", "C1 (F)")]:
    figures.rc({"team datapack": team[key], "new OCV set": sets["new OCV set"][key], "team OCV set": sets["team OCV set"][key]},
               col, lab, FIG / f"{key}_all_sets.png")
whole = val[(val.part == "whole file") & (val.soc_start == "SOC 1 at full charge")]
for T, tr in rpt_traces.items():
    st = whole[(whole.test == "RPT") & (whole["T"] == T)].set_index("set").rmse_mV.to_dict()
    figures.whole_file(tr, order, f"RPT {T} °C, held-out cell CS2D7485, entire file"
                       + (" (60 °C is outside the fitted range of every set)" if T == 60 else ""),
                       FIG / f"rpt_{T}C_whole_file.png", st)
figures.rpt_summary(whole[whole.test == "RPT"], order, FIG / "rpt_summary.png")
figures.hppc_overview(hppc_traces, order, {k: whole[whole.file == f].set_index("set").rmse_mV.to_dict()
                                            for k, f in hppc_file.items()}, FIG / "hppc_whole_files.png")

prov = dict(
    run_date=DATE, config=cfg,
    inputs=[dict(test=t, T_C=T, file=str(p), bytes=p.stat().st_size) for t, T, p in load.all_files()],
    team_datapack=str(load.DATAPACK),
    new_ocv=dict(gitt_points=ocv_info["gitt_points"].to_dict("records"),
                 hppc_rest_points_used=ocv_info["hppc_points_used"].to_dict("records"),
                 hppc_soc_offsets=ocv_info["hppc_alignment"].round(4).to_dict("records"),
                 branch_values_before_meeting_at_ends=ocv_info["branch_ends_before_meeting"].round(4).to_dict("records"),
                 monotonicity_violations=ocv_info["monotonicity_violations"], q_test_ah=ocv_info["q_test_ah"],
                 normal_charge_end_offset_ah=OFFSET_AH,
                 entropy_candidate=dict(table=ent_candidate.to_dict("records"), slopes=ocv_info["entropy_slopes"])),
    rc={n: dict(tau1_s=S["tau1"], hppc_soc_offsets=S["axes"].round(4).to_dict("records"),
                windows=S["fits"].assign(excluded=ecm.excluded(S["fits"])).round(7).to_dict("records"),
                per_cell_tables=S["per_cell"].round(9).to_dict("records"),
                activation_energies_eV=S["ea"].round(4).to_dict("records"), soc_end_log_fallback=S["fallback"],
                rows_measured=[15, 25, 45], rows_arrhenius=[T for T in load.R_T_GRID if T not in (15, 25, 45)],
                temperature_gate=dict(max_logged_change_K=float(np.nanmax(S["fits"][~ecm.excluded(S["fits"])].T_rise_K)),
                                      max_self_heating_bound_K=float(S["fits"][~ecm.excluded(S["fits"])].heat_bound_K.max())))
        for n, S in sets.items()},
    validation=dict(metrics=val.round(4).to_dict("records"), soc_reanchoring=jumps.round(5).to_dict("records")),
)
json.dump(prov, open(OUT / "provenance.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

# ---------------------------------------------------------------- summary for README
for n, S in sets.items():
    print(f"{n}: tau1 = {S['tau1']:.1f} s; SOC-end log fallback: {S['fallback'] or 'none'}; "
          f"gate {prov['rc'][n]['temperature_gate']}")
print("new OCV monotonicity violations:", ocv_info["monotonicity_violations"] or "none")
r = val[val.test == "RPT"]
print("\nRPT whole file, RMSE / P99 / within 50 mV:")
print(r[r.part == "whole file"].pivot_table(index=["set", "soc_start"], columns="T",
      values=["rmse_mV", "p99_mV", "within_50mV"]).round(2).to_string())
print("\nRPT by part (RMSE mV, SOC 1 start):")
print(r[(r.soc_start == "SOC 1 at full charge") & (r.part != "whole file")].pivot_table(
      index=["part", "set"], columns="T", values="rmse_mV").round(1).to_string())
h = val[(val.test == "HPPC") & (val.part == "whole file")]
print("\nHPPC whole file RMSE mV:\n", h.pivot_table(index="file", columns="set", values="rmse_mV").round(1).to_string())
print("\nSOC found at each re-anchoring:", jumps.soc_before_reanchoring.min().round(4), "-", jumps.soc_before_reanchoring.max().round(4))
print("normal-charge offset (Ah):", {k: round(v, 3) for k, v in OFFSET_AH.items()})
print("wrote", OUT)
