"""Rebuild everything:  python3 run.py

1. New OCV from GITT (+ HPPC rests).
2. From HPPC, with the new OCV and with the team's current OCV subtracted: a 1RC (R0, R1/C1) and a 2RC
   (the same R0; a slow branch in r1/c1 and a fast one in r2/c2).
3. Four drop-in cell folders (same files and cell.yaml format as cell-performance/cells/Veken170Ah-HT-B):
     outputs/<run_date>/new_ocv/Veken170Ah-HT-B/       new ocv.csv + r0/r1/c1 fitted on it (1RC)
     outputs/<run_date>/team_ocv/Veken170Ah-HT-B/      team ocv.csv (unchanged) + r0/r1/c1 fitted on it (1RC)
     outputs/<run_date>/new_ocv_2rc/Veken170Ah-HT-B/   new ocv.csv + r0/r1/c1/r2/c2 fitted on it (2RC)
     outputs/<run_date>/team_ocv_2rc/Veken170Ah-HT-B/  team ocv.csv (unchanged) + r0/r1/c1/r2/c2 fitted on it (2RC)
   Every other file is copied unchanged from the team datapack.
4. Whole-file simulation of every RPT and HPPC record with the team datapack and all four sets; figures;
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
BRANCHES = ecm.BRANCHES_2RC[cfg["two_rc"]["branch_1"]]
FAST_IS = "2" if cfg["two_rc"]["branch_1"] == "slow" else "1"

# ---------------------------------------------------------------- 1-2. OCV and the four fits
print("1/4 new OCV ...", flush=True)
ocv_new, ent_candidate, ocv_info = ocv.build(cfg)
print("2/4 1RC and 2RC fits ...", flush=True)
sets, variants = {}, {}
for name, ocv_tab, leave_out, leave_out_2rc in [("new OCV", ocv_new, [], []),
                                                ("team OCV", team["ocv"], cfg["hppc"]["team_ocv_leave_out"],
                                                 cfg["two_rc"]["team_ocv_leave_out"])]:
    prep = ecm.prepare(cfg, ocv_tab, leave_out)
    fits, tau1 = ecm.fit_1rc(prep)
    full, per_cell, ea, fallback = ecm.tables(fits)
    sets[f"{name} set"] = dict(ocv=ocv_tab, fits=fits, tau1=tau1, axes=prep[2], per_cell=per_cell, ea=ea,
                               fallback=fallback, leave_out=leave_out, **ecm.to_datapack(full))
    if leave_out_2rc != leave_out:
        prep = ecm.prepare(cfg, ocv_tab, leave_out_2rc)
    fits2, taus2, diag2 = ecm.fit_2rc(prep, cfg["two_rc"])
    full2, per_cell2, ea2, fallback2 = ecm.tables(fits2, BRANCHES)
    sets[f"{name} 2RC set"] = dict(ocv=ocv_tab, fits=fits2, taus=taus2, diag=diag2, axes=prep[2], per_cell=per_cell2,
                                   ea=ea2, fallback=fallback2, leave_out=leave_out_2rc, **ecm.to_datapack(full2))
    # tested alternatives (validation only, no folder): the fast time constant the pulse windows alone prefer,
    # with the joint slow one; and the same fast one added on top of the 1RC's own time constant
    tf_p = diag2["by_window_kind"]["pulse"]["best_fast_s_at_chosen_slow"]
    for label, pair in [(f"{name} 2RC, fast tau from pulses", (tf_p, taus2["slow"])),
                        (f"{name} 2RC, fast tau from pulses + 1RC tau1", (tf_p, tau1))]:
        fits3, taus3, diag3 = ecm.fit_2rc(prep, cfg["two_rc"], fixed=pair)
        full3, per_cell3, ea3, fallback3 = ecm.tables(fits3, BRANCHES)
        variants[label] = dict(ocv=ocv_tab, fits=fits3, taus=taus3, diag=diag3, axes=prep[2], per_cell=per_cell3,
                               ea=ea3, fallback=fallback3, leave_out=leave_out_2rc, **ecm.to_datapack(full3))
ORDER = ["team datapack", "new OCV set", "team OCV set", "new OCV 2RC set", "team OCV 2RC set"]
FOLDER = {"new OCV set": "new_ocv", "team OCV set": "team_ocv", "new OCV 2RC set": "new_ocv_2rc",
          "team OCV 2RC set": "team_ocv_2rc"}

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
ROWS = ("15/25/45 C measured (Mean = median of 2 cells, Error = std of 2 cells); -30 to 10 C per-SOC Arrhenius through "
        "the measured rows.")


def block(key, description, origin):
    lines = [f"  {key}:", f"    file: {key}.csv", "    description: >"]
    words, line = description.split(), "     "
    for w in words:
        if len(line) + len(w) + 1 > 100:
            lines.append(line); line = "     "
        line += " " + w
    lines += [line, "    source:", f"      origin: {origin}", f'      retrieved: "{DATE}"']
    return "\n".join(lines) + "\n"


def _block_pattern(key):
    return re.compile(rf"  {key}:\n    file: {key}\.csv\n.*?retrieved: \"[^\"]*\"\n", re.S)


def replace_block(text, key, new):
    pat = _block_pattern(key)
    if not pat.search(text):
        raise ValueError(f"{key} block not found in cell.yaml")
    return pat.sub(lambda m: new, text, count=1)


def insert_after_block(text, key, new):
    m = _block_pattern(key).search(text)
    if not m:
        raise ValueError(f"{key} block not found in cell.yaml")
    return text[:m.end()] + "\n" + new + text[m.end():]


def left_out_text(leave_out):
    if not leave_out:
        return ""
    socs = " and ".join(f"{100 * x['soc']:.0f}" for x in leave_out)
    return (f" The 45 C move{'s' if len(leave_out) > 1 else ''} ending at {socs} % SOC {'are' if len(leave_out) > 1 else 'is'} "
            f"left out (the cells' rested voltage at 25 % sits 9-10 mV above this OCV); {socs} % "
            f"{'are' if len(leave_out) > 1 else 'is'} interpolated from the neighbours.")


def rc_blocks_1rc(text, S, on):
    text = replace_block(text, "r0", block("r0",
        "Ohmic resistance: voltage change at the first complete sample (0.1 s) of the HPPC discharge pulses, median over "
        f"pulse currents. {ROWS} Fitted with {on} subtracted.", ORIGIN_HPPC))
    text = replace_block(text, "r1", block("r1",
        f"1RC resistance of the minutes-scale response, from HPPC SOC moves + 2 h rests with one time constant for all "
        f"SOC and temperatures (tau1 = {S['tau1']:.0f} s), R0 held. Rows as for r0. Fitted with {on} subtracted."
        + left_out_text(S["leave_out"]), ORIGIN_HPPC))
    text = replace_block(text, "c1", block("c1", f"C1 = tau1 / R1 with tau1 = {S['tau1']:.0f} s. Fitted with {on} subtracted.",
                                           ORIGIN_HPPC))
    return text


def rc_blocks_2rc(text, S, on):
    tf, ts = S["taus"]["fast"], S["taus"]["slow"]
    slow, fast = ("r1", "r2") if FAST_IS == "2" else ("r2", "r1")
    tau = {"r1": ts if FAST_IS == "2" else tf, "r2": tf if FAST_IS == "2" else ts}
    text = replace_block(text, "r0", block("r0",
        "Ohmic resistance: voltage change at the first complete sample (0.1 s) of the HPPC discharge pulses, median over "
        f"pulse currents. {ROWS} Fitted with {on} subtracted.", ORIGIN_HPPC))
    desc = {slow: f"2RC model, slow branch: resistance of the minutes-scale response, time constant {ts:.0f} s for all SOC "
                  f"and temperatures. Fitted jointly with the fast branch ({fast}/{fast.replace('r', 'c')}) on the HPPC "
                  f"pulses, their rests and the SOC moves + 2 h rests, R0 held. Rows as for r0. Fitted with {on} subtracted."
                  + left_out_text(S["leave_out"]),
            fast: f"2RC model, fast branch: resistance of the seconds-scale response, time constant {tf:.1f} s for all SOC "
                  f"and temperatures. Fitted jointly with the slow branch ({slow}/{slow.replace('r', 'c')}); see {slow}. "
                  f"Rows as for r0. Fitted with {on} subtracted."}
    text = replace_block(text, "r1", block("r1", desc["r1"], ORIGIN_HPPC))
    text = replace_block(text, "c1", block("c1", f"C1 = tau1 / R1 with tau1 = {tau['r1']:.4g} s (2RC model, see r1).", ORIGIN_HPPC))
    text = insert_after_block(text, "c1", block("r2", desc["r2"], ORIGIN_HPPC) + "\n"
                              + block("c2", f"C2 = tau2 / R2 with tau2 = {tau['r2']:.4g} s (2RC model, see r2).", ORIGIN_HPPC))
    return text


yaml_src = (load.DATAPACK / "cell.yaml").read_text()
status_old = ("# -B DATAPACK STATUS (2026-09-09): nominal values are MEASURED -B data\n"
              "# (mass / specific heat / capacity); OCV (RPT, 2026-08-04) and R0 / R1 / C1\n"
              "# (CAPE-215 HPPC, 2026-05-20) are -B files. ")
assert status_old in yaml_src, "cell.yaml status comment changed upstream"
for set_name, folder in FOLDER.items():
    S = sets[set_name]
    two = "taus" in S
    new_ocv = set_name.startswith("new")
    dest = OUT / folder / "Veken170Ah-HT-B"
    shutil.copytree(load.DATAPACK, dest, ignore=shutil.ignore_patterns("Claude outputs", ".DS_Store"))
    keys = ["r0", "r1", "c1"] + (["r2", "c2"] if two else [])
    for key in keys:
        load.write_table(S[key], dest / f"{key}.csv", key[0])
    on = "the new OCV (ocv.csv in this folder)" if new_ocv else "the team ocv.csv (unchanged, in this folder)"
    text = (rc_blocks_2rc if two else rc_blocks_1rc)(yaml_src, S, on)
    rc_what = "R0 / R1 / C1 / R2 / C2" if two else "R0 / R1 / C1"
    how = f"(HPPC {'2RC ' if two else ''}reparam fitted on that OCV, {DATE})"
    if new_ocv:
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
        ocv_what = f"OCV (GITT reparam, {DATE})"
    else:
        ocv_what = "OCV (unchanged, 2026-08-04)"
    status = (f"# -B DATAPACK STATUS ({DATE}): nominal values are MEASURED -B data\n"
              f"# (mass / specific heat / capacity); {ocv_what} and {rc_what}\n"
              f"# {how} are -B files. ")
    if two:
        status += ("\n# r2 / c2 are the second RC branch of a 2RC model; a simulator that reads only r0 / r1 / c1\n"
                   "# ignores them. ")
    text = text.replace(status_old, status)
    yaml.safe_load(text)                                       # still valid YAML
    (dest / "cell.yaml").write_text(text)
    changed = {f.name for f in dest.iterdir()
               if not (load.DATAPACK / f.name).exists() or not filecmp.cmp(f, load.DATAPACK / f.name, shallow=False)}
    expected = {f"{k}.csv" for k in keys} | {"cell.yaml"} | ({"ocv.csv"} if new_ocv else set())
    assert changed == expected, f"{folder}: unexpected changed files {changed ^ expected}"
for a, b in [("new_ocv", "new_ocv_2rc"), ("team_ocv", "team_ocv_2rc")]:
    assert filecmp.cmp(OUT / a / "Veken170Ah-HT-B" / "r0.csv", OUT / b / "Veken170Ah-HT-B" / "r0.csv", shallow=False), \
        f"r0.csv differs between {a} and {b}"

# ---------------------------------------------------------------- 4. validation, figures, provenance
print("3/4 whole-file simulations ...", flush=True)
team_rc = ecm.from_datapack(team["r0"], team["r1"], team["c1"])
sim = {"team datapack": (team["ocv"], team_rc)}
sim.update({n: (S["ocv"], ecm.from_datapack(*[S[k] for k in ["r0", "r1", "c1", "r2", "c2"] if k in S]))
            for n, S in (sets | variants).items()})
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
figures.ocv(ocv_new, team["ocv"], FIG / "ocv_new_vs_team.png")
figures.entropy(ent_candidate, team["entropy"], FIG / "entropy_candidate_vs_team.png")
for key, col, lab in [("r0", "Mean (mOhms)", "R0 (mΩ)"), ("r1", "Mean (mOhms)", "R1 (mΩ)"), ("c1", "Mean (Farad)", "C1 (F)"),
                      ("r2", "Mean (mOhms)", "R2 (mΩ)"), ("c2", "Mean (Farad)", "C2 (F)")]:
    tabs = {n: (team[key] if n == "team datapack" else sets[n][key]) for n in ORDER if n == "team datapack" and key in team or n in sets and key in sets[n]}
    figures.rc(tabs, col, lab, FIG / f"{key}_all_sets.png")
whole = val[(val.part == "whole file") & (val.soc_start == "SOC 1 at full charge")]
for T, tr in rpt_traces.items():
    st = whole[(whole.test == "RPT") & (whole["T"] == T)].set_index("set").rmse_mV.to_dict()
    figures.whole_file(tr, ORDER, f"RPT {T} °C, held-out cell CS2D7485, entire file"
                       + (" (60 °C is outside the fitted range of every set)" if T == 60 else ""),
                       FIG / f"rpt_{T}C_whole_file.png", st)
figures.rpt_summary(whole[whole.test == "RPT"], ORDER, FIG / "rpt_summary.png")
figures.hppc_overview(hppc_traces, ORDER, {k: whole[whole.file == f].set_index("set").rmse_mV.to_dict()
                                            for k, f in hppc_file.items()}, FIG / "hppc_whole_files.png")


def rc_prov(S):
    f = S["fits"]
    used = f[~ecm.excluded(f)]
    out = dict(hppc_soc_offsets=S["axes"].round(4).to_dict("records"), leave_out=S["leave_out"],
               windows=f.assign(excluded=ecm.excluded(f)).round(7).to_dict("records"),
               per_cell_tables=S["per_cell"].round(9).to_dict("records"),
               activation_energies_eV=S["ea"].round(4).to_dict("records"), soc_end_log_fallback=S["fallback"],
               rows_measured=[15, 25, 45], rows_arrhenius=[T for T in load.R_T_GRID if T not in (15, 25, 45)],
               temperature_gate=dict(max_logged_change_K=float(np.nanmax(used.T_rise_K)),
                                     max_self_heating_bound_K=float(used.heat_bound_K.max())))
    if "taus" in S:
        out.update(model="2RC", tau_fast_s=S["taus"]["fast"], tau_slow_s=S["taus"]["slow"],
                   branches={"r1/c1": cfg["two_rc"]["branch_1"], "r2/c2": "fast" if FAST_IS == "2" else "slow"},
                   joint_fit=S["diag"])
    else:
        out.update(model="1RC", tau1_s=S["tau1"])
    return out


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
    rc={n: rc_prov(S) for n, S in sets.items()},
    rc_tested_alternatives={n: rc_prov(S) for n, S in variants.items()},
    validation=dict(metrics=val.round(4).to_dict("records"), soc_reanchoring=jumps.round(5).to_dict("records")),
)
json.dump(prov, open(OUT / "provenance.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

# ---------------------------------------------------------------- summary for README
for n, S in (sets | variants).items():
    gate = (prov["rc"] | prov["rc_tested_alternatives"])[n]["temperature_gate"]
    if "taus" in S:
        d = S["diag"]
        print(f"{n}: tau_fast = {S['taus']['fast']:.1f} s, tau_slow = {S['taus']['slow']:.0f} s; joint-fit RMS "
              f"{d['rms_mV_best_pair']:.3f} mV; points: {d['n_points_complete']} complete, {d['n_points_held']} with one "
              f"branch held; notes: {S['fallback'] or 'none'}; gate {gate}")
        if "grid_best" in d:
            print(f"   {d['used']} from grid best {d['grid_best']}, grid edge: {d['at_grid_edge']}; best single tau "
                  f"{d['best_single_tau_s']:.0f} s gives {d['rms_mV_best_single_tau']:.3f} mV; by window kind: {d['by_window_kind']}")
        f = S["fits"]; u = f[~ecm.excluded(f)]
        for kind in ["pulse", "move"]:
            g = u[u.kind == kind]
            print(f"   {kind} windows: n {len(g)}, RMSE {g.rmse_mV.mean():.2f} mV, mean residual during the step "
                  f"{g.resid_during_mV.mean():+.2f} mV, first 10 s of rest {np.nanmean(g.resid_first_10s_of_rest_mV):+.2f} mV")
    else:
        print(f"{n}: tau1 = {S['tau1']:.1f} s; notes: {S['fallback'] or 'none'}; gate {gate}")
        f = S["fits"]; g = f[(~ecm.excluded(f)) & (f.kind == "move")]
        print(f"   move windows: n {len(g)}, RMSE {g.rmse_mV.mean():.2f} mV, mean residual during the move "
              f"{g.resid_during_move_mV.mean():+.2f} mV, first 10 s of rest {np.nanmean(g.resid_first_10s_of_rest_mV):+.2f} mV")
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
