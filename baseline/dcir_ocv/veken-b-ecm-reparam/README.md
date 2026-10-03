# Veken 170Ah-HT -B: OCV, 1RC and 2RC parameters from GITT and HPPC

This folder rebuilds the -B cell's OCV and RC tables from the raw lab files and compares them with the current team
datapack (`cell-performance/cells/Veken170Ah-HT-B`). Built 2026-10-01/02 by DJ with Claude. Every decision below
was made with DJ unless it says *my choice*.

Which OCV to use is the team's call, and the 1RC and 2RC are both on the table, so there are **four complete
parameter sets**. Each is a copy of the team's cell folder with the new files in place, ready to drop into
`cell-performance/cells/`:

| Folder (`outputs/2026-10-02/…`) | OCV | Model | Files that differ from the team folder |
|---|---|---|---|
| `new_ocv/Veken170Ah-HT-B/` | **new OCV** built here from the GITT (below) | 1RC, τ1 = 178 s | `ocv.csv`, `r0.csv`, `r1.csv`, `c1.csv`, `cell.yaml` |
| `team_ocv/Veken170Ah-HT-B/` | **team `ocv.csv`**, unchanged | 1RC, τ1 = 210 s | `r0.csv`, `r1.csv`, `c1.csv`, `cell.yaml` |
| `new_ocv_2rc/Veken170Ah-HT-B/` | new OCV | 2RC: slow branch τ1 = 1579 s in `r1`/`c1`, fast branch τ2 = 34 s in `r2`/`c2` | `ocv.csv`, `r0.csv`, `r1.csv`, `c1.csv`, `r2.csv`, `c2.csv`, `cell.yaml` |
| `team_ocv_2rc/Veken170Ah-HT-B/` | team `ocv.csv`, unchanged | 2RC: slow branch τ1 = 1664 s in `r1`/`c1`, fast branch τ2 = 39 s in `r2`/`c2` | `r0.csv`, `r1.csv`, `c1.csv`, `r2.csv`, `c2.csv`, `cell.yaml` |

- **Other files.** All other files are byte-for-byte copies of the team folder, including `entropy.csv` in every set.
- **`cell.yaml`.** It is the team file with only two kinds of change: the status comment, and the `model_parameters`
  entries for the files that changed (description, origin, retrieved date). The 2RC sets add `r2` and `c2` entries
  shaped like `r1` and `c1`.
- **`r2.csv` / `c2.csv`** have the same grid, columns and units as `r1.csv` / `c1.csv`. The team's simulator reads
  only `r0`/`r1`/`c1` today; it would have to be extended to use the second branch. In the 2RC sets `r1`/`c1` is the
  slow branch, so a simulator that ignores `r2`/`c2` still gets the minutes-scale response (which branch is "1" is a
  naming choice, open to the team).
- **R1, C1, R2, C2 depend on the OCV.** The fit subtracts the OCV from the measured voltage, so each set's RC tables
  belong with its own `ocv.csv` and shouldn't be mixed across sets. R0 doesn't depend on the OCV or the model and is
  the same in all four sets.

**Results on the held-out RPT cell (whole files, RMSE at 15 / 25 / 45 / 60 °C):**

| Set | RMSE (mV) |
|---|---|
| team datapack | 27 / 25 / 23 / 21 |
| new-OCV set (1RC) | 24 / 25 / 28 / 30 |
| team-OCV set (1RC) | 21 / 23 / 23 / 21 |
| new-OCV 2RC set | 23 / 25 / 32 / 49 |
| team-OCV 2RC set | 29 / 36 / 28 / 44 |

**The 2RC does not beat the 1RC on the held-out cell.** It is marginally better on the 10–150 s pulses and worse
on the hour-long charges and discharges, where its slow branch (about 1600 s, fitted on the 2 h rests) keeps
building polarization that the cell doesn't show. Adding a fast branch to the 1RC's own time constant changes the
whole-file numbers by at most 0.3 mV. Details under "2RC".

---

## Where things are

```
README.md        this document
run.py           python3 run.py  rebuilds everything below (parses raw files into ~/.cache on first run)
config.yaml      the values decided with DJ (threshold, capacities, temperature cap, OCV options, fit bounds)
reparam/
  load.py        raw-file locations, Neware CSV loader, datapack file reader/writer
  ocv.py         OCV and dU/dT from GITT (+ HPPC rests)
  ecm.py         RC models, HPPC windows, 1RC and 2RC fits, gridding, cold-temperature fill
  validate.py    whole-file simulation of RPT and HPPC records
  figures.py     the plots
tests/           python3 -m pytest tests   (unit tests + synthetic 1RC and 2RC HPPC files with known answers)
outputs/2026-10-02/
  new_ocv/Veken170Ah-HT-B/         drop-in cell folder, new OCV, 1RC
  team_ocv/Veken170Ah-HT-B/        drop-in cell folder, team OCV, 1RC
  new_ocv_2rc/Veken170Ah-HT-B/     drop-in cell folder, new OCV, 2RC
  team_ocv_2rc/Veken170Ah-HT-B/    drop-in cell folder, team OCV, 2RC
  provenance.json                  every rest point, every fit window, every validation number, for all sets
  figures/                     listed under "Figures"
```

Raw data (read-only): `~/Documents/Projects/Battery Modeling Data/Veken/Veken-B-POR/` (GITT, HPPC, RPT).

---

## Which OCV each set uses

- **New-OCV set.** GITT rest voltages of the same cells the team used, read as the average of every sample where the
  voltage drifts less than 1 mV/h.
  - **SOC axis.** SOC = 1 − ΔQ / Q(T), counted from the rested state after a slow (C/20) charge to 3.45 V.
  - **SOC 0.** The rested voltage after a normal-rate discharge to 1.5 V.
  - **Gaps between GITT points.** Filled with the HPPC 2 h rest voltages of the same cells at 15 and 45 °C. At
    25 °C they're filled with the 15/45 °C curve shape, bent through the 25 °C GITT points.
  - **Charge curve.** The discharge curve plus the measured charge/discharge gap.
  - Details are under "OCV" in "Every decision, in plain words".
- **Team-OCV set.** The team's `ocv.csv` exactly as it is in `cell-performance`. Where it comes from is described in
  "Where the team's OCV comes from".

---

## Results

### Held-out cell (CS2D7485, RPT files), entire file simulated

Every row of each RPT file is simulated: capacity cycles, the pulse-and-move block and the final cycles.
Error = model − measured voltage. 60 °C is outside the fitted range of every set.

| Parameter set | RMSE (mV), 15 / 25 / 45 / 60 °C | 99th percentile \|error\| (mV) | Share within 50 mV |
|---|---|---|---|
| Team datapack | 26.9 / 25.1 / 22.6 / 21.2 | 148 / 121 / 122 / 92 | 0.97 / 0.97 / 0.95 / 0.97 |
| New-OCV set (1RC) | 24.4 / 24.7 / 28.2 / 30.3 | 132 / 90 / 123 / 101 | 0.97 / 0.93 / 0.93 / 0.87 |
| New-OCV set, starting at the SOC a normal-rate charge reaches (see "SOC = 100 %") | 21.6 / 19.4 / 30.5 / 43.7 | 110 / 82 / 102 / 154 | 0.97 / 0.97 / 0.87 / 0.83 |
| Team-OCV set (1RC) | 20.5 / 23.4 / 22.8 / 21.4 | 98 / 79 / 94 / 87 | 0.98 / 0.95 / 0.94 / 0.97 |
| New-OCV 2RC set | 23.4 / 24.7 / 32.1 / 48.6 | 77 / 90 / 115 / 214 | 0.95 / 0.93 / 0.90 / 0.86 |
| Team-OCV 2RC set | 28.7 / 36.3 / 27.9 / 44.1 | 145 / 170 / 114 / 235 | 0.96 / 0.94 / 0.92 / 0.89 |

Mean RMSE over 15–45 °C: team datapack 24.9 mV, new-OCV set 25.8 mV, team-OCV set 22.2 mV, new-OCV 2RC set
26.7 mV, team-OCV 2RC set 31.0 mV.

Where the errors are (RMSE mV, 15 / 25 / 45 / 60 °C):

| Part of the test | Team datapack | New-OCV set | Team-OCV set | New-OCV 2RC set | Team-OCV 2RC set |
|---|---|---|---|---|---|
| 10–60 s pulses | 50 / 53 / 22 / 28 | 19 / 28 / 32 / 33 | 42 / 47 / 24 / 28 | 21 / 27 / 33 / 34 | 41 / 45 / 24 / 28 |
| 150 s pulses | 60 / 52 / 35 / 55 | 33 / 34 / 44 / 69 | 66 / 54 / 35 / 56 | 29 / 30 / 44 / 69 | 64 / 53 / 34 / 54 |
| Constant-current discharges (hours) | 23 / 16 / 24 / 19 | 27 / 20 / 20 / 35 | 16 / 9 / 20 / 20 | 24 / 19 / 36 / 63 | 24 / 32 / 29 / 48 |
| Constant-current charges (hours) | 16 / 15 / 14 / 15 | 16 / 25 / 32 / 21 | 13 / 23 / 22 / 15 | 19 / 25 / 22 / 21 | 23 / 28 / 22 / 34 |
| Constant-power discharges (minutes) | 68 / 59 / 61 / 36 | 53 / 38 / 25 / 47 | 36 / 23 / 40 / 31 | 47 / 44 / 58 / 86 | 41 / 52 / 46 / 56 |
| Constant-power charges (minutes) | 92 / 79 / 35 / 50 | 31 / 39 / 57 / 44 | 73 / 71 / 44 / 49 | 45 / 43 / 46 / 44 | 92 / 91 / 51 / 67 |
| SOC 10–90 % | 8 / 5 / 8 / 12 | 10 / 12 / 15 / 23 | 8 / 10 / 12 / 13 | 17 / 15 / 21 / 29 | 15 / 16 / 16 / 18 |
| SOC below 10 % | 73 / 67 / 62 / 53 | 64 / 57 / 71 / 64 | 54 / 59 / 58 / 52 | 50 / 51 / 75 / 120 | 72 / 92 / 69 / 121 |
| SOC above 90 % | 12 / 9 / 6 / 8 | 8 / 27 / 8 / 15 | 9 / 10 / 9 / 9 | 9 / 27 / 9 / 13 | 9 / 11 / 5 / 6 |
| Rests | 65 / 44 / 29 / 58 | 63 / 53 / 44 / 42 | 65 / 43 / 28 / 60 | 41 / 44 / 52 / 71 | 87 / 84 / 58 / 100 |

What this says:
- **Team-OCV set.** Its OCV is identical to the team datapack's, so every difference from the datapack comes from
  the new R0/R1/C1. It is lower than the datapack at 15 and 25 °C (27 → 21, 25 → 23 mV) and within 0.3 mV at
  45/60 °C. The biggest gain is on the constant-power discharges, the minutes-long load.
- **New-OCV set.** It has the lowest 1RC error on the pulses and constant-power charges at 15/25 °C. Its error is
  higher between 10 and 90 % SOC, above 90 % SOC at 25 °C, and over whole files at 45/60 °C. See "How the two OCVs
  differ".
- **2RC sets.** Against the 1RC set on the same OCV, the 2RC is 1–4 mV better on the 10–150 s pulses and on the
  constant-power charges with the new OCV, and 4–30 mV worse on the hour-long constant-current charges and
  discharges, below 10 % SOC and at 60 °C. The whole-file figures show where: the slow branch (τ ≈ 1600 s) keeps
  building polarization through the hour-long steps and overshoots at the cut-off voltages. Both 2RC sets drop the
  share of rows within 50 mV.
- **Coulomb counting drifts slowly.** Every time the cell is fully charged again, the counted SOC reads
  1.001–1.004 instead of 1. Between full charges the cell takes in 0.1–0.4 % of capacity more than it gives out
  (side reactions / self-discharge). provenance.json lists each re-anchoring.

### HPPC files, entire file

These files are in-sample for all four new sets. Their RC tables are fitted on them, and the new OCV also uses
their 2 h rests at 15/45 °C. RMSE in mV:

| File | Team datapack | New-OCV set | Team-OCV set | New-OCV 2RC set | Team-OCV 2RC set |
|---|---|---|---|---|---|
| 15 °C CS2D0052 / CS2D9049 | 18.5 / 14.5 | 20.9 / 16.1 | 18.3 / 14.2 | 20.8 / 16.1 | 18.2 / 14.3 |
| 25 °C CS2D7013 / CS2D9805 | 19.0 / 19.9 | 26.5 / 27.7 | 18.2 / 19.0 | 26.7 / 27.8 | 18.8 / 19.5 |
| 45 °C CS2D2520 / CS2D8533 | 39.9 / 46.0 | 45.5 / 50.7 | 39.9 / 46.0 | 45.5 / 50.7 | 39.9 / 45.9 |

- **1RC and 2RC are within 0.6 mV of each other on every HPPC file.** The HPPC-file error is dominated by SOC and
  OCV, not by the RC branches that were fitted on these very files.
- **45 °C.** All sets drift upward late in the test, at low SOC (`hppc_whole_files.png`). The same drift in
  every set points to SOC or OCV at 45 °C rather than to any RC model.
- **25 °C.** The 100 h rest at full charge is the other large error for all sets: the model holds the OCV while the
  cell keeps relaxing.

### Where the team's OCV comes from

Checked against the raw files and the `cell-performance` git history:

- **The test files are separate tests.**
  - The RPT files are cell CS2D7485 only, and that cell appears in no GITT or HPPC file.
  - The GITT and HPPC files are cells CS2D0052, 9049, 7013, 9805, 2520 and 8533.
  - Dates (2026): 25 °C HPPC Feb 20 – Mar 2, RPT Feb 27 – Mar 30, 25 °C GITT Mar 12–16, 15/45 °C GITT Jul 17–30,
    15/45 °C HPPC Aug 4–11.
- **The team's `ocv.csv` history:**
  - Jun 23: a 1 % table at 15/25/45 °C.
  - Jul 6: "switch -B ocv curves to use RPT-derived curves", 13 points per temperature at 5–10 % spacing.
  - Aug 5: the current 1 % table.
- **The interior (1–99 % SOC) is GITT data.**
  - It matches the GITT voltages at the scheduled end of each rest, on the same cells, within 0.3–0.8 mV RMS on both
    branches at 15, 25 and 45 °C.
  - That's measured at 10–90 % SOC, after a SOC shift of 0.1–0.9 % that differs by branch.
- **The SOC 0 and SOC 1 rows are from the RPT test cell.**
  - At 15 and 25 °C they're carried unchanged from the Jul 6 RPT-derived version: SOC 0 = 1.858 / 1.843 V and
    SOC 1 = 3.358 / 3.361 V.
  - These match the RPT test cell's 3 h rested voltages after its normal-rate discharge (1.852 / 1.837 V) and
    charge (3.359 / 3.361 V).
  - At 45 °C, SOC 1 (3.377 V) is also from that version. SOC 0 (1.624 V) matches the 45 °C GITT bottom after the
    slow discharge (1.626 V).
- **The "RPT" source in `cell.yaml`** applies to those SOC 0 / SOC 1 rows only.
- **What this means for validation.** Between 1 and 99 % SOC the RPT test cell is independent of both OCVs. At
  SOC 0 and SOC 1 (15 / 25 °C) the team's values are that cell's own measurements.

### How the two OCVs differ

| | Team | New | Effect |
|---|---|---|---|
| Reading each GITT rest | voltage at the scheduled end of the rest | average of every sample where the voltage drifts less than 1 mV/h | at 20–80 % SOC the two readings differ by 0.2–0.5 mV. Charge/discharge gap at 50 % SOC: 17 / 13 / 19 mV (new reading), 16 / 13 / 18 mV (end of rest), 13 / 15 / 13 mV (team) at 15 / 25 / 45 °C |
| SOC axis | interior sits 0.1–0.9 % SOC from the new axis, by a different amount for each branch (exact rule not identified) | SOC = 1 − ΔQ / Q(T) from the slow-charge 100 %, the same for both branches | about 0.9 % SOC between branches is worth about 5 mV at 50 % SOC: most of the gap difference at 15/45 °C |
| SOC 1 | RPT test cell, rested after a normal-rate charge | GITT cells, rested after a slow (C/20) charge | new is +22 / +58 / +15 mV at 15 / 25 / 45 °C |
| SOC 0 | RPT test cell, rested after a normal-rate discharge (15/25 °C); GITT slow-discharge bottom (45 °C) | GITT cells, rested after a normal-rate discharge | 1.95 / 1.91 / 1.75 V new vs 1.86 / 1.84 / 1.62 V team |
| Filling GITT gaps | none | HPPC 2 h rests (15/45 °C); 15/45 °C curve shape (25 °C) | at 45 °C the HPPC cells' rested voltage at 25 % SOC sits 9–10 mV above the team curve and within 1.3 mV of the new one |

### New-OCV choices that were tested

Each choice was changed one at a time and scored with a 1RC refitted on that OCV, over the whole RPT files. RMSE mV
at 15 / 25 / 45 / 60 °C; "mean" is over 15–45 °C.

| Variant | RMSE | Mean | In the new-OCV set? |
|---|---|---|---|
| 25 °C gaps filled with HPPC rests, branches averaged at SOC 0 | 28 / 39 / 27 / 35 | 31.5 | no |
| … with the 25 °C curve below 20 % SOC from the 15/45 °C shape | 28 / 26 / 27 / 31 | 27.2 | — |
| … with SOC 0 = rested voltage after a normal-rate discharge to 1.5 V | 24 / 32 / 28 / 33 | 28.3 | — |
| both of the above | 24 / 25 / 28 / 30 | 25.8 | **yes** |
| both, plus 100 % = rested after a normal-rate charge (as the team) | 26 / 20 / 30 / 47 | 25.6 | **no (DJ)**: 100 % stays at the slow charge, the closest state to equilibrium; a normal-rate 100 % would put a slowly charged cell above 100 % |
| a narrower charge/discharge gap | — | — | not tested (DJ): the measured GITT gap is kept |

### Figures (`outputs/2026-10-02/figures/`)

| File | Shows |
|---|---|
| `rpt_summary.png` | whole-file RMSE and 99th-percentile error on the RPT cell, all five sets, by temperature |
| `rpt_15C_whole_file.png` … `rpt_60C_whole_file.png` | measured vs all five sets over the whole RPT file, error with ±50 mV band, current |
| `hppc_whole_files.png` | model − measured over each whole HPPC file, all five sets |
| `ocv_new_vs_team.png` | charge, discharge, average OCV at 15/25/45 °C, new vs team, with the difference |
| `entropy_candidate_vs_team.png` | the team's dU/dT and the candidate from the new OCV (in neither set) |
| `r0_all_sets.png`, `r1_all_sets.png`, `c1_all_sets.png` | each table at all 8 temperatures for all five sets; grey panels are extrapolated rows |
| `r2_all_sets.png`, `c2_all_sets.png` | the fast branch of the two 2RC sets |

---

## How the R0 / R1 / C1 tables differ from the team datapack, and why

Percentages are medians over 20–80 % SOC (new-OCV set / team-OCV set, both 1RC). The 2RC sets have the same R0; their
`r1` (slow branch, τ1 = 1579 / 1664 s) is 31–36 % below the team R1 at 15/25 °C and 7–13 % below at 45 °C, and their
`r2` (fast branch, τ2 = 34 / 39 s) has no team counterpart.

| Quantity | Team datapack | Both new sets | Difference | Main cause |
|---|---|---|---|---|
| R0 | voltage change on the first logged row of each pulse. That row is captured mid-switch and shows about 30 % of the real drop | voltage change at the first complete sample, 0.1 s into the pulse (same in all sets) | 0.34 vs 0.080 mΩ at 25 °C; 0.57 vs 0.067 at 15 °C; 0.19 vs 0.074 at 45 °C; about 8× in the cold rows | method; same 25 °C raw data |
| τ1 | fixed near 362 s at every SOC and temperature | one τ1 fitted to all SOC moves at all temperatures: 178 s (new-OCV set), 210 s (team-OCV set) | — | method; τ1 depends on the OCV subtracted in the fit |
| R1 | 10 s pulse resistance − R0 | fitted per SOC and temperature with R0 and τ1 held | 25 °C −26 / −24 %; 15 °C −15 / −13 %; 45 °C −10 / −3 %; −30 °C −27 / −38 % | method: the team R1 is the 10 s response, the new R1 the minutes-long one |
| C1 | from the team τ1 and R1 | τ1 / R1 | 25 °C −34 / −24 %; 15 °C −42 / −34 %; 45 °C −46 / −41 % | follows τ1 and R1 |
| 25 °C row | not measured: it equals the Arrhenius interpolation of the team's 15 and 45 °C rows | measured | see above | convention |
| Lowest SOC rows | team rows 0 / 0.025 / 0.05 / 0.075 hold the HPPC points measured at 2.5 / 5 / 7.5 / 10 % | true SOC; SOC 0 extrapolated | low-SOC values shifted by one HPPC step | convention |
| Cold rows (−30 to 10 °C) | R0 follows the -A cell's temperature trend | Arrhenius through the measured 15/25/45 °C rows, separately at each SOC | R0 ≈8× higher, R1 lower (above) | method |
| Entropy | approximated from RPT OCV | not changed: all sets keep the team's `entropy.csv` | — | the candidate from the new OCV is in provenance.json only (see "Entropy") |

---

## Every decision, in plain words

### Reading the data
- **Time inside a second.** The files print time to the whole second but log pulses every 0.1 s, so ten rows
  share each printed second. They're spaced evenly inside it (DJ). Using the charge counter instead would add up to
  ±21 ms of rounding at 8.5 A.
- **Rows the logger didn't capture.** The first row of every step is logged mid-switch (it shows the full current
  but only about 30 % of the voltage step). When the second row falls within about 0.5 s and the step is logged
  coarser than 0.1 s, that row repeats the first row's voltage exactly: it happens on every SOC move and every rest
  after one, in the HPPC and the RPT files alike, and never on the 0.1 s pulse logging. Both rows are left out of
  every fit and every validation statistic (*my choice*).
- **Test pauses.** Three files stopped logging for a while: 28 h in both 15 °C GITT files, 23 h in the 45 °C GITT,
  5 min in both 45 °C HPPC files. The missing time is added back from the wall clock so relaxation is timed
  correctly (*my choice*).
- **Cell temperature.**
  - The 15 °C files have two cell probes (T1a, T2a) and one ambient probe that is written to two columns
    (T1b = T2b). Cell temperature = average of the two cell probes (DJ).
  - The 15 °C cells actually sat at 15.2–19.3 °C. Fits use the logged value, but the table rows still say 15 °C.
  - The 25/45 °C and RPT files have no temperature channel, so they use the chamber setpoint.
- **Temperature gate.**
  - A rest or fit window is dropped if the cell temperature changes by more than 2 K within it (DJ's rule).
  - Because the probes jitter ±0.3 K from sample to sample, the change is measured between 10-minute averages
    (*my choice*, approved).
  - This dropped 3 of the 172 GITT rests, all on the 15 °C charge curve during a day-long chamber warm-up. No HPPC
    window came near the limit (largest change 0.9 K).
- **Self-heating where temperature isn't logged.** Upper bound ∫I²R dt / (m·c_p) with m = 4.69 kg and
  c_p = 917 J/kg·K (DJ's rule). The largest bound in any window used is 0.56 K.

### OCV (new-OCV set)
- **Equilibrium voltage of a rest.**
  - The average of every sample where the voltage changes by less than 1 mV/h (DJ).
  - The rate comes from a straight-line fit over the preceding 173 s (1 s logging) or 1200 s (300 s logging).
    Those lengths keep the noise error at or below 0.2 mV/h, given the measured voltage noise of 37 and 48 µV
    (*my choice*, following the spec's rule).
  - The result is within 0.2–0.5 mV of the voltage at the scheduled end of the rest, which is what the team uses.
- **Rests that never get there.** 16 of 172 GITT rests never relax to 1 mV/h; 12 of them are the 2 h rests at
  25 °C. They use their final voltage and are flagged in provenance.json (DJ).
- **SOC.** SOC = 1 − (charge removed) / Q(T), with Q(T) = 171.69 / 172.4 / 177.73 / 179.08 Ah at 15 / 25 / 45 /
  60 °C, the values in cell.yaml (DJ: keep the team's capacities).
- **SOC = 100 %.** The rested state after a slow (C/20) charge to 3.45 V: the closest thing to equilibrium in these
  data (DJ).
  - A normal-rate (≈0.24C) charge to 3.45 V stops 0.76 / 0.95 / 1.42 Ah short of 100 % at 15 / 25 / 45 °C, i.e. at
    SOC ≈ 0.996 / 0.994 / 0.992.
  - Starting the RPT simulation there instead of at 1 gives 21 / 19 / 31 / 44 mV (Results table).
- **SOC = 0 %.** At SOC 0 both curves equal the rested voltage after a normal-rate discharge to 1.5 V, from each
  GITT test's own conditioning discharge (DJ). The slow discharge's points past that state aren't used, so the
  discharge curve rises steadily from there.
- **Filling the GITT gaps.** The GITT steps are 10 % apart in mid-SOC, and at 25 °C everywhere.
  - **15 and 45 °C:** the discharge curve adds the voltages at the end of the HPPC test's 2 h rests on the same
    cells, every 2.5–5 % SOC (DJ).
    - Each HPPC file's SOC count is shifted by one constant (0.3–0.5 % SOC, either direction) so its rests sit on
      the same cell's GITT curve. This cuts the mismatch from 3–6 mV to 1–3 mV (*my choice*).
    - An HPPC point is used only where no GITT point lies within 1 % SOC (*my choice*).
  - **25 °C:** the 15/45 °C discharge curve shape, bent through the 25 °C GITT points (DJ; tested above).
- **Charge curve.** Discharge curve + the charge/discharge gap. At 15 and 45 °C the gap is measured by the GITT.
  At 25 °C the GITT is too coarse, so the gap is taken from 15 and 45 °C and corrected to the 25 °C GITT values
  between 10 and 90 % SOC (DJ).
- **Smooth curve.** Points are joined with a shape-preserving interpolation (it can't overshoot between points)
  onto 1 % SOC steps (spec). All curves rise steadily with SOC.
- **Top of the table.** At SOC 1 both branches are set to their average, as in the team file (DJ).
- **Two cells per temperature.** Averaged (DJ). The 45 °C GITT has one cell.
- **No coulombic-efficiency correction.** The slow discharge pulls 1.9–3.6 Ah more than the charge put in, because
  it ends at 1.5 V at a much lower current than the conditioning discharge. That's a different end state, not
  lost charge (*my choice*).

### Entropy (candidate, in neither set)
- **What it is.** The slope of the new OCV against temperature across the 15/25/45 °C cells at each SOC, on the same
  SOC axis as the new `ocv.csv` (DJ: "keep it consistent"), using the logged 15 °C temperature.
- **Values.** −0.50 mV/K at 50 % SOC (team −0.29), and −1.3 to −7 mV/K below 20 %. Values below 10 % and above
  90 % SOC are low-confidence.
- **Why it isn't the entropic coefficient.** Q(T) grows 3.6 % from 15 to 45 °C, so this slope includes capacity
  growth; it's not the fixed-charge entropic coefficient. On a fixed-charge axis mid-SOC would be about −0.1 to
  −0.3 mV/K.
- **The team's own files disagree.** The team's `ocv.csv` implies the same steep slopes as this candidate, which
  doesn't match the team's `entropy.csv`.
- **Where to find it.** In provenance.json and `entropy_candidate_vs_team.png`. All sets ship the team's
  `entropy.csv`.

### 1RC (both 1RC sets, same procedure)
- **R0.** Voltage change divided by current at the first complete sample (0.1 s) of the discharge pulses, median
  over the pulse currents (DJ). The same R0 is used by the 2RC sets.
  - R0 doesn't depend on pulse current (0.340–0.341 mΩ from 56 to 170 A at 25 °C).
- **Which response the RC branch represents.** The cell shows responses on three time scales, about 5–10 s, about
  50 s and about 1600 s (see "2RC"); a 1RC holds one compromise.
  - The minutes-scale one was chosen, with one time constant for every SOC and temperature (DJ).
  - τ1 comes from all SOC-move windows fitted together; R1 is fitted per SOC and temperature with R0 held.
  - Alternatives tested on the held-out cell (new-OCV set):
    - 10 s pulses only (τ ≈ 6–15 s) under-predicted sustained load;
    - one τ1 per temperature (221 / 130 / 1255 s) gave cold rows that run the wrong way with temperature;
    - a free τ1 per window scattered from 8 to 3400 s.
- **OCV inside the fits.** Each set's own discharge curve, subtracted sample by sample along the counted SOC.
  - The HPPC history is all discharge moves, and 10 s pulses barely move the hysteresis (*my choice*).
  - Each HPPC file's SOC count is shifted by one constant so its 2 h rests sit on that OCV (*my choice*). The shift
    is −0.4 to +0.7 % SOC with the new OCV and −0.4 to +0.5 % with the team OCV.
- **Fit details** (*my choice*; each was needed to recover known values from the synthetic test file):
  - each window's voltage zero level is fitted rather than pinned to one noisy sample;
  - the RC state entering a window is simulated from the last long rest;
  - samples are weighted by their time spacing;
  - current is held at each step's mean.
- **Windows left out:**
  - **All sets** (*my choice*):
    - pulses at the 100 % point and the first SOC move. They start fresh off the charge, mid-way from the charge
      to the discharge OCV branch.
    - pulses at the last (2.5 %) point. That point runs charge-first and gave R1 of 7–38 mΩ at 45 °C.
  - **Team-OCV set only, the 45 °C move from 25 to 20 % SOC on both cells (DJ).**
    - At 25 % the cells' 2 h rested voltage sits 9–10 mV above the team OCV. Every other rest from 12 to 40 % is
      within about 4 mV, except the lowest one at 12–13 % (7–8.5 mV).
    - The fit absorbs that mismatch by driving R1 to 0 (the non-negative fit's boundary), against 0.34–0.41 mΩ for
      the neighbouring moves.
    - With the move left out, the 45 °C value at 20 % is interpolated from its neighbours (0.37 mΩ) and the −30 °C
      value is 2.9 mΩ. The zero rule under "Tables" gives the same table, so the leave-out changes no number.
  - **Windows used.** The new-OCV set uses 138 SOC-move windows and the team-OCV set 136. Both use 414 pulse
    windows.
- **Residual shape.** Per move window the single branch leaves 2.3 / 2.2 mV RMSE and +1.0 / +0.8 mV mean residual
  during the move (new-OCV / team-OCV set). The rests after the moves are logged every 300 s, so nothing faster
  than that is visible in them.

### 2RC (both 2RC sets, same procedure)
- **Model.** R0 as in the 1RC, plus a fast and a slow RC branch, each with one time constant for every SOC and
  temperature; R_fast and R_slow vary with SOC and temperature (DJ).
- **Joint fit (DJ).** Every HPPC window takes part at once: the discharge pulses (10 s + 60 s rest, 0.1 s
  logging), the charge pulses (10 s + 900 s rest, 0.1 s logging) and the SOC moves (5–10 min + 2 h rest, logged
  every 30–300 s). The two time constants are chosen on a log grid of 80 values from 0.5 to 5000 s, by the summed
  weighted error over all windows after the best non-negative (R_fast, R_slow) pair has been solved at each SOC
  point, then refined between grid points. The 1RC's fit details (free zero level per window, RC history from the
  last long rest, time-spacing weights, step-mean current, excluded windows, SOC shift) are kept unchanged.
- **Per SOC point.** R_fast and R_slow are solved together from the point's pulses and move, with R0 held. R_fast
  is kept from points that have usable pulses and R_slow from points that have a usable move; a point missing one
  kind (the first and the last point) holds the branch it cannot see at the nearest complete point's value.
- **What came out.**
  - New-OCV set: τ_fast = 34 s, τ_slow = 1579 s. Team-OCV set: τ_fast = 39 s, τ_slow = 1664 s. Neither sits at
    a grid edge, and the error surface has a single minimum.
  - The windows disagree about the time constants. The pulse windows alone prefer about 5 s + 47 s; the move
    windows alone prefer 37–47 s + 1560–1750 s. The joint fit lands between them because the moves carry more
    time weight. A 2RC holds two of the cell's three time scales.
  - On the fit windows, the pair leaves 0.42 / 0.40 mV RMS against 0.96 / 0.91 mV for the best single time
    constant (135 / 151 s). Per window: pulses 0.6–0.7 mV RMSE, moves 1.4 mV.
  - R_fast is small and noisy at 45 °C (0.03–0.2 mΩ, alternating between neighbouring points), and two points at
    45 °C / 10 % SOC (new OCV) and one at 95 % (team OCV, alternative only) fitted to exactly 0 and are
    interpolated. Its activation energies span −0.02 to 0.8 eV, so the cold rows of `r2` are not to be trusted;
    the slow branch's activation energies also go negative at some SOC (−0.3 to 0.4 eV).
- **Team-OCV 2RC set only, two 45 °C moves left out (ending at 25 % and at 20 % SOC).** With two branches the
  OCV mismatch at 25 % spoils both moves around it: the move ending at 20 % drives R_fast to 0, the one ending at
  25 % drives R_slow to 0.04 mΩ against 0.1–0.4 mΩ at the neighbours (and through the cold rows to 63 mΩ at
  −30 °C). Both points are interpolated from 30 and 15 %. The team-OCV 2RC set therefore uses 134 move windows;
  both 2RC sets use 828 pulse windows (discharge and charge).
- **Tested alternatives** (validation only, no folder; whole-file RPT RMSE at 15 / 25 / 45 / 60 °C):

  | Variant | New OCV | Team OCV |
  |---|---|---|
  | joint fit (the 2RC sets) | 23.4 / 24.7 / 32.1 / 48.6 | 28.7 / 36.3 / 27.9 / 44.1 |
  | τ_fast = 8 s (what the pulses prefer at the joint τ_slow), joint τ_slow | 24.6 / 27.5 / 32.3 / 48.4 | 29.8 / 38.9 / 27.8 / 43.5 |
  | τ_fast = 8 s added on top of the 1RC's τ1 (178 / 210 s) | 24.8 / 24.4 / 28.5 / 30.5 | 20.6 / 23.2 / 22.9 / 21.4 |
  | the 1RC set itself | 24.4 / 24.7 / 28.2 / 30.3 | 20.5 / 23.4 / 22.8 / 21.4 |

  Adding a fast branch to the 1RC changes the whole-file error by at most 0.3 mV: the 10 s response is a small
  part of the whole-file error, and the gain the 2RC sets show on the pulses is paid for several times over on the
  hour-long steps by the 1600 s branch.
- **Naming.** `r1`/`c1` hold the slow branch and `r2`/`c2` the fast one (`config.yaml`, `two_rc.branch_1`), so a
  simulator that reads only `r1`/`c1` still gets the minutes-scale response. Open to the team.

### Tables
- **Grid.** The team's 24 SOC points and 8 temperatures. The SOC 0 and 1 rows are extrapolated in a straight line
  from the two nearest measured points (DJ).
- **Mean / Error.** Mean = median of the two cells, Error = standard deviation of the two cells (spec; sample
  standard deviation is *my choice*). C1 = τ1 / R1 (and C2 = τ2 / R2) after combining the cells.
- **A resistance that fitted to exactly 0** is the non-negative fit's boundary, not a measurement: that point is
  interpolated from its neighbours and listed in provenance.json (*my choice*).
- **Cold rows.**
  - At each SOC, ln R0 and ln R1 vs 1/T are fitted through 15 (at the logged ≈16.7 °C), 25 and 45 °C and evaluated
    at −30 to 10 °C (spec).
  - Activation energies (1RC sets): R0 0.27–0.31 eV and R1 0.14–0.46 eV (new-OCV set); R0 0.28–0.32 eV and R1
    0.17–0.62 eV (team-OCV set). The 2RC sets' are under "2RC".
  - Error on these rows = value × the median relative Error of the measured rows at that SOC (*my choice*).

### Validation
- Entire RPT and HPPC files are simulated row by row with the same simulator for all sets, with one RC branch or
  two as the set's files say.
- SOC is counted from the end of each rested full charge, where it is set to 1. Rows before the first full charge
  are counted back from it. SOC is never inferred from voltage (*my choice* for whole files, following the spec's
  rule for the start).
- The OCV branch switches to the most recent current direction the moment the current changes (*my choice*).
- 60 °C for all sets: OCV extended linearly in temperature from the 25 and 45 °C rows; R0, R1, τ1 (and R2, τ2) by
  Arrhenius through the 15/25/45 °C rows (*my choice*).
- The rows the logger didn't capture (see "Reading the data") are left out of the error statistics.

---

## Checks

`python3 -m pytest tests` runs 28 tests.

- **Synthetic 1RC file.** An HPPC-shaped Neware file with:
  - known R0, R1 and τ1 = 360 s;
  - realistic noise and logging;
  - mid-switch first rows;
  - a known SOC shift.

  The loader, SOC shift and fitter recover the known values. Without noise they come back to within 3 %. With
  noise, R0 is within 1 % and R1/τ1 are within the fit's own error bars.
- **Synthetic 2RC file.** The same file with a slow branch (0.25 mΩ, 360 s) and a fast one (0.12 mΩ, 8 s). The
  joint fit recovers both time constants and both resistances within 3 % without noise and 5 % with noise, keeps
  one pair of resistances per SOC point, and beats the best single time constant on the same windows by more than
  a factor of 2.
- **Unit tests:**
  - sign flip;
  - sub-second time;
  - duplicated temperature headers;
  - RC step response, and the two-branch simulation as the sum of two step responses;
  - the unusable-row rule (mid-switch first row, stale second row);
  - relaxation detection on a known exponential;
  - interpolation without overshoot;
  - per-sample OCV subtraction;
  - Arrhenius recovery;
  - round-trip of all five team files through the reader and writer.
- **Drop-in folders.** Each run checks:
  - each `cell.yaml` loads as YAML;
  - only the listed files differ from the team folder.
