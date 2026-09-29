# EVE

# EVE A3-155Ah — CaLT ModelingInput Builder

**Notebook:** `metadata_builder_EVE-A3-155Ah_CaLT.ipynb`

Converts raw EVE A3-155Ah calendar aging (CaLT) source data into a formatted Excel workbook that matches the 31-column reference schema used by Peak Energy's calendar fade model pipeline.

---

## What the notebook does

1. **Reads** the translated EVE source Excel file, one sheet per storage condition.
2. **Parses** per-cell, per-period measurement rows, extracting retained and recovered capacity/energy values.
3. **Computes** all derived metrics: relative retained/recovered capacity, self-discharge, relative energy, etc.
4. **Prepends** a synthetic Day 0 baseline row for each cell (all relative metrics = 1.0, self-discharge = 0 mAh).
5. **Applies** an optional 99.18% correction factor to `Rel_Ret_Cap` and `Rel_Dis_Cap` for the three storage conditions specified in `_CAP_CORR_SHEETS`.
6. **Writes** a formatted Excel workbook containing a `Metadata` sheet followed by one sheet per cell.

---

## Input data format

**File:** `NF155L storage data (26.5.8)_translated.xlsx`

The source file contains one worksheet per storage condition (e.g. `25℃-100%SOC Storage`). Each sheet follows the layout below. All column indices are **0-based**.

### Sheet layout

| 0-based col | Label in source | Content | Notes |
|---|---|---|---|
| 1 | Period | Storage period label (`1M`, `2M`, … or `7d*N`) | Spans multiple rows per group; forward-filled by the parser |
| 2 | No. | Cell serial number | String; header rows contain `"No."` or `"Data"` and are skipped |
| 6 | ① Capacity | Before-test (initial) discharge capacity | **Labeled "mAh" in source — actual unit is Ah** (see note below) |
| 7 | ① Energy | Before-test discharge energy | Wh — correctly labeled |
| 11 | ③ Capacity | Residual/retained discharge capacity | **Labeled "mAh" in source — actual unit is Ah**; populated at 100% SOC only; `/` or blank at 0%/50% SOC |
| 12 | ③ Energy | Residual/retained discharge energy | Wh — correctly labeled |
| 13 | ④ Capacity | Recovered discharge capacity | **Labeled "mAh" in source — actual unit is Ah** |
| 14 | ④ Energy | Recovered discharge energy | Wh — correctly labeled |

Columns not listed above are not read by the parser.

### Unit mislabeling — capacity columns

Columns 6, 11, and 13 are labeled **"mAh"** in the source file header, but all values are in **Ah** (e.g., 160.4 for a nominally 155 Ah cell). This is confirmed by the energy columns: 449 Wh ÷ 160 Ah ≈ 2.81 V average discharge voltage — physically consistent.

**The parser reads these values directly in Ah and applies no unit conversion.** If the source file is ever corrected to true mAh values, divide by 1000 at columns 6, 11, and 13 in `_parse_sheet`.

### Period label format

The parser converts period labels to days using the following rules:

| Label pattern | Example | Converted to |
|---|---|---|
| `7d*N` | `7d*4` | 7 × N days = 28 days |
| `NM` | `3M` | N × 30.44 days ≈ 91.3 days |

### What is NOT in the source

The following data required by the reference schema is absent from the EVE source file and is set to NaN throughout:

- Recovered **charge** capacity and energy (only recovered discharge is available)
- Coulombic efficiency, energy efficiency, energy inefficiency (require charge data)
- Relative energy efficiency / inefficiency
- Time-weighted and capacity-weighted ΔV
- Temperature columns (ambient, max charge, max discharge)
- DCIR (10 s, 50% SOC; SOC-weighted)
- Thickness

Self-discharge is approximated using `shift(1) of recovered discharge capacity − retained discharge` (in Ah, converted to mAh for output). The reference formula uses a measured SOC-adjustment charge that is not present in the source.

Retained capacity (col 11) is only populated at 100% SOC conditions; it is NaN at 0% and 50% SOC.

---

## Configuration

All user-editable settings are in the **CONFIG cell** (`b2c3d4e5`):

```python
_INPUT_FILE   # Absolute path to the source Excel file
_OUTPUT_FILE  # Absolute path for the output workbook

_SHEET_MAP    # List of (source_sheet_name, temperature_°C, average_soc) tuples
              # Order determines Cell ID numbering in the output

_CAP_CORR         # Correction factor scalar (currently 0.9918)
_CAP_CORR_SHEETS  # Set of sheet names that receive the correction factor
_CAP_CORR_NOTE    # Text written to the Specific Notes column for corrected sheets

_ANNOTATIONS  # Dict keyed by serial number for per-cell overrides:
              #   {'type_note': '', 'specific_notes': '', 'general_notes': ''}
```

---

## Output format

**File:** `EVE-A3-155Ah_CaLT_MasterCode_v{N}.xlsx`

### Metadata sheet

One row per cell. Columns:

| Column | Content |
|---|---|
| Cell ID | Sequential integer (order follows `_SHEET_MAP`) |
| Temperature | Storage temperature (°C) |
| Test Type | `CaLT` |
| Calendar: Average SOC | Fractional SOC (e.g. `1.0`, `0.5`, `0.0`) |
| Cycling: P-rate | `0` (not applicable for calendar aging) |
| Serial number | Cell serial from source file |
| Type note | From `_ANNOTATIONS` |
| Specific notes | Correction factor note (for corrected sheets) or from `_ANNOTATIONS` |
| General notes | From `_ANNOTATIONS` |

### Per-cell sheets

One sheet per cell, labelled `Cell N`. Each sheet contains the 31 reference columns:

| # | Column | Internal name | Unit |
|---|---|---|---|
| 1 | Time (days) | `Time_days` | days |
| 2 | Cycle | `Cycle` | — (NaN) |
| 3 | Retained discharge capacity | `Retained_Cap_Ah` | Ah |
| 4 | Recovered charge capacity | `Rec_Chg_Cap_Ah` | Ah (NaN) |
| 5 | Recovered discharge capacity | `Rec_Dis_Cap_Ah` | Ah |
| 6 | Relative retained discharge capacity | `Rel_Ret_Cap` | — |
| 7 | Relative discharge capacity | `Rel_Dis_Cap` | — |
| 8 | Self discharge | `Self_Dis_Ah` | mAh |
| 9 | Coulombic efficiency | `CE` | — (NaN) |
| 10 | Retained discharge energy | `Retained_E_Wh` | Wh |
| 11 | Recovered charge energy | `Rec_Chg_E_Wh` | Wh (NaN) |
| 12 | Recovered discharge energy | `Rec_Dis_E_Wh` | Wh |
| 13 | Relative retained discharge energy | `Rel_Ret_E` | — |
| 14 | Relative discharge energy | `Rel_Dis_E` | — |
| 15–18 | EE, EI, Rel_EE, Rel_EI | — | NaN |
| 19–22 | ΔV columns | — | mV (NaN) |
| 23–25 | Temperature columns | — | °C (NaN) |
| 26–29 | DCIR columns | — | mΩ (NaN) |
| 30–31 | Thickness, Relative thickness | — | mm (NaN) |

Row 1 of each cell sheet is a synthetic Day 0 baseline (Time = 0, all relative metrics = 1.0, Self_Dis_Ah = 0.0). Subsequent rows correspond to measurement periods in source order.

---

## 99.18% correction factor

Three storage conditions require a correction factor on relative capacity metrics, per the source spreadsheet methodology:

- `25℃-100%SOC Storage`
- `45℃-100%SOC Storage`
- `25℃-0%SOC Storage`

For these sheets, `Rel_Ret_Cap` and `Rel_Dis_Cap` are multiplied by `0.9918` at every measured time point. The Day 0 baseline row always remains at `1.0` regardless of the correction. A note is written to the `Specific notes` column of the Metadata sheet for each affected cell.

---

## Known artefact — Rel_Dis_Cap > 1 at 0% SOC

Some cells stored at 0% SOC show `Rel_Dis_Cap` slightly above 1.0 at early checkup periods. This is a measurement artefact: the before-test reference capacity (col 6) was measured at 0.5P, while checkup recovered capacity (col 13) was measured at 0.25P. The lower C-rate yields marginally higher measured capacity, causing the ratio to exceed unity for cells with minimal degradation. No physical capacity gain is implied.

---

## Dependencies

```
numpy
pandas
openpyxl
```

All available in the base Anaconda environment. Run with:

```
C:\Users\DanWindsor\anaconda3\python.exe -m jupyter notebook
```

