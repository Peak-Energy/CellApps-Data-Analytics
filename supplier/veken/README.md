# VEKEN

# metadata_builder — CyLT Master Code Generator

Jupyter notebook that ingests per-cycle CSV exports from long-term cycling (CyLT) tests, computes derived metrics, and writes a formatted multi-sheet Excel workbook.

**Current version:** `20260813_metadata_builder_v4.ipynb`

---

## Input Data

### CSV files
One CSV per cell/channel group, placed in the directory specified by `_CSV_DIR`.  
Each file must contain the following columns (produced by the MCM data conversion pipeline):

| Column | Description |
|---|---|
| `Cycle` | Cycle index (integer, 1-based) |
| `Charge_Capacity_Ah` | Charge capacity for that cycle (Ah) |
| `Discharge_Capacity_Ah` | Discharge capacity for that cycle (Ah) |
| `Charge_Energy_Wh` | Charge energy (Wh) |
| `Discharge_Energy_Wh` | Discharge energy (Wh) |
| `Energy_Efficiency_pct` | Round-trip energy efficiency expressed as a percentage (e.g. 93.5) |
| `Cell_Max_Charge_Temp_C` | Peak cell temperature during charge (°C) |
| `Cell_Max_Discharge_Temp_C` | Peak cell temperature during discharge (°C) |

### Configuration (`_CELLS` list)
Each entry in `_CELLS` describes one cell/test and carries:

| Key | Description |
|---|---|
| `csv` | Filename of the CSV (relative to `_CSV_DIR`) |
| `temperature` | Chamber temperature (°C) |
| `test_type` | Test identifier — `'CyLT'` |
| `p_rate` | Cycling rate: `0.25` or `0.5` (used to compute elapsed time) |
| `serial` | Cell serial number / label used as the sheet title |
| `file_folder` | Optional source file or folder reference |
| `type_note` | Optional test-type annotation |
| `specific_notes` | Cell-specific notes |
| `general_notes` | Batch-level notes |

---

## What the Code Does

### 1. Derived metric calculation (`_calc`)
For each CSV the notebook computes:

| Derived column | Formula |
|---|---|
| `Time_days` | `Cycle × cycle_duration_days` (see [Time approximation](#time-days-approximation) below) |
| `Coulombic_efficiency` | `Discharge_Capacity_Ah / Charge_Capacity_Ah` |
| `Energy_efficiency` | `Energy_Efficiency_pct / 100` |
| `Energy_inefficiency` | `1 − Energy_efficiency` |
| `Relative_discharge_capacity` | `Discharge_Capacity_Ah / (cycle-1 value)` |
| `Relative_discharge_energy` | `Discharge_Energy_Wh / (cycle-1 value)` |
| `Relative_energy_efficiency` | `Energy_efficiency / (cycle-1 value)` |
| `Relative_energy_inefficiency` | `Energy_inefficiency / (cycle-1 value)` |

All relative metrics are normalised to the **first cycle (index 1)**.

#### Time (days) approximation
Elapsed time is estimated from the cycle protocol — no timestamp column is required:

| P-rate | Charge | Rest (post-charge) | Discharge | Rest (post-discharge) | Total/cycle |
|---|---|---|---|---|---|
| 0.25P | 4 h | 10 min | 4 h | 30 min | **8 h 40 min ≈ 0.361 days** |
| 0.5P  | 2 h | 10 min | 2 h | 30 min | **4 h 40 min ≈ 0.194 days** |

Formula: `cycle_days = (2 / p_rate + 10/60 + 30/60) / 24`

### 2. Excel workbook generation
A single `.xlsx` file is written to `_OUTPUT_FILE`, containing:

#### `Metadata` sheet
One row per cell with the fields from `_CELLS` (cell number, serial, temperature, test type, P-rate, notes). Column A is frozen; rows alternate light-grey / white.

#### Per-cell sheets (`Cell 1`, `Cell 2`, …)
One sheet per entry in `_CELLS`. Columns written (in order):

1. Time (days) *(approximated)*
2. Cycle
3. Charge capacity (Ah)
4. Discharge capacity (Ah)
5. Relative discharge capacity
6. Coulombic efficiency
7. Charge energy (Wh)
8. Discharge energy (Wh)
9. Relative discharge energy
10. Energy efficiency
11. Energy inefficiency
12. Relative energy efficiency
13. Relative energy inefficiency
14. dV charge (V) — *placeholder, filled `NAN`*
15. dV discharge (V) — *placeholder*
16. dV charge cumulative (V) — *placeholder*
17. dV discharge cumulative (V) — *placeholder*
18. DCR charge (mΩ) — *placeholder*
19. DCR discharge (mΩ) — *placeholder*
20. DCR charge cumulative (mΩ) — *placeholder*
21. DCR discharge cumulative (mΩ) — *placeholder*
22. Ambient temperature (°C) — *placeholder*
23. Max charge temp (°C)
24. Max discharge temp (°C)
25. Thickness (mm) — *placeholder*
26. Relative thickness — *placeholder*

Columns marked *placeholder* are reserved for future data sources and output `NAN` until populated.

### 3. Styling
All sheets use a consistent Peak Energy style:
- **Dark blue** (`#1F4E79`) — section header rows
- **Medium blue** (`#2E75B6`) — column header rows
- **Alternating light-grey / white** — data rows
- Arial font throughout; frozen panes on columns A–B

---

## Output

```
_OUTPUT_FILE  (e.g. 20260813_DV2_CyLT_MasterCode_v7.xlsx)
├── Metadata         ← one row per cell
├── Cell 1           ← serial: A-DV2A184-0.25P
├── Cell 2
│   ...
└── Cell N
```

A summary is printed to the notebook output after each cell is processed, including row count and cycle-2 discharge capacity.

---

## Dependencies

```
numpy
pandas
openpyxl
pathlib (stdlib)
```


#Parsing the Veken provided spreadsheet into csv files
# Cycle Summary Export — `20260813_cycle_summary_export.ipynb`

Reads a Chinese-language Cycle Summary `.xlsx` file, translates all column headers to English, and exports user-defined `.csv` files based on Excel column ranges.

---

## Requirements

- Anaconda Python 3.x (`C:\Users\DanWindsor\anaconda3\python.exe`)
- `pandas`, `openpyxl` (included with Anaconda)
- Jupyter Notebook or JupyterLab

---

## How to use

Run all cells **top-to-bottom** after editing Cell 4. Do not edit any other cells.

### Cell layout

| Cell | Purpose |
|------|---------|
| Cell 1 | Title / instructions (Markdown) |
| Cell 2 | Imports (`os`, `re`, `pandas`) |
| Cell 3 | Configuration guide (Markdown) |
| **Cell 4** | **User configuration — edit this cell** |
| Cell 5 | Helper function header (Markdown) |
| Cell 6 | Helper functions — do not edit |
| Cell 7 | Run header (Markdown) |
| Cell 8 | Executes exports and saves CSVs |
| Cell 9 | Column layout reference header (Markdown) |
| Cell 10 | Prints block-to-column-letter map for all sheets |

---

## Configuration (Cell 4)

### Paths

```python
INPUT_FILE = (
    r"C:\Users\DanWindsor\OneDrive - Peak Energy\Documents"
    r"\MCM_DataConversion\20260813_CyLT_MCM_HT-B"
    r"\20260804_Attachment 4.  Cycle  Summary.xlsx"
)

OUTPUT_DIR = (
    r"C:\Users\DanWindsor\OneDrive - Peak Energy\Documents"
    r"\MCM_DataConversion\20260813_CyLT_MCM_HT-B\CSV_Exports"
)
```

- `INPUT_FILE` — full path to the source `.xlsx`. **The file must be closed in Excel before running** — an open file will cause a `PermissionError`.
- `OUTPUT_DIR` — folder where all `.csv` files are saved. Created automatically if it does not exist.

### Exports list

Each entry in `EXPORTS` defines one output CSV:

```python
{
    "filename":     "25C_0.25P_Cell184_A-H.csv",   # output file name
    "sheet":        "25°C 0.25P",                   # exact sheet name from the xlsx
    "col_range":    "A:H",                           # Excel column range (see below)
    "data_columns": None,                            # None = all 8 columns
},
```

| Key | Required | Description |
|-----|----------|-------------|
| `filename` | Yes | Output CSV file name |
| `sheet` | Yes | Exact sheet name from the xlsx (case-sensitive, including trailing spaces) |
| `col_range` | Yes | Excel column range covering the cell block(s), e.g. `"A:H"` |
| `data_columns` | No | List of column names to keep; `None` = all columns |

---

## Finding the right column range

Each cell occupies **8 consecutive columns** (1 cycle column + 7 measurement columns). Blocks repeat across the sheet:

| Block | Columns | Description |
|-------|---------|-------------|
| 1st cell | A:H | |
| 2nd cell | I:P | |
| 3rd cell | Q:X | |
| 4th cell | Y:AF | |
| … | … | +8 columns per cell |

To export multiple adjacent cells in one CSV, widen the range: `"A:P"` captures two cells.

Run **Cell 10** to print the full block-to-column-letter map for every sheet without opening Excel.

---

## Sheet names

| Sheet | Notes |
|-------|-------|
| `25°C 0.25P` | |
| `25°C 0.5P` | |
| `45°C 0.5P` | |
| `60°C 0.25P` | |
| `60°C 0.5P ` | Trailing space in the original sheet name |

---

## Output columns

Each exported CSV contains a `Cell_ID` column followed by the selected data columns:

| Column | Description |
|--------|-------------|
| `Cell_ID` | Battery cell identifier read from the xlsx |
| `Cycle` | Cycle number |
| `Charge_Capacity_Ah` | Charge capacity (Ah) |
| `Discharge_Capacity_Ah` | Discharge capacity (Ah) |
| `Charge_Energy_Wh` | Charge energy (Wh) |
| `Discharge_Energy_Wh` | Discharge energy (Wh) |
| `Energy_Efficiency_pct` | Energy efficiency (%) |
| `Cell_Max_Charge_Temp_C` | Cell max charge temperature (°C) |
| `Cell_Max_Discharge_Temp_C` | Cell max discharge temperature (°C) |

---

## Common errors

| Error | Cause | Fix |
|-------|-------|-----|
| `PermissionError: [Errno 13]` | Source xlsx is open in Excel | Close the file in Excel and re-run |
| `ERROR: sheet not found` | Sheet name typo or missing trailing space | Check the exact name with Cell 10 or in Excel |
| `WARNING: no data extracted` | `col_range` points outside a valid block | Run Cell 10 to verify the correct column letters |
| `WARNING: data_columns not found` | Column name typo in `data_columns` | Check spelling against the Output columns table above |
