[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://lmdutils-qfyodkmxyhworjpq24nanz.streamlit.app/)

# LMD_utils
Collection of scripts and code for working in spatial proteomics and laser micro-dissection microscopy (LMD) workflows. Useful for sample sorting before laser-dissection and for general documentation purposes.

---

## Full Workflow

### Step 1 — Re-classify GeoJSON (`import_json.py` or online app)

Prepares the QuPath GeoJSON export for the Coscia Lab converter by copying `properties.name` into `properties.classification.name` for each annotation.

**Input:** `*.geojson` exported from QuPath  
**Output:** `*_corrected.geojson` — ready for the Coscia Lab converter

**Usage (online):** Upload GeoJSON in the [online app](https://lmdutils-qfyodkmxyhworjpq24nanz.streamlit.app/), click Re-classify, download corrected file.  
**Usage (local):** Set `INPUT_FILE` in `import_json.py` and run.

---

### Step 2 — Convert GeoJSON to XML (external: Coscia Lab converter)

1. Go to [https://qupath-to-lmd-mdcberlin.streamlit.app/](https://qupath-to-lmd-mdcberlin.streamlit.app/)
2. Upload the corrected GeoJSON from Step 1
3. Check calibration points and QCs
4. Select 96-well plate (no need to update the well map — Step 3 handles sorting)
5. Confirm plate layout, process file, download zipped output

---

### Step 3 — Process LMD Collection (`process_lmd_collection.py` or online app)

Sorts samples alphabetically, assigns well positions A1–H12, updates the XML, and generates output files.

**Input** (from the Coscia Lab converter zip):
- `*.xml` — shapes with CapID + XY coordinates
- `samples_and_wells.json` — `{ "SampleName": "WellPosition" }`

**Output:**
- `*_sorted.xml` — shapes sorted by new well; `<TransferID>` added; `<CapID>` updated
- `*_96wellplate.csv` — 96-well layout, samples assigned alphabetically
- `*_sample_list.csv` — ROI number, sample name, well ID, comments, processed
- `samples_and_wells_updated.json` — updated well mapping
- `*_platemap.png` — color-coded plate map

**Usage (online):** Upload zip in the [online app](https://lmdutils-qfyodkmxyhworjpq24nanz.streamlit.app/), click Process, download outputs.  
**Usage (local):** Set `COLLECTION_FOLDER` in `process_lmd_collection.py` and run. Outputs written to `lmd_outputs/` subfolder.

---

### Step 4 — Generate MS Sample Queue (`create_ms_queue.py`)

Generates the mass spectrometry instrument queue from the sample list CSV produced in Step 3.

**Input:** `*_sample_list.csv` (Core, ROI, Well_ID, Dropout columns)

**Output:**
- `*_queue.csv` / `*_queue.xlsx` — instrument queue with K562, Supermix, and Blank controls
- `*_slot1.csv` / `*_slot1.png` — 96-well layout for samples (Slot1)
- `*_slot2.csv` / `*_slot2.png` — 96-well layout for controls (Slot2: K562→row A, Supermix→row C, Blank→row E+)

**Queue logic:**
- Start of each core: K562, Supermix, Blank
- Every group of ≤6 samples: [samples], Blank
- Samples with Dropout=Y excluded from queue, shown on plate map in grey with red label
- Spare vials (max(3, 10%) per control type) placed in Slot2 plate but not in queue

**Usage (local):** Set `INPUT_CSV` and `OUTPUT_DIR` at the top of `create_ms_queue.py`, configure LC/MS methods, run script.

---

## Scripts

| Script | Description |
|---|---|
| `import_json.py` | Re-classify GeoJSON annotations (Step 1, offline) |
| `process_lmd_collection.py` | Process LMD collection folder (Step 3, offline) |
| `process_lmd_collection_streamlit.py` | Online app for Steps 1 + 3 |
| `create_ms_queue.py` | Generate MS sample queue (Step 4, offline) |
| `sort_XML_ROI_by96well.py` | Legacy: sort XML by CapID only (superseded by Step 3) |

---

### Links
- Coscia Lab QuPath to XML converter: [https://qupath-to-lmd-mdcberlin.streamlit.app/](https://qupath-to-lmd-mdcberlin.streamlit.app/)
- Coscia Lab GitHub: [https://github.com/CosciaLab/Qupath_to_LMD](https://github.com/CosciaLab/Qupath_to_LMD)
