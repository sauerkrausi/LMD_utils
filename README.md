| [![L2-LMD](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://lmdutils-jvcvkwqujnfbvutyhfs5l.streamlit.app/) | [![L2-msQ](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://l2msqueue.streamlit.app/) |
|:---:|:---:|
| L2-LMD | L2-msQ |

# LMD_utils
Scripts and Streamlit apps for spatial proteomics laser micro-dissection (LMD) workflows.

The primary tool is the **L2-LMD integrated app** (`L2_LMD/`), which covers the full pipeline from QuPath export to MS instrument queue in four tabs. A standalone MS queue generator (`L2_msQ/`) is also available for cases where only the queue step is needed.

---

## L2-LMD — Integrated App (`L2_LMD/`)

**Run locally:**
```bash
cd L2_LMD
pip install -r requirements.txt
streamlit run L2_LMD_app.py
```

### Tab 1 — Reclassify GeoJSON
Copies `properties.name` into `properties.classification.name` for each QuPath annotation.
Required if annotations were drawn without a class assigned. Can be skipped via toggle if annotations are already classified.

**Input:** `*.geojson` from QuPath export
**Output:** `*_reclassified.geojson`

---

### Tab 2 — Convert to LMD XML
Converts reclassified GeoJSON to LMD XML using [py-lmd](https://github.com/MannLabs/py-lmd) (MannLabs, Apache-2.0).

- Select 3 calibration points; coverage check confirms ROI centroids lie within the calibration triangle
- Define sample groups / supergroups from name prefixes; add custom prefixes to subdivide
- Well assignment across 96-well or 384-well plates (user-selectable)
  - 384-well: configurable margin, space between rows, and space between columns for easier pipetting
  - Group-aware bin-packing keeps sample groups together while balancing plates
- Plate map PNG preview and cutting list CSV available before XML conversion

**Input:** `*_reclassified.geojson`
**Output:** `*_PlateN.xml` per plate, `samples_and_wells.json`, plate map PNGs, cutting list CSV

---

### Tab 3 — Post-Cutting QC
Mark dropout ROIs after stereomicroscope inspection. Dropout ROIs are excluded from the MS queue and from re-run XMLs.

**Input:** Sample list piped from Tab 2 (or CSV upload)
**Output:** `*_sample_list.csv` (with Dropout Y/N), `*_rerun.zip` (XMLs with dropouts removed)

---

### Tab 4 — MS Sample Queue
Generates the Bruker timsTOF instrument queue. Run groups re-indexed from 1 after dropout removal.

- K562, Supermix, Blank controls with configurable loads and spare vials (max(3, 10%) per type)
- User-defined blank interval (samples per block) and optional randomized run order within groups
- Group size summary, divide-group helper, and per-plate slot assignment
- Dropout wells shown in grey with red label on plate map

**Input:** `*_sample_list.csv` (from Tab 3 or upload)
**Output:** `*_queue.xlsx`, per-slot CSVs and plate map PNGs (zip)

---

## L2 MS Queue — Standalone App (`L2_msQ/`)

Lightweight standalone queue generator. No GeoJSON or XML required.

**Input:** CSV with `sample_name` and `well` columns (optional `group` column)
**Features:** Auto-group detection from filename, group editor, divide-group helper, randomize toggle + block size
**Output:** `*_queue.xlsx`, plate map PNGs (zip)

```bash
cd L2_msQ
pip install -r requirements.txt
streamlit run app.py
```

---

## Legacy Scripts

These standalone scripts are superseded by the L2-LMD app but retained for offline/scripted use.

| Script | Description |
|---|---|
| `import_json.py` | Re-classify GeoJSON annotations (Tab 1, offline) |
| `process_lmd_collection.py` | Process LMD collection folder (offline) |
| `process_lmd_collection_streamlit.py` | Older combined Streamlit app (Steps 1+3) |
| `create_ms_queue_streamlit.py` | Older standalone MS queue Streamlit app |
| `sort_XML_ROI_by96well.py` | Legacy: sort XML by CapID only |

---

## Links
- py-lmd (MannLabs): [https://github.com/MannLabs/py-lmd](https://github.com/MannLabs/py-lmd) — Apache-2.0
- Coscia Lab converter: [https://qupath-to-lmd-mdcberlin.streamlit.app/](https://qupath-to-lmd-mdcberlin.streamlit.app/)
- Coscia Lab GitHub: [https://github.com/CosciaLab/Qupath_to_LMD](https://github.com/CosciaLab/Qupath_to_LMD)