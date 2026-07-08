# L2-LMD — Lange Lab LMD Platform

Integrated Streamlit app for the full laser microdissection (LMD) spatial proteomics workflow,
from QuPath annotation export to mass spectrometry instrument queue.

**Run locally:**
```bash
cd L2_LMD
pip install -r requirements.txt
streamlit run app.py
```

---

## Workflow Overview

```
QuPath export (.geojson)
        |
        v
[ Tab 1: Reclassify ]
  Copy annotation name into classification.name
  (required if QuPath annotations lack classification)
        |
        v  (.geojson, reclassified)
[ Tab 2: Convert ]
  Select calibration points, assign wells, generate LMD XML
  Powered by py-lmd (MannLabs)
        |
        v  (.xml + samples_and_wells.json)  <-- or upload zip from Coscia Lab converter
[ Tab 3: Process ]
  Sort XML by well order, generate 96-well plate map + sample list CSV
        |
        v  (sample_list.csv)
[ Tab 4: MS Queue ]
  Generate instrument queue (XLSX + plate maps) with K562, Supermix, Blank controls
```

---

## Tabs

### Tab 1 — Reclassify GeoJSON
Copies `properties.name` into `properties.classification.name` for each annotation.
Required before conversion if QuPath annotations were drawn without a class assigned.
Optional if annotations already have classification set.

**Input:** `*.geojson` from QuPath export  
**Output:** `*_reclassified.geojson`

### Tab 2 — Convert to LMD XML
Converts the reclassified GeoJSON to LMD XML using the `py-lmd` library.

- Select 3 calibration points from the GeoJSON
- Wells assigned alphabetically (A1, A2, ...) with optional randomization
- QC image shows annotation coverage within calibration triangle

**Input:** `*_reclassified.geojson` (from Tab 1 or upload)  
**Output:** `*.xml`, `samples_and_wells.json`, `collection.png`

### Tab 3 — Process LMD Collection
Sorts the LMD XML by well order, renumbers shapes, assigns new well positions,
and generates all downstream files.

**Input:** XML + JSON from Tab 2 (piped directly or upload zip from Coscia Lab converter)  
> Note: zip upload is provided as fallback for users coming from the
> [Coscia Lab QuPath-to-LMD converter](https://qupath-to-lmd-mdcberlin.streamlit.app/)

**Output:** `*_sorted.xml`, `*_96wellplate.csv`, `*_sample_list.csv`, `*_platemap.png`

### Tab 4 — MS Sample Queue
Generates the Bruker timsTOF instrument queue from the sample list CSV.

- K562, Supermix, Blank controls with configurable loads
- Spare vials (max(3, 10%) per control type)
- One-slot mode: fit samples + controls into a single 96-well plate when possible
- Dropout handling: excluded from queue, shown in grey with red label on plate map

**Input:** `*_sample_list.csv` (from Tab 3 or upload)  
**Output:** `*_queue.xlsx`, `*_slot1.csv/png`, `*_slot2.csv/png`

---

## File Structure

```
L2_LMD/
  app.py                  # main Streamlit app (4 tabs)
  requirements.txt
  README.md
  utils/
    __init__.py
    geojson_utils.py      # Tab 1: reclassify logic
    convert_utils.py      # Tab 2: py-lmd Collection wrapper, calibration, well assignment
    plate_utils.py        # shared: well layout helpers, plate visualization
    process_utils.py      # Tab 3: sort XML, generate plate CSV + sample list
    ms_queue_utils.py     # Tab 4: queue builder, slot2 layout, controls
```

---

## Attribution and Dependencies

**Core XML generation:**  
[py-lmd](https://github.com/MannLabs/py-lmd) — Wallmann, Madler, Schmacke et al., MannLabs / Hornung Lab.
Apache-2.0 license.
> Schmacke et al. (2023) SPARCS, a platform for genome-scale CRISPR screening for spatial cellular phenotypes. bioRxiv. https://doi.org/10.1101/2023.06.01.542416

**Inspiration and workflow reference:**  
[Qupath_to_LMD](https://github.com/CosciaLab/Qupath_to_LMD) — Coscia Lab, MDC Berlin.

---

## Links

- Coscia Lab online converter (alternative for Tab 2): https://qupath-to-lmd-mdcberlin.streamlit.app/
- py-lmd documentation: https://mannlabs.github.io/py-lmd/
- py-lmd GitHub: https://github.com/MannLabs/py-lmd
