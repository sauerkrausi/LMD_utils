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
  Select calibration points
  Define sample groups / supergroups
  Assign wells across 96-well plates (group-aware bin-packing)
  Generate LMD XML(s), plate maps, cutting list
  Powered by py-lmd (MannLabs)
        |
        v  (.xml per plate + samples_and_wells.json)
[ Tab 3: Post-Cutting QC ]
  Mark dropout ROIs after stereomicroscope inspection
  Generate dropout-free XMLs for LMD re-runs on fresh sections
        |
        v  (sample_list.csv, dropouts excluded)
[ Tab 4: MS Queue ]
  Generate instrument queue (XLSX + plate maps) with K562, Supermix, Blank controls
  Run groups re-indexed from 1 after dropout removal
```

---

## Tabs

### Tab 1: Reclassify GeoJSON
Copies `properties.name` into `properties.classification.name` for each annotation.
Required before conversion if QuPath annotations were drawn without a class assigned.

**Input:** `*.geojson` from QuPath export  
**Output:** `*_reclassified.geojson`

---

### Tab 2: Convert to LMD XML
Converts the reclassified GeoJSON to LMD XML using the `py-lmd` library.

**Calibration:** Select 3 point annotations as calibration markers. Coverage check confirms
all ROI centroids lie within the calibration triangle.

**Sample Groups:** Auto-detected from name prefixes (first hyphen-delimited token, e.g. `AS20`, `CH22`, `H20`).
Groups control which ROIs are kept together on the same plate and flow through to the MS queue.
- Edit the **Group** dropdown per supergroup to reassign
- Add custom prefixes (e.g. `H20-015885`) to subdivide a supergroup — the longer prefix wins
  in matching, and the original supergroup disappears automatically once all its ROIs are captured
- Supergroups with 0 ROIs (fully captured by a longer prefix) are hidden automatically
- Per-ROI fine-tuning available in the expandable ROI-level table

**Well Assignment:** Alphabetical with optional randomization. Multi-plate support with
group-aware bin-packing: sample groups stay together while plates are balanced as evenly as possible.

**Input:** `*_reclassified.geojson` (from Tab 1 or upload)  
**Output (immediate):** Plate map PNGs, cutting list CSV  
**Output (after Convert):** `*_PlateN.xml` per plate, `samples_and_wells.json`, all-plates zip

---

### Tab 3: Post-Cutting QC
Mark ROIs that failed cutting after inspecting the plate under the stereomicroscope.
Dropout ROIs are excluded from the MS queue (Tab 4) and from the re-run XMLs.

**Input:** Sample list piped from Tab 2 (or CSV upload)  
**Output:**
- `*_sample_list.csv` — updated with Dropout Y/N, piped to Tab 4
- `*_rerun.zip` — XMLs with dropout ROIs removed, for re-running LMD on a fresh section

---

### Tab 4: MS Sample Queue
Generates the Bruker timsTOF instrument queue from the sample list CSV.
Dropout ROIs are automatically excluded. Run groups are re-indexed from 1.

- K562, Supermix, Blank controls with configurable loads
- Spare vials (max(3, 10%) per control type)
- Group-aware: each group gets controls at the start, then a Blank after every 6 samples
- Per-plate slot assignment (Plate1: Slot1+Slot2, Plate2: Slot3+Slot4, ...)
- Dropout wells shown in grey with red label on plate map

**Input:** `*_sample_list.csv` (from Tab 3 or upload)  
**Output:** `*_queue.xlsx`, per-slot CSVs and plate map PNGs (zip)

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
    convert_utils.py      # Tab 2: py-lmd wrapper, calibration, group editor, well assignment
    process_utils.py      # Tab 3: dropout editor, dropout-free XML generation
    ms_queue_utils.py     # Tab 4: queue builder, slot layout, controls
```

---

## Attribution and Dependencies

**Core XML generation:**  
[py-lmd](https://github.com/MannLabs/py-lmd) — Wallmann, Madler, Schmacke et al., MannLabs / Hornung Lab.
Apache-2.0 license.
> Schmacke et al. (2023) SPARCS, a platform for genome-scale CRISPR screening for spatial cellular phenotypes. bioRxiv. https://doi.org/10.1101/2023.06.01.542416

**Inspiration and workflow reference:**  
[Qupath_to_LMD](https://github.com/CosciaLab/Qupath_to_LMD) — Coscia Lab, MDC Berlin.
GPL-3.0 license.

---

## Links

- Coscia Lab online converter (alternative for Tab 2): https://qupath-to-lmd-mdcberlin.streamlit.app/
- py-lmd documentation: https://mannlabs.github.io/py-lmd/
- py-lmd GitHub: https://github.com/MannLabs/py-lmd
