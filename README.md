# spatial_misc
collection of scripts and code for working in spatial proteomics and laser micro-dissection microscopy (LMD) workflows

---

## process_lmd_collection.py

Post-processes the output folder from the online LMD converter.

**Input** (must be in the collection folder):
- `*.xml` — shapes with CapID + XY coordinates
- `samples_and_wells.json` — `{ "SampleName": "WellPosition" }`

**Output** (written to `lmd_outputs/`):
- `*_sorted.xml` — shapes sorted by new well; `<TransferID>` added; `<CapID>` updated
- `*_96wellplate.csv` — 96-well layout, samples assigned alphabetically
- `*_sample_list.csv` — ROI number, sample name, well ID, comments, processed
- `samples_and_wells_updated.json` — updated well mapping

**Usage:**


---

## sort_XML_ROI_by96well.py

Sorts an XML file by CapID in 96-well order (A1, A2, ..., H12) and renumbers shapes.

**Input:** XML file with Shape_N elements containing CapID + XY coordinates  
**Output:** `*_sorted.xml` in the same directory

**Usage:** Set `INPUT_PATH` at the top of the script and run.

> Note: `process_lmd_collection.py` supersedes this script for full workflows.