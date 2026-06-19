[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://lmdutils-qfyodkmxyhworjpq24nanz.streamlit.app/)

# LMD_utils
Collection of scripts and code for working in spatial proteomics and laser micro-dissection microscopy (LMD) workflows. Might be useful for sample sorting before laser-dissection and for general documentation purposes.

---

## process_lmd_collection.py

Post-processes the output folder from the online Coscia Lab QuPath to XML converter.

**Input** (must be in the collection zip folder):
- `*.xml` - shapes with CapID + XY coordinates
- `samples_and_wells.json` - `{ "SampleName": "WellPosition" }`

**Output** (written to `lmd_outputs/` or available for download when using online app):
- `*_sorted.xml` - shapes sorted by new well; `<TransferID>` added; `<CapID>` updated
- `*_96wellplate.csv` - 96-well layout, samples assigned alphabetically
- `*_sample_list.csv` - ROI number, sample name, well ID, comments, processed
- `samples_and_wells_updated.json` - updated well mapping
- `platemap.png` - color coded plate map 

**Usage:**
1. use CosciaLab QuPath to LMD online tool to convert geojson file into xml file   
Link: (https://qupath-to-lmd-mdcberlin.streamlit.app/)
2. Check that calibration points and QCs are all working. 
3. Make sure to select 96 well plate. No need to update the well map, since the python script will sort by the name of the QuPath annotation.
Confirm plate layout  >> process file >> download zipped filed
4. Extract zipped folder and copy folder path.
5. Add folder path into python file under COLLECTION_FOLDER
6. Save and run python `process_lmd_collection.py` script
7. New subfolder called `lmd_outputs` will be created containing
> - `*_96wellplate.csv`:  96 well format with sample names, sorted in alphabetical order
> - `*_sample_list.csv:  csv file containing samples in well order (A1-H12) as well as their sample name (== TransferID)
> - `*_sorted.xml`: well order sorted (A1-H12) xml file for import to LMD microscope
> - `samples_and_wells_updated.json`: json format of 96 well format




---

## sort_XML_ROI_by96well.py

Sorts an XML file by CapID in 96-well order (A1, A2, ..., H12) and renumbers shapes.

**Input:** XML file with Shape_N elements containing CapID + XY coordinates  
**Output:** `*_sorted.xml` in the same directory

**Usage:** Set `INPUT_PATH` at the top of the script and run.

> Note: `process_lmd_collection.py` supersedes this script for full workflows.

---
### Links
Coscia Lab QuPath to XML github repo: https://github.com/CosciaLab/Qupath_to_LMD
