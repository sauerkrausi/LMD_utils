#!/usr/bin/env python3
"""
process_lmd_collection.py
=========================
Post-processes the output folder from the online LMD converter.

REQUIRED INPUT FILES (must be inside the collection folder)
------------------------------------------------------------
  *.xml                  XML file exported by the online converter
                         (shapes with CapID + XY coordinates)
  samples_and_wells.json JSON mapping  { "SampleName": "WellPosition" }
                         e.g. { "Tumor_ROI1": "A1", ... }

WHAT IT DOES
------------
  1. Sorts samples alphabetically, assigns new wells A1, A2, ..., H12.
     Sorts XML shapes by new well order and updates each <CapID> to the
     new well. Injects <TransferID> (= sample name) before <CapID>.
  2. Generates a 96-well plate CSV with samples sorted alphabetically
     (A1 = first alphabetically, A2 = second, ...).
  3. Generates a sample-list CSV with columns:
       ROI_number | TransferID | well_ID | comments | processed
  4. Generates an updated samples_and_wells.json reflecting the new
     alphabetical well assignments.

OUTPUTS  (written to <collection_folder>/lmd_outputs/)
------------------------------------------------------
  *_sorted.xml
  *_96wellplate.csv
  *_sample_list.csv
  samples_and_wells_updated.json

USAGE
-----
  python process_lmd_collection.py <path/to/collection_folder>
  python process_lmd_collection.py          # uses COLLECTION_FOLDER below
"""

import sys
import json
import csv
import re
import xml.etree.ElementTree as ET
from pathlib import Path

# ---- CONFIG (used when no CLI arg is given) ----
COLLECTION_FOLDER = Path("/Users/felix/Desktop/tmpLMD/VG/ID-0842 - 1875 A2_collection")
OUTPUT_SUBDIR = "lmd_outputs"

# ---- WELL LAYOUT ----
ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))
WELLS_96 = [f"{r}{c}" for r in ROWS for c in COLS]  # A1, A2, ..., H12

def well_sort_key(well_str):
    well_str = well_str.strip()
    return (ord(well_str[0].upper()) - ord('A'), int(well_str[1:]))

# ---- XML HELPERS ----
def indent_xml(elem, level=0):
    """Add pretty-print indentation in-place (pure stdlib, works on Python 3.8+)."""
    pad = "\n" + "  " * level
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = pad + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = pad
        for child in elem:
            indent_xml(child, level + 1)
        if not child.tail or not child.tail.strip():
            child.tail = pad
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = pad
    if not level:
        elem.tail = "\n"

# ---- MAIN ----
def process_collection(folder: Path):
    folder = Path(folder).resolve()
    out_dir = folder / OUTPUT_SUBDIR
    out_dir.mkdir(exist_ok=True)

    # Locate files
    xml_files = [f for f in folder.glob("*.xml") if "_sorted" not in f.stem]
    json_file = folder / "samples_and_wells.json"

    if not xml_files:
        raise FileNotFoundError("No XML file found in folder (excluding *_sorted.xml)")
    if not json_file.exists():
        raise FileNotFoundError("samples_and_wells.json not found in folder")

    xml_path = xml_files[0]
    stem = xml_path.stem
    print(f"Input XML : {xml_path.name}")
    print(f"Input JSON: {json_file.name}")
    print(f"Output dir: {out_dir}")

    # Load JSON: sample_name -> original_well
    with open(json_file, encoding='utf-8') as f:
        sample_to_orig_well = json.load(f)

    # Invert: original_well -> sample_name
    orig_well_to_sample = {v: k for k, v in sample_to_orig_well.items()}

    # ------------------------------------------------------------------ #
    # Build alphabetical well assignments first (needed for XML sort)     #
    # ------------------------------------------------------------------ #
    all_samples_alpha = sorted(sample_to_orig_well.keys(), key=str.casefold)

    if len(all_samples_alpha) > 96:
        print(f"WARNING: {len(all_samples_alpha)} samples exceed 96-well capacity; extras omitted.")

    sample_to_new_well = {
        sample: WELLS_96[i]
        for i, sample in enumerate(all_samples_alpha)
        if i < 96
    }

    # ------------------------------------------------------------------ #
    # 1. Sort XML by new well, update CapID, inject TransferID            #
    # ------------------------------------------------------------------ #
    with open(xml_path, encoding='utf-8-sig') as f:
        tree = ET.parse(f)
    root = tree.getroot()

    shape_pat = re.compile(r'^Shape_\d+$')
    shapes     = [el for el in root if shape_pat.match(el.tag)]
    non_shapes = [el for el in root if not shape_pat.match(el.tag)]

    def shape_new_well_key(el):
        orig_cap = el.find('CapID').text.strip()
        sample = orig_well_to_sample.get(orig_cap, "")
        new_well = sample_to_new_well.get(sample, "Z99")  # unmapped shapes sort last
        return well_sort_key(new_well)

    shapes_sorted = sorted(shapes, key=shape_new_well_key)

    # Rebuild root: non-shapes first, then renumbered shapes
    for el in list(root):
        root.remove(el)
    for el in non_shapes:
        root.append(el)

    # Update ShapeCount
    sc = root.find('ShapeCount')
    if sc is not None:
        sc.text = str(len(shapes_sorted))

    for new_idx, el in enumerate(shapes_sorted, start=1):
        el.tag = f"Shape_{new_idx}"
        orig_cap = el.find('CapID').text.strip()
        sample_name = orig_well_to_sample.get(orig_cap, "")
        new_well = sample_to_new_well.get(sample_name, "")

        # Update CapID to new well position
        cap_el = el.find('CapID')
        cap_el.text = new_well

        # Insert <TransferID> immediately before <CapID>
        children = list(el)
        cap_pos = children.index(cap_el)
        transfer_el = ET.Element('TransferID')
        transfer_el.text = sample_name
        el.insert(cap_pos, transfer_el)

        root.append(el)

    indent_xml(root)
    sorted_xml_path = out_dir / (stem + "_sorted.xml")
    tree.write(sorted_xml_path, encoding='UTF-8', xml_declaration=True)
    print(f"\n[1] Sorted XML written : {sorted_xml_path.name}")

    # ------------------------------------------------------------------ #
    # 2. 96-well plate CSV — samples sorted alphabetically                #
    # ------------------------------------------------------------------ #

    # Build grid
    grid = {r: {c: "" for c in COLS} for r in ROWS}
    for sample, well in sample_to_new_well.items():
        row = well[0]
        col = int(well[1:])
        grid[row][col] = sample

    wellplate_path = out_dir / (stem + "_96wellplate.csv")
    with open(wellplate_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow([""] + COLS)
        for r in ROWS:
            writer.writerow([r] + [grid[r][c] for c in COLS])
    print(f"[2] 96-well plate CSV  : {wellplate_path.name}")

    # ------------------------------------------------------------------ #
    # 3. Sample list CSV — one row per ROI in XML sort order              #
    # ------------------------------------------------------------------ #
    sample_list_path = out_dir / (stem + "_sample_list.csv")
    with open(sample_list_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["ROI_number", "TransferID / Sample Name", "well_ID", "comments", "processed"])
        for roi_num, el in enumerate(shapes_sorted, start=1):
            sample = el.find('TransferID').text or ""
            new_well = el.find('CapID').text or ""
            writer.writerow([roi_num, sample, new_well, "", ""])
    print(f"[3] Sample list CSV    : {sample_list_path.name}")

    # ------------------------------------------------------------------ #
    # 4. Updated samples_and_wells.json — new alphabetical well positions  #
    # ------------------------------------------------------------------ #
    updated_json_path = out_dir / "samples_and_wells_updated.json"
    with open(updated_json_path, 'w', encoding='utf-8') as f:
        json.dump(sample_to_new_well, f, indent=4, ensure_ascii=False)
    print(f"[4] Updated JSON       : {updated_json_path.name}")

    print(f"\nDone. {len(shapes_sorted)} ROIs processed.")

if __name__ == "__main__":
    folder = Path(sys.argv[1]) if len(sys.argv) > 1 else COLLECTION_FOLDER
    process_collection(folder)
