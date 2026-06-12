# Sort ROI XML shapes by CapID in 96-well plate order (A1-A12, B1-B12, ...)
# Input: XML file with Shape_N elements containing CapID and XY coordinates
# Output: _sorted version in the same directory, with shapes renumbered and header intact

import xml.etree.ElementTree as ET
import re
from pathlib import Path

INPUT_PATH = "or can do /path/to/xmlfile.xml"

def capid_sort_key(cap_id):
    row = ord(cap_id[0].upper()) - ord('A')
    col = int(cap_id[1:])
    print(row, col)
    return (row, col)

def sort_rois(input_path):
    input_path = Path(input_path)
    output_path = input_path.with_name(input_path.stem + "_sorted" + input_path.suffix)

    # utf-8-sig handles BOM-prefixed UTF-8 files
    with open(input_path, encoding='utf-8-sig') as f:
        tree = ET.parse(f)
    root = tree.getroot()

    shape_pattern = re.compile(r'^Shape_\d+$')
    shapes = [(el.tag, el) for el in root if shape_pattern.match(el.tag)]
    non_shapes = [(el.tag, el) for el in root if not shape_pattern.match(el.tag)]

    shapes_sorted = sorted(shapes, key=lambda x: capid_sort_key(x[1].find('CapID').text.strip()))

    for el in list(root):
        root.remove(el)

    for tag, el in non_shapes:
        root.append(el)

    for new_idx, (_, el) in enumerate(shapes_sorted, start=1):
        el.tag = f"Shape_{new_idx}"
        root.append(el)

    tree.write(output_path, encoding='UTF-8', xml_declaration=True)
    print(f"Sorted file written to: {output_path}")

if __name__ == "__main__":
    sort_rois(INPUT_PATH)