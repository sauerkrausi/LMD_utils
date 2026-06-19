"""
process_lmd_collection_streamlit.py
====================================
Streamlit UI for process_lmd_collection.py

Upload the .zip from the online LMD converter → get sorted XML + CSVs + plate maps.

Deploy: streamlit run process_lmd_collection_streamlit.py
"""

import io
import json
import csv
import re
import zipfile
import xml.etree.ElementTree as ET
import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# ============================================================
# WELL LAYOUT
# ============================================================
ROWS   = list("ABCDEFGH")
COLS   = list(range(1, 13))
WELLS_96 = [f"{r}{c}" for r in ROWS for c in COLS]

def well_sort_key(well_str):
    well_str = well_str.strip()
    if not well_str or len(well_str) < 2:
        return (99, 99)
    try:
        return (ord(well_str[0].upper()) - ord('A'), int(well_str[1:]))
    except ValueError:
        return (99, 99)

def indent_xml(elem, level=0):
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

# ============================================================
# CORE PROCESSING (in-memory)
# ============================================================
def process_collection_data(xml_bytes: bytes, json_bytes: bytes, stem: str) -> dict:
    """
    Returns dict with keys:
      sorted_xml, wellplate_csv, sample_list_csv, updated_json (all bytes)
      grid (dict for plate visualization)
      n_rois (int)
      warnings (list of str)
    """
    warnings = []

    sample_to_orig_well = json.loads(json_bytes.decode('utf-8'))
    orig_well_to_sample = {v: k for k, v in sample_to_orig_well.items()}

    all_samples_alpha = sorted(sample_to_orig_well.keys(), key=str.casefold)
    if len(all_samples_alpha) > 96:
        warnings.append(f"{len(all_samples_alpha)} samples exceed 96-well capacity; extras omitted.")

    sample_to_new_well = {
        sample: WELLS_96[i]
        for i, sample in enumerate(all_samples_alpha)
        if i < 96
    }

    # Parse XML
    tree = ET.parse(io.StringIO(xml_bytes.decode('utf-8-sig')))
    root = tree.getroot()

    shape_pat  = re.compile(r'^Shape_\d+$')
    shapes     = [el for el in root if shape_pat.match(el.tag)]
    non_shapes = [el for el in root if not shape_pat.match(el.tag)]

    def shape_new_well_key(el):
        orig_cap = el.find('CapID').text.strip()
        sample   = orig_well_to_sample.get(orig_cap, "")
        new_well = sample_to_new_well.get(sample, "Z99")
        return well_sort_key(new_well)

    shapes_sorted = sorted(shapes, key=shape_new_well_key)

    for el in list(root):
        root.remove(el)
    for el in non_shapes:
        root.append(el)

    sc = root.find('ShapeCount')
    if sc is not None:
        sc.text = str(len(shapes_sorted))

    for new_idx, el in enumerate(shapes_sorted, start=1):
        el.tag = f"Shape_{new_idx}"
        orig_cap    = el.find('CapID').text.strip()
        sample_name = orig_well_to_sample.get(orig_cap, "")
        new_well    = sample_to_new_well.get(sample_name, "")

        cap_el      = el.find('CapID')
        cap_el.text = new_well

        children    = list(el)
        cap_pos     = children.index(cap_el)
        transfer_el = ET.Element('TransferID')
        transfer_el.text = sample_name
        el.insert(cap_pos, transfer_el)
        root.append(el)

    indent_xml(root)
    xml_buf = io.BytesIO()
    tree.write(xml_buf, encoding='UTF-8', xml_declaration=True)
    sorted_xml_bytes = xml_buf.getvalue()

    # 96-well plate CSV
    grid = {r: {c: "" for c in COLS} for r in ROWS}
    for sample, well in sample_to_new_well.items():
        grid[well[0]][int(well[1:])] = sample

    wellplate_buf = io.StringIO()
    w = csv.writer(wellplate_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [grid[r][c] for c in COLS])

    # Sample list CSV
    sample_list_buf = io.StringIO()
    w = csv.writer(sample_list_buf)
    w.writerow(["ROI_number", "TransferID / Sample Name", "well_ID", "comments", "processed"])
    for roi_num, el in enumerate(shapes_sorted, start=1):
        sample  = el.find('TransferID').text or ""
        new_well = el.find('CapID').text or ""
        w.writerow([roi_num, sample, new_well, "", ""])

    # Updated JSON
    updated_json_str = json.dumps(sample_to_new_well, indent=4, ensure_ascii=False)

    return {
        "sorted_xml":       sorted_xml_bytes,
        "wellplate_csv":    wellplate_buf.getvalue().encode('utf-8'),
        "sample_list_csv":  sample_list_buf.getvalue().encode('utf-8'),
        "updated_json":     updated_json_str.encode('utf-8'),
        "grid":             grid,
        "sample_to_new_well": sample_to_new_well,
        "n_rois":           len(shapes_sorted),
        "n_samples":        len(all_samples_alpha),
        "warnings":         warnings,
    }

# ============================================================
# PLATE VISUALIZATION
# ============================================================
def plot_plate_png(grid, title) -> bytes:
    """Render 96-well plate colored by sample group prefix. Returns PNG bytes."""
    # Extract group label (everything before last underscore segment for coloring)
    all_labels = sorted({grid[r][c] for r in ROWS for c in COLS if grid[r][c]})
    # Color by unique prefix (first token before first underscore) for visual grouping
    def group_of(label):
        return label.split("_")[0] if label else ""

    groups = sorted({group_of(l) for l in all_labels if l})
    palette = cm.tab20
    group_color = {g: palette(i / max(len(groups), 1)) for i, g in enumerate(groups)}

    fig, ax = plt.subplots(figsize=(14, 7))
    ax.set_xlim(-0.8, 12.5)
    ax.set_ylim(-0.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)

    for r_idx, r in enumerate(ROWS):
        for c_idx, c in enumerate(COLS):
            x     = c_idx
            y     = 7 - r_idx
            label = grid[r][c]
            color = group_color.get(group_of(label), "whitesmoke") if label else "whitesmoke"
            edge  = "#999999" if not label else "#333333"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.7, zorder=2))
            if label:
                short = label.replace("_", "\n")
                ax.text(x, y, short, ha="center", va="center",
                        fontsize=4, zorder=3, color="black")

    for r_idx, r in enumerate(ROWS):
        ax.text(-0.65, 7 - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")

    patches = [mpatches.Patch(color=group_color[g], label=g) for g in groups]
    if patches:
        ax.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc="upper left",
                  fontsize=7, title="Sample group", title_fontsize=8)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()

# ============================================================
# STREAMLIT UI
# ============================================================
st.set_page_config(page_title="LMD Collection Processor", layout="wide")

# Footer
st.markdown(
    "<div style='position:fixed;bottom:10px;right:16px;font-size:12px;color:#888;'>"
    "Questions / issues: <a href='https://github.com/sauerkrausi/LMD_utils' target='_blank'>"
    "github.com/sauerkrausi/LMD_utils</a></div>",
    unsafe_allow_html=True
)
st.title("LMD Collection Processor")
st.markdown("Upload the **`.zip`** exported from the online LMD converter. "
            "The app sorts ROIs alphabetically, assigns new well positions, "
            "and generates all output files for download.")

uploaded = st.file_uploader("Upload collection .zip", type=["zip"])

if uploaded:
    with zipfile.ZipFile(uploaded) as z:
        names     = z.namelist()
        xml_names = [n for n in names if n.endswith(".xml") and "_sorted" not in n
                     and not n.startswith("__MACOSX")]
        json_names = [n for n in names if n.endswith("samples_and_wells.json")]

    if not xml_names:
        st.error("No XML file found in the zip (expected one *.xml file).")
        st.stop()
    if not json_names:
        st.error("samples_and_wells.json not found in the zip.")
        st.stop()

    xml_name  = xml_names[0]
    json_name = json_names[0]
    stem      = xml_name.split("/")[-1].replace(".xml", "")

    col1, col2 = st.columns(2)
    col1.success(f"XML: `{xml_name.split('/')[-1]}`")
    col2.success(f"JSON: `{json_name.split('/')[-1]}`")

    if st.button("Process", type="primary"):
        with st.spinner("Processing..."):
            with zipfile.ZipFile(uploaded) as z:
                xml_bytes  = z.read(xml_name)
                json_bytes = z.read(json_name)

            results = process_collection_data(xml_bytes, json_bytes, stem)

        for w in results["warnings"]:
            st.warning(w)

        st.success(f"Done — {results['n_rois']} ROIs / {results['n_samples']} samples")

        # Plate visualization
        st.subheader("96-well plate layout")
        png_bytes = plot_plate_png(results["grid"], f"{stem} — alphabetical well assignment")
        st.image(png_bytes, use_container_width=True)

        # Downloads
        st.subheader("Downloads")
        c1, c2, c3, c4 = st.columns(4)
        c1.download_button("Sorted XML",        results["sorted_xml"],
                           file_name=f"{stem}_sorted.xml",       mime="application/xml")
        c2.download_button("96-well plate CSV", results["wellplate_csv"],
                           file_name=f"{stem}_96wellplate.csv",  mime="text/csv")
        c3.download_button("Sample list CSV",   results["sample_list_csv"],
                           file_name=f"{stem}_sample_list.csv",  mime="text/csv")
        c4.download_button("Updated JSON",      results["updated_json"],
                           file_name="samples_and_wells_updated.json", mime="application/json")

        c5, _ = st.columns([1, 3])
        c5.download_button("Plate map PNG",     png_bytes,
                           file_name=f"{stem}_platemap.png",     mime="image/png")
