"""
process_lmd_collection_streamlit.py
====================================
Two-tab Streamlit app for the LMD workflow.

TAB 1 — Re-classify GeoJSON
  Upload a QuPath GeoJSON export, copy properties.name into
  properties.classification.name for each feature, download corrected
  GeoJSON for use with the Coscia Lab online converter.

TAB 2 — Process LMD Collection
  Upload the .zip from the Coscia Lab QuPath-to-LMD online converter.
  Sorts ROIs alphabetically, assigns new well positions A1..H12,
  updates CapID, injects TransferID, outputs all files for download.

Deploy:  streamlit run process_lmd_collection_streamlit.py
"""

import io
import json
import csv
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# ============================================================
# PAGE CONFIG
# ============================================================
st.set_page_config(page_title="LMD Collection Processor", layout="wide")

# Fixed footer with GitHub link
st.markdown(
    "<div style='position:fixed;bottom:10px;right:16px;font-size:12px;color:#888;'>"
    "Questions / issues: <a href='https://github.com/sauerkrausi/LMD_utils' target='_blank'>"
    "github.com/sauerkrausi/LMD_utils</a></div>",
    unsafe_allow_html=True
)

# ============================================================
# WELL LAYOUT CONSTANTS
# ============================================================
ROWS     = list("ABCDEFGH")
COLS     = list(range(1, 13))
WELLS_96 = [f"{r}{c}" for r in ROWS for c in COLS]  # A1, A2, ..., H12

# ============================================================
# HELPERS
# ============================================================
def well_sort_key(well_str):
    """Convert well ID (e.g. 'B3') to a sortable (row, col) tuple."""
    well_str = well_str.strip()
    if not well_str or len(well_str) < 2:
        return (99, 99)
    try:
        return (ord(well_str[0].upper()) - ord('A'), int(well_str[1:]))
    except ValueError:
        return (99, 99)

def indent_xml(elem, level=0):
    """Recursively add pretty-print indentation to an XML element tree."""
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
# CORE PROCESSING  (all in-memory, no disk I/O)
# ============================================================
def process_collection_data(xml_bytes: bytes, json_bytes: bytes, stem: str) -> dict:
    """
    Takes raw file bytes, returns a dict of output bytes + metadata.

    Steps:
      1. Parse JSON mapping (sample -> original well)
      2. Build alphabetical well assignments (sample -> new well A1..H12)
      3. Parse XML, sort shapes by new well, update CapID, inject TransferID
      4. Serialise sorted XML to bytes
      5. Build 96-well plate CSV
      6. Build sample list CSV (one row per ROI in XML order)
      7. Build updated JSON with new well assignments
    """
    warnings = []

    # 1. Load sample -> original well mapping and invert it
    sample_to_orig_well  = json.loads(json_bytes.decode('utf-8'))
    orig_well_to_sample  = {v: k for k, v in sample_to_orig_well.items()}

    # 2. Alphabetical well assignment
    # Respect custom start position from JSON (e.g. E1 instead of A1)
    all_samples_alpha = sorted(sample_to_orig_well.keys(), key=str.casefold)

    existing_wells = sorted(sample_to_orig_well.values(), key=well_sort_key)
    start_well     = existing_wells[0] if existing_wells else "A1"
    start_idx      = WELLS_96.index(start_well) if start_well in WELLS_96 else 0

    if start_idx > 0:
        warnings.append(f"Custom start position detected: wells assigned from {start_well} "
                        f"(as defined in samples_and_wells.json).")

    available = 96 - start_idx
    if len(all_samples_alpha) > available:
        warnings.append(f"{len(all_samples_alpha)} samples exceed available wells "
                        f"from {start_well} ({available} slots); extras omitted.")

    sample_to_new_well = {
        sample: WELLS_96[start_idx + i]
        for i, sample in enumerate(all_samples_alpha)
        if start_idx + i < 96
    }

    # 3. Parse XML and sort shapes
    tree = ET.parse(io.StringIO(xml_bytes.decode('utf-8-sig')))
    root = tree.getroot()

    shape_pat  = re.compile(r'^Shape_\d+$')
    shapes     = [el for el in root if shape_pat.match(el.tag)]
    non_shapes = [el for el in root if not shape_pat.match(el.tag)]

    def shape_new_well_key(el):
        orig_cap = el.find('CapID').text.strip()
        sample   = orig_well_to_sample.get(orig_cap, "")
        new_well = sample_to_new_well.get(sample, "Z99")  # unmapped shapes sort last
        return well_sort_key(new_well)

    shapes_sorted = sorted(shapes, key=shape_new_well_key)

    # Rebuild XML root: header elements first, then renumbered shapes
    for el in list(root):
        root.remove(el)
    for el in non_shapes:
        root.append(el)

    sc = root.find('ShapeCount')
    if sc is not None:
        sc.text = str(len(shapes_sorted))

    for new_idx, el in enumerate(shapes_sorted, start=1):
        el.tag       = f"Shape_{new_idx}"
        orig_cap     = el.find('CapID').text.strip()
        sample_name  = orig_well_to_sample.get(orig_cap, "")
        new_well     = sample_to_new_well.get(sample_name, "")

        # Overwrite CapID with new well position
        cap_el      = el.find('CapID')
        cap_el.text = new_well

        # Insert <TransferID> immediately before <CapID>
        children    = list(el)
        cap_pos     = children.index(cap_el)
        transfer_el = ET.Element('TransferID')
        transfer_el.text = sample_name
        el.insert(cap_pos, transfer_el)
        root.append(el)

    # 4. Serialise sorted XML
    indent_xml(root)
    xml_buf = io.BytesIO()
    tree.write(xml_buf, encoding='UTF-8', xml_declaration=True)

    # 5. 96-well plate CSV
    grid = {r: {c: "" for c in COLS} for r in ROWS}
    for sample, well in sample_to_new_well.items():
        grid[well[0]][int(well[1:])] = sample

    wellplate_buf = io.StringIO()
    w = csv.writer(wellplate_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [grid[r][c] for c in COLS])

    # 6. Sample list CSV (ROI order matches XML shape order)
    sample_list_buf = io.StringIO()
    w = csv.writer(sample_list_buf)
    w.writerow(["cut order", "ROI", "Well_ID", "Dropout {Y/N}", "comments", "processed"])
    for roi_num, el in enumerate(shapes_sorted, start=1):
        sample   = el.find('TransferID').text or ""
        new_well = el.find('CapID').text or ""
        w.writerow([roi_num, sample, new_well, "", "", ""])

    # 7. Updated JSON
    updated_json_str = json.dumps(sample_to_new_well, indent=4, ensure_ascii=False)

    return {
        "sorted_xml":         xml_buf.getvalue(),
        "wellplate_csv":      wellplate_buf.getvalue().encode('utf-8'),
        "sample_list_csv":    sample_list_buf.getvalue().encode('utf-8'),
        "updated_json":       updated_json_str.encode('utf-8'),
        "grid":               grid,
        "n_rois":             len(shapes_sorted),
        "n_samples":          len(all_samples_alpha),
        "warnings":           warnings,
        "stem":               stem,
    }

# ============================================================
# PLATE VISUALIZATION
# ============================================================
def plot_plate_png(grid, title) -> bytes:
    """
    Renders a 96-well plate as a PNG.
    Wells are colored by sample group (first token before '_').
    Returns PNG as bytes.
    """
    all_labels = sorted({grid[r][c] for r in ROWS for c in COLS if grid[r][c]})

    def group_of(label):
        return label.split("_")[0] if label else ""

    groups      = sorted({group_of(l) for l in all_labels if l})
    palette     = cm.tab20
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
                ax.text(x, y, label.replace("_", "\n"), ha="center", va="center",
                        fontsize=4, zorder=3, color="black")

    # Row (A-H) and column (1-12) axis labels
    for r_idx, r in enumerate(ROWS):
        ax.text(-0.65, 7 - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")

    # Legend: one entry per sample group
    patches = [mpatches.Patch(color=group_color[g], label=g) for g in groups]
    if patches:
        ax.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc="upper left",
                  fontsize=7, title="Sample group", title_fontsize=8)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()

def build_download_zip(results: dict, png_bytes: bytes) -> bytes:
    """Bundle all output files into a single zip for 'Download all'."""
    stem = results["stem"]
    buf  = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}_sorted.xml",              results["sorted_xml"])
        z.writestr(f"{stem}_96wellplate.csv",         results["wellplate_csv"])
        z.writestr(f"{stem}_sample_list.csv",         results["sample_list_csv"])
        z.writestr("samples_and_wells_updated.json",  results["updated_json"])
        z.writestr(f"{stem}_platemap.png",            png_bytes)
    return buf.getvalue()

# ============================================================
# GEOJSON RE-CLASSIFIER
# ============================================================
def reclassify_geojson(geojson_bytes: bytes) -> tuple:
    """
    For each feature: copy properties.name -> properties.classification.name.
    Returns (corrected_bytes, n_fixed, name_list, rejected_list).
    """
    geojson  = json.loads(geojson_bytes.decode('utf-8'))
    n_fixed  = 0
    names    = []
    rejected = []

    for feature in geojson.get("features", []):
        props = feature.get("properties", {})
        name  = props.get("name", "").strip()
        if not name:
            rejected.append("(unnamed feature)")
            continue
        classification = props.setdefault("classification", {})
        if classification.get("name") != name:
            classification["name"] = name
            n_fixed += 1
        names.append(name)

    corrected = json.dumps(geojson, indent=2).encode('utf-8')
    return corrected, n_fixed, names, rejected

# ============================================================
# STREAMLIT UI
# ============================================================
st.title("LMD Utils")

# ============================================================
# SECTION 1 — Re-classify GeoJSON
# ============================================================
st.header("1. Re-classify GeoJSON")
if True:
    st.markdown(
        "Upload a QuPath GeoJSON export. The app copies **`properties.name`** into "
        "**`properties.classification.name`** for each ROI. Re-classification of annotation-names to classes of the same name isrequired for the "
        "[Coscia Lab QuPath to XML converter](https://qupath-to-lmd-mdcberlin.streamlit.app/)."
    )

    geo_uploaded = st.file_uploader("Upload GeoJSON", type=["geojson", "json"], key="geo_up")

    # Session state for tab 1
    for k in ("geo_result", "geo_last"):
        if k not in st.session_state:
            st.session_state[k] = None

    if geo_uploaded and geo_uploaded.name != st.session_state.geo_last:
        st.session_state.geo_result = None
        st.session_state.geo_last   = geo_uploaded.name

    if geo_uploaded:
        if st.button("Re-classify", type="primary", key="geo_btn"):
            with st.spinner("Processing..."):
                corrected, n_fixed, names, rejected = reclassify_geojson(geo_uploaded.getvalue())
                st.session_state.geo_result = {
                    "bytes":    corrected,
                    "n_fixed":  n_fixed,
                    "names":    names,
                    "rejected": rejected,
                    "stem":     Path(geo_uploaded.name).stem,
                }

    if st.session_state.geo_result:
        r = st.session_state.geo_result
        st.success(f"Done — {len(r['names'])} ROIs, {r['n_fixed']} classification(s) updated.")

        if r["rejected"]:
            st.warning(f"{len(r['rejected'])} feature(s) skipped (no name): {r['rejected']}")

        with st.expander(f"ROI names ({len(r['names'])})", expanded=False):
            st.write(sorted(r["names"]))

        st.download_button(
            "⬇ Download corrected GeoJSON",
            r["bytes"],
            file_name=f"{r['stem']}_reclassified.geojson",
            mime="application/geo+json",
            type="primary",
        )

        st.info("Next step: upload the corrected GeoJSON to the "
                "[Coscia Lab converter](https://qupath-to-lmd-mdcberlin.streamlit.app/), "
                "then bring the downloaded zip to Step 2.")

st.divider()

# ============================================================
# SECTION 2 — Process LMD Collection
# ============================================================
st.header("2. Process LMD Collection")
if True:
    st.markdown(
        "Upload the **`.zip`** exported from the "
        "[Coscia Lab QuPath to XML converter](https://qupath-to-lmd-mdcberlin.streamlit.app/). "
        "The app sorts ROIs alphabetically, assigns new well positions, "
        "and generates all output files for download."
    )

    uploaded = st.file_uploader("Upload collection .zip", type=["zip"], key="lmd_up")

    # Session state for tab 2
    for k in ("results", "png", "zip_all", "last_file"):
        if k not in st.session_state:
            st.session_state[k] = None

    if uploaded and uploaded.name != st.session_state.last_file:
        st.session_state.results   = None
        st.session_state.png       = None
        st.session_state.zip_all   = None
        st.session_state.last_file = uploaded.name

    if uploaded:
        with zipfile.ZipFile(uploaded) as z:
            names      = z.namelist()
            xml_names  = [n for n in names if n.endswith(".xml")
                          and "_sorted" not in n and not n.startswith("__MACOSX")]
            json_names = [n for n in names if n.endswith("samples_and_wells.json")]

        if not xml_names:
            st.error("No XML file found in the zip.")
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

        if st.button("Process", type="primary", key="lmd_btn"):
            with st.spinner("Processing..."):
                with zipfile.ZipFile(uploaded) as z:
                    xml_bytes  = z.read(xml_name)
                    json_bytes = z.read(json_name)

                results  = process_collection_data(xml_bytes, json_bytes, stem)
                png      = plot_plate_png(results["grid"], f"{stem} — well assignment")
                zip_all  = build_download_zip(results, png)

                st.session_state.results  = results
                st.session_state.png      = png
                st.session_state.zip_all  = zip_all

    if st.session_state.results:
        results = st.session_state.results
        stem    = results["stem"]

        for w in results["warnings"]:
            st.warning(w)

        st.success(f"Done — {results['n_rois']} ROIs / {results['n_samples']} samples")

        st.subheader("96-well plate layout")
        st.image(st.session_state.png, use_container_width=True)

        st.subheader("Downloads")
        st.download_button(
            "⬇ Download all (zip)",
            st.session_state.zip_all,
            file_name=f"{stem}_lmd_outputs.zip",
            mime="application/zip",
            type="primary",
        )

        st.markdown("---")

        c1, c2, c3, c4, c5 = st.columns(5)
        c1.download_button("Sorted XML",        results["sorted_xml"],
                           file_name=f"{stem}_sorted.xml",      mime="application/xml")
        c2.download_button("96-well plate CSV", results["wellplate_csv"],
                           file_name=f"{stem}_96wellplate.csv", mime="text/csv")
        c3.download_button("Sample list CSV",   results["sample_list_csv"],
                           file_name=f"{stem}_sample_list.csv", mime="text/csv")
        c4.download_button("Updated JSON",      results["updated_json"],
                           file_name="samples_and_wells_updated.json", mime="application/json")
        c5.download_button("Plate map PNG",     st.session_state.png,
                           file_name=f"{stem}_platemap.png",    mime="image/png")
