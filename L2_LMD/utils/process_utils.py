"""
process_utils.py
================
Tab 3 -- Process LMD Collection.

Input (piped from Tab 2 or uploaded as zip):
  - XML bytes
  - samples_and_wells dict {sample_name: well_id}

Steps:
  1. Sort samples alphabetically, assign new wells A1..H12
  2. Sort XML shapes by new well, update CapID, inject TransferID
  3. Generate 96-well plate CSV + plate map PNG
  4. Generate sample list CSV (columns: cut order, ROI, Well_ID, Dropout {Y/N}, comments, processed)

Pipes to Tab 4: sample_list CSV bytes + stem via session_state.t3_sample_list / t3_stem
"""

import io
import json
import csv
import re
import zipfile
import xml.etree.ElementTree as ET

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import pandas as pd
import streamlit as st

ROWS     = list("ABCDEFGH")
COLS     = list(range(1, 13))
WELLS_96 = [f"{r}{c}" for r in ROWS for c in COLS]


# ============================================================
# HELPERS
# ============================================================
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
# CORE LOGIC
# ============================================================
def process_collection(xml_bytes: bytes, saw_dict: dict, stem: str) -> dict:
    """
    All in-memory. Returns dict with output bytes, grid, metadata.
    saw_dict: {sample_name: original_well}
    """
    warnings = []

    orig_well_to_sample = {v: k for k, v in saw_dict.items()}
    all_samples_alpha   = sorted(saw_dict.keys(), key=str.casefold)

    existing_wells = sorted(saw_dict.values(), key=well_sort_key)
    start_well     = existing_wells[0] if existing_wells else "A1"
    start_idx      = WELLS_96.index(start_well) if start_well in WELLS_96 else 0

    if start_idx > 0:
        warnings.append(f"Custom start position detected: wells assigned from {start_well}.")

    available = 96 - start_idx
    if len(all_samples_alpha) > available:
        warnings.append(f"{len(all_samples_alpha)} samples exceed available wells from "
                        f"{start_well} ({available} slots); extras omitted.")

    sample_to_new_well = {
        s: WELLS_96[start_idx + i]
        for i, s in enumerate(all_samples_alpha)
        if start_idx + i < 96
    }

    # Parse + sort XML
    tree       = ET.parse(io.StringIO(xml_bytes.decode("utf-8-sig")))
    root       = tree.getroot()
    shape_pat  = re.compile(r'^Shape_\d+$')
    shapes     = [el for el in root if shape_pat.match(el.tag)]
    non_shapes = [el for el in root if not shape_pat.match(el.tag)]

    def shape_key(el):
        cap = el.find("CapID")
        orig = cap.text.strip() if cap is not None and cap.text else ""
        sample = orig_well_to_sample.get(orig, "")
        return well_sort_key(sample_to_new_well.get(sample, "Z99"))

    shapes_sorted = sorted(shapes, key=shape_key)

    for el in list(root):
        root.remove(el)
    for el in non_shapes:
        root.append(el)

    sc = root.find("ShapeCount")
    if sc is not None:
        sc.text = str(len(shapes_sorted))

    for new_idx, el in enumerate(shapes_sorted, start=1):
        el.tag      = f"Shape_{new_idx}"
        cap_el      = el.find("CapID")
        orig_cap    = cap_el.text.strip() if cap_el is not None and cap_el.text else ""
        sample_name = orig_well_to_sample.get(orig_cap, "")
        new_well    = sample_to_new_well.get(sample_name, "")

        if cap_el is not None:
            cap_el.text = new_well

        children    = list(el)
        cap_pos     = children.index(cap_el) if cap_el is not None else 0
        transfer_el = ET.Element("TransferID")
        transfer_el.text = sample_name
        el.insert(cap_pos, transfer_el)
        root.append(el)

    indent_xml(root)
    xml_buf = io.BytesIO()
    tree.write(xml_buf, encoding="UTF-8", xml_declaration=True)

    # 96-well plate grid + CSV
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
    w.writerow(["cut order", "ROI", "Well_ID", "Dropout {Y/N}", "comments", "processed"])
    for roi_num, el in enumerate(shapes_sorted, start=1):
        t_el   = el.find("TransferID")
        c_el   = el.find("CapID")
        sample = t_el.text if t_el is not None and t_el.text else ""
        well   = c_el.text if c_el is not None and c_el.text else ""
        w.writerow([roi_num, sample, well, "", "", ""])

    updated_json = json.dumps(sample_to_new_well, indent=4, ensure_ascii=False).encode("utf-8")

    return {
        "sorted_xml":      xml_buf.getvalue(),
        "wellplate_csv":   wellplate_buf.getvalue().encode("utf-8"),
        "sample_list_csv": sample_list_buf.getvalue().encode("utf-8"),
        "updated_json":    updated_json,
        "grid":            grid,
        "n_rois":          len(shapes_sorted),
        "n_samples":       len(all_samples_alpha),
        "warnings":        warnings,
        "stem":            stem,
    }


def plot_plate_png(grid: dict, title: str) -> bytes:
    all_labels  = sorted({grid[r][c] for r in ROWS for c in COLS if grid[r][c]})
    groups      = sorted({l.split("_")[0] for l in all_labels})
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
            grp   = label.split("_")[0] if label else ""
            color = group_color.get(grp, "whitesmoke") if label else "whitesmoke"
            edge  = "#999999" if not label else "#333333"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.7, zorder=2))
            if label:
                ax.text(x, y, label.replace("_", "\n"), ha="center", va="center",
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


def build_zip(results: dict, png_bytes: bytes) -> bytes:
    stem = results["stem"]
    buf  = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}_sorted.xml",             results["sorted_xml"])
        z.writestr(f"{stem}_96wellplate.csv",        results["wellplate_csv"])
        z.writestr(f"{stem}_sample_list.csv",        results["sample_list_csv"])
        z.writestr("samples_and_wells_updated.json", results["updated_json"])
        z.writestr(f"{stem}_platemap.png",           png_bytes)
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_process_tab():
    st.header("Process LMD Collection")
    st.caption(
        "Sorts ROIs alphabetically, assigns wells A1..H12, updates XML CapID, "
        "generates plate map and sample list. "
        "Accepts output piped from Tab 2 or a zip from the "
        "[Coscia Lab converter](https://qupath-to-lmd-mdcberlin.streamlit.app/) (fallback)."
    )

    # Source: piped from Tab 2 or zip upload
    t2_xml = st.session_state.get("t2_xml")
    t2_saw = st.session_state.get("t2_saw")
    t2_stem = st.session_state.get("t2_stem")

    source = None

    if t2_xml is not None and t2_saw is not None:
        st.info("Using XML and well map piped from Tab 2.")
        use_pipe = st.checkbox("Use piped data from Tab 2", value=True)
        if use_pipe:
            source = "pipe"

    uploaded_zip = st.file_uploader(
        "Or upload collection zip (from Coscia Lab converter or Tab 2 download)",
        type=["zip"], key="proc_upload"
    )

    if uploaded_zip is not None:
        source = "zip"

    if source is None:
        st.stop()

    # Extract bytes
    if source == "pipe":
        xml_bytes = t2_xml
        saw_dict  = t2_saw
        stem      = t2_stem or "collection"
    else:
        with zipfile.ZipFile(uploaded_zip) as z:
            names      = z.namelist()
            xml_names  = [n for n in names if n.endswith(".xml")
                          and "_sorted" not in n and not n.startswith("__MACOSX")]
            json_names = [n for n in names if n.endswith("samples_and_wells.json")]

        if not xml_names:
            st.error("No XML file found in zip.")
            return
        if not json_names:
            st.error("samples_and_wells.json not found in zip.")
            return

        with zipfile.ZipFile(uploaded_zip) as z:
            xml_bytes = z.read(xml_names[0])
            saw_dict  = json.loads(z.read(json_names[0]).decode("utf-8"))

        stem = xml_names[0].split("/")[-1].replace(".xml", "")

    # Clear stale results when input changes
    input_key = (id(xml_bytes), stem)
    if st.session_state.get("proc_last_key") != input_key:
        st.session_state.proc_result = None
        st.session_state.proc_last_key = input_key

    st.metric("Samples in JSON", len(saw_dict))

    if st.button("Process", type="primary"):
        with st.spinner("Processing..."):
            result    = process_collection(xml_bytes, saw_dict, stem)
            png_bytes = plot_plate_png(result["grid"], f"{stem} - well assignment")
            zip_bytes = build_zip(result, png_bytes)

        st.session_state.proc_result    = result
        st.session_state.proc_png       = png_bytes
        st.session_state.proc_zip       = zip_bytes
        st.session_state.proc_last_key  = input_key

        # Pipe to Tab 4
        st.session_state.t3_sample_list = result["sample_list_csv"]
        st.session_state.t3_stem        = stem

        st.success(f"{result['n_rois']} ROIs processed.")

    result = st.session_state.get("proc_result")
    if result:
        for w in result["warnings"]:
            st.warning(w)

        st.subheader("96-well plate layout")
        st.image(st.session_state.proc_png, use_container_width=True)

        st.subheader("Downloads")
        stem_out = result["stem"]

        dl0, dl1, dl2, dl3, dl4 = st.columns(5)
        dl0.download_button("All (zip)", st.session_state.proc_zip,
                            file_name=f"{stem_out}_lmd_outputs.zip",
                            mime="application/zip", type="primary")
        dl1.download_button("Sorted XML", result["sorted_xml"],
                            file_name=f"{stem_out}_sorted.xml", mime="application/xml")
        dl2.download_button("96-well CSV", result["wellplate_csv"],
                            file_name=f"{stem_out}_96wellplate.csv", mime="text/csv")
        dl3.download_button("Sample list", result["sample_list_csv"],
                            file_name=f"{stem_out}_sample_list.csv", mime="text/csv")
        dl4.download_button("Plate PNG", st.session_state.proc_png,
                            file_name=f"{stem_out}_platemap.png", mime="image/png")

        st.caption("Sample list is available in Tab 4 (MS Queue) without re-uploading.")
