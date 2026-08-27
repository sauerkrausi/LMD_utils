"""
process_utils.py
================
Tab 3 -- Process LMD Collection (multi-plate aware).

Input (piped from Tab 2 or uploaded as zip):
  - t2_plates: [(plate_label, xml_bytes), ...]
  - t2_saw:    {sample_name: "Plate1_A1"}  (or legacy {sample_name: "A1"})

Steps per plate:
  1. Sort samples alphabetically, assign new wells A1..H12
  2. Sort XML shapes by new well, update CapID, inject TransferID
  3. Generate 96-well plate CSV + plate map PNG
  4. Generate sample list CSV (cut order, Plate, ROI, Well_ID, Dropout {Y/N}, comments, processed)

Downloads: zip of all per-plate outputs + combined sample list
Pipes to Tab 4: combined sample_list CSV bytes via session_state.t3_sample_list / t3_stem
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


def _split_saw_by_plate(saw_dict: dict) -> dict:
    """
    Split {name: 'Plate1_A1'} into {'Plate1': {name: 'A1'}, 'Plate2': {name: 'A1'}, ...}.
    Falls back to single plate if values have no Plate prefix.
    """
    is_multi = any(
        re.match(r'^Plate\d+_', v) for v in saw_dict.values()
    )
    if not is_multi:
        return {"Plate1": saw_dict}
    result = {}
    for name, full_well in saw_dict.items():
        m = re.match(r'^(Plate\d+)_(.+)$', full_well)
        if m:
            plate_label, well = m.group(1), m.group(2)
            result.setdefault(plate_label, {})[name] = well
        else:
            result.setdefault("Plate1", {})[name] = full_well
    return result


# ============================================================
# CORE LOGIC (single plate)
# ============================================================
def process_collection(xml_bytes: bytes, saw_dict: dict, stem: str) -> dict:
    """
    All in-memory. saw_dict: {sample_name: well_id} with plain wells (no Plate prefix).
    Returns dict with output bytes, grid, metadata.
    """
    warnings = []

    orig_well_to_sample = {v: k for k, v in saw_dict.items()}
    all_samples_alpha   = sorted(saw_dict.keys(), key=str.casefold)

    existing_wells = sorted(saw_dict.values(), key=well_sort_key)
    start_well     = existing_wells[0] if existing_wells else "A1"
    start_idx      = WELLS_96.index(start_well) if start_well in WELLS_96 else 0

    if start_idx > 0:
        warnings.append(f"Custom start position: wells assigned from {start_well}.")

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
    tree      = ET.parse(io.StringIO(xml_bytes.decode("utf-8-sig")))
    root      = tree.getroot()
    shape_pat = re.compile(r'^Shape_\d+$')
    shapes    = [el for el in root if shape_pat.match(el.tag)]
    non_shapes = [el for el in root if not shape_pat.match(el.tag)]

    def shape_key(el):
        cap    = el.find("CapID")
        orig   = cap.text.strip() if cap is not None and cap.text else ""
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

        existing_tid = el.find("TransferID")
        if existing_tid is not None:
            existing_tid.text = sample_name
        else:
            children    = list(el)
            cap_pos     = children.index(cap_el) if cap_el is not None else 0
            transfer_el = ET.Element("TransferID")
            transfer_el.text = sample_name
            el.insert(cap_pos + 1, transfer_el)
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

    # Sample list rows (no plate column here; caller adds it)
    sample_rows = []
    for roi_num, el in enumerate(shapes_sorted, start=1):
        t_el   = el.find("TransferID")
        c_el   = el.find("CapID")
        sample = t_el.text if t_el is not None and t_el.text else ""
        well   = c_el.text if c_el is not None and c_el.text else ""
        sample_rows.append([roi_num, sample, well, "", "", ""])

    updated_json = json.dumps(sample_to_new_well, indent=4, ensure_ascii=False).encode("utf-8")

    return {
        "sorted_xml":    xml_buf.getvalue(),
        "wellplate_csv": wellplate_buf.getvalue().encode("utf-8"),
        "sample_rows":   sample_rows,
        "updated_json":  updated_json,
        "grid":          grid,
        "n_rois":        len(shapes_sorted),
        "n_samples":     len(all_samples_alpha),
        "warnings":      warnings,
        "stem":          stem,
    }


def process_all_plates(plates: list, saw_dict: dict, stem: str) -> dict:
    """
    plates: [(plate_label, xml_bytes), ...]
    saw_dict: {name: 'Plate1_A1'} or {name: 'A1'} (single-plate legacy)
    Returns combined results dict.
    """
    per_plate_saw = _split_saw_by_plate(saw_dict)
    all_results   = {}

    for plate_label, xml_bytes in plates:
        plate_saw = per_plate_saw.get(plate_label, {})
        if not plate_saw:
            continue
        r = process_collection(xml_bytes, plate_saw, f"{stem}_{plate_label}")
        all_results[plate_label] = r

    return all_results


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


def build_combined_sample_list(all_results: dict) -> bytes:
    """Combine per-plate sample rows into one CSV with a Plate column."""
    buf = io.StringIO()
    w   = csv.writer(buf)
    w.writerow(["Plate", "cut order", "ROI", "Well_ID", "Dropout {Y/N}", "comments", "processed"])
    for plate_label in sorted(all_results.keys()):
        for row in all_results[plate_label]["sample_rows"]:
            w.writerow([plate_label] + row)
    return buf.getvalue().encode("utf-8")


def build_multi_zip(all_results: dict, png_map: dict, stem: str) -> bytes:
    """Package all per-plate outputs into one zip."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for plate_label, r in all_results.items():
            ps = f"{stem}_{plate_label}"
            z.writestr(f"{ps}_sorted.xml",      r["sorted_xml"])
            z.writestr(f"{ps}_96wellplate.csv", r["wellplate_csv"])
            z.writestr(f"{ps}_updated_saw.json", r["updated_json"])
            if plate_label in png_map:
                z.writestr(f"{ps}_platemap.png", png_map[plate_label])
        combined_csv = build_combined_sample_list(all_results)
        z.writestr(f"{stem}_sample_list.csv", combined_csv)
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_process_tab():
    st.header("Process LMD Collection")
    st.caption(
        "Sorts ROIs alphabetically per plate, assigns wells A1..H12, updates XML CapID, "
        "generates plate maps and sample list. "
        "Accepts output piped from Tab 2 or a zip upload."
    )

    # Source: piped from Tab 2 or zip upload
    t2_plates = st.session_state.get("t2_plates")   # [(plate_label, xml_bytes)]
    t2_saw    = st.session_state.get("t2_saw")       # {name: "Plate1_A1"}
    t2_stem   = st.session_state.get("t2_stem")

    source = None
    uploaded_zip = st.file_uploader(
        "Upload collection zip (from Tab 2 download)",
        type=["zip"], key="proc_upload"
    )

    if uploaded_zip is not None:
        source = "zip"
    elif t2_plates is not None and t2_saw is not None:
        st.info(f"Using {len(t2_plates)} plate(s) piped from Tab 2.")
        source = "pipe"

    if source is None:
        st.info("Convert GeoJSON in Tab 2 first, or upload a zip here.")
        st.stop()

    # Extract plates + saw
    if source == "pipe":
        plates   = t2_plates
        saw_dict = t2_saw
        stem     = t2_stem or "collection"
    else:
        with zipfile.ZipFile(uploaded_zip) as z:
            names      = z.namelist()
            xml_names  = sorted([
                n for n in names if n.endswith(".xml")
                and "_sorted" not in n and not n.startswith("__MACOSX")
            ])
            json_names = [n for n in names if n.endswith("samples_and_wells.json")]

        if not xml_names:
            st.error("No XML file(s) found in zip.")
            return
        if not json_names:
            st.error("samples_and_wells.json not found in zip.")
            return

        with zipfile.ZipFile(uploaded_zip) as z:
            saw_dict = json.loads(z.read(json_names[0]).decode("utf-8"))
            plates   = []
            for xn in xml_names:
                # Try to extract plate label from filename e.g. stem_Plate1.xml
                m = re.search(r'(Plate\d+)\.xml$', xn, re.IGNORECASE)
                label = m.group(1) if m else xn.replace(".xml", "").split("/")[-1]
                plates.append((label, z.read(xn)))

        stem = xml_names[0].split("/")[-1].replace(".xml", "").rsplit("_", 1)[0]

    # Cache key
    input_key = (stem, len(plates), len(saw_dict))
    if st.session_state.get("proc_last_key") != input_key:
        st.session_state.proc_all_results = None
        st.session_state.proc_last_key    = input_key

    c1, c2 = st.columns(2)
    c1.metric("Samples in JSON", len(saw_dict))
    c2.metric("Plates", len(plates))

    if st.button("Process all plates", type="primary"):
        with st.spinner(f"Processing {len(plates)} plate(s)..."):
            all_results = process_all_plates(plates, saw_dict, stem)
            png_map     = {pl: plot_plate_png(r["grid"], f"{pl} — well assignment")
                           for pl, r in all_results.items()}
            combined_csv = build_combined_sample_list(all_results)
            zip_bytes    = build_multi_zip(all_results, png_map, stem)

        st.session_state.proc_all_results = all_results
        st.session_state.proc_png_map     = png_map
        st.session_state.proc_combined    = combined_csv
        st.session_state.proc_zip         = zip_bytes
        st.session_state.proc_last_key    = input_key
        st.session_state.t3_sample_list   = combined_csv
        st.session_state.t3_stem          = stem

        total_rois = sum(r["n_rois"] for r in all_results.values())
        st.success(f"{total_rois} ROIs processed across {len(all_results)} plate(s).")

    all_results = st.session_state.get("proc_all_results")
    if not all_results:
        return

    # Warnings
    for pl, r in all_results.items():
        for w in r["warnings"]:
            st.warning(f"{pl}: {w}")

    # Plate selector
    plate_labels = sorted(all_results.keys())
    sel_plate    = st.selectbox("View plate", plate_labels, key="proc_plate_sel")
    result       = all_results[sel_plate]
    png_map      = st.session_state.get("proc_png_map", {})

    st.subheader(f"{sel_plate} — 96-well layout")
    if sel_plate in png_map:
        st.image(png_map[sel_plate], use_container_width=True)

    st.subheader("Downloads")
    stem_out = stem

    dl0, dl1, dl2, dl3 = st.columns(4)
    dl0.download_button("All plates (zip)", st.session_state.proc_zip,
                        file_name=f"{stem_out}_lmd_outputs.zip",
                        mime="application/zip", type="primary")
    dl1.download_button(f"{sel_plate} sorted XML", result["sorted_xml"],
                        file_name=f"{stem_out}_{sel_plate}_sorted.xml", mime="application/xml")
    dl2.download_button(f"{sel_plate} 96-well CSV", result["wellplate_csv"],
                        file_name=f"{stem_out}_{sel_plate}_96wellplate.csv", mime="text/csv")
    dl3.download_button("Combined sample list", st.session_state.proc_combined,
                        file_name=f"{stem_out}_sample_list.csv", mime="text/csv")

    # --------------------------------------------------------
    # DROPOUT EDITOR
    # --------------------------------------------------------
    st.divider()
    st.subheader("Stereomicroscope QC — Mark Dropouts")
    st.caption(
        "After inspecting plates under the stereomicroscope, "
        "tick Dropout for any failed wells. "
        "Dropouts are excluded from the MS queue. "
        "Alternatively, upload an updated sample list CSV."
    )

    raw_csv = st.session_state.get("t3_sample_list") or st.session_state.get("proc_combined")
    df      = pd.read_csv(io.BytesIO(raw_csv))
    df.columns = [c.strip() for c in df.columns]

    dropout_col = "Dropout {Y/N}"
    if dropout_col not in df.columns:
        df[dropout_col] = False
    else:
        df[dropout_col] = df[dropout_col].apply(
            lambda v: True if str(v).strip().upper() == "Y" else False
        )

    # Filter to selected plate for display
    if "Plate" in df.columns:
        df_view = df[df["Plate"] == sel_plate].copy()
    else:
        df_view = df.copy()

    edited_df_view = st.data_editor(
        df_view,
        column_config={dropout_col: st.column_config.CheckboxColumn("Dropout?", default=False)},
        disabled=[c for c in df_view.columns if c != dropout_col],
        hide_index=True,
        use_container_width=True,
        key=f"proc_dropout_editor_{sel_plate}",
    )

    updated_upload = st.file_uploader(
        "Or upload updated sample list CSV", type=["csv"], key="proc_updated_csv"
    )

    apply_col, dl_col = st.columns([1, 1])

    if apply_col.button("Apply dropout changes", type="primary", key="proc_apply_dropout"):
        if updated_upload is not None:
            updated_csv = updated_upload.getvalue()
        else:
            # Merge edited plate rows back into full df
            if "Plate" in df.columns:
                df.loc[df["Plate"] == sel_plate, dropout_col] = edited_df_view[dropout_col].values
            else:
                df[dropout_col] = edited_df_view[dropout_col].values
            df[dropout_col] = df[dropout_col].apply(lambda v: "Y" if v else "")
            buf = io.StringIO()
            df.to_csv(buf, index=False)
            updated_csv = buf.getvalue().encode("utf-8")

        st.session_state.t3_sample_list = updated_csv
        n_drop = pd.read_csv(io.BytesIO(updated_csv))[dropout_col].apply(
            lambda v: str(v).strip().upper() == "Y"
        ).sum()
        st.success(f"Dropout list updated — {int(n_drop)} dropout(s). Tab 4 will exclude these.")

    current_csv = st.session_state.get("t3_sample_list") or raw_csv
    dl_col.download_button(
        "Download updated sample list", current_csv,
        file_name=f"{stem_out}_sample_list.csv", mime="text/csv",
    )

    st.caption("Tab 4 (MS Queue) uses this sample list — apply dropout changes before generating the queue.")
