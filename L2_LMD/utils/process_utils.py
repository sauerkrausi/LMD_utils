"""
process_utils.py
================
Tab 3 -- Post-cutting QC (dropout editor + re-run XML generation).

Input (piped from Tab 2 or CSV upload):
  - t2_plates: [(plate_label, xml_bytes), ...]  -- original XMLs
  - t2_saw:    {sample_name: "Plate1_A1"}       -- well assignments
  - t2_stem:   str                              -- file stem
  - t2_groups: {roi_name: group_label}          -- group assignments

Outputs:
  - Updated sample list CSV (with Dropout Y/N filled in)
  - Dropout-free XMLs zip (original XMLs with dropout ROIs removed, for LMD re-run)
  - Pipes t3_sample_list to Tab 4 (MS queue excludes dropouts automatically)
"""

import io
import csv
import json
import re
import zipfile
import xml.etree.ElementTree as ET

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


def _split_saw_by_plate(saw_dict: dict) -> dict:
    """Split {name: 'Plate1_A1'} into {'Plate1': {name: 'A1'}, ...}."""
    is_multi = any(re.match(r'^Plate\d+_', v) for v in saw_dict.values())
    if not is_multi:
        return {"Plate1": saw_dict}
    result = {}
    for name, full_well in saw_dict.items():
        m = re.match(r'^(Plate\d+)_(.+)$', full_well)
        if m:
            result.setdefault(m.group(1), {})[name] = m.group(2)
        else:
            result.setdefault("Plate1", {})[name] = full_well
    return result


def build_sample_list_from_t2(saw_dict: dict, groups_dict: dict, stem: str) -> bytes:
    """Build sample list CSV directly from Tab 2 well assignments."""
    per_plate = _split_saw_by_plate(saw_dict)
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(["Plate", "cut order", "ROI", "Well_ID", "Group",
                "Dropout {Y/N}", "comments", "processed"])
    for plate in sorted(per_plate.keys()):
        pairs = sorted(per_plate[plate].items(), key=lambda x: well_sort_key(x[1]))
        for i, (roi, well) in enumerate(pairs, start=1):
            grp = (groups_dict or {}).get(roi, roi.split("-")[0] if "-" in roi else roi)
            w.writerow([plate, i, roi, well, grp, "", "", ""])
    return buf.getvalue().encode("utf-8")


def build_dropout_free_xmls(plates: list, dropout_rois: set) -> list:
    """Remove dropout ROIs from each plate XML. Returns [(label, xml_bytes), ...]."""
    results = []
    shape_pat = re.compile(r'^Shape_\d+$')

    for plate_label, xml_bytes in plates:
        tree = ET.parse(io.StringIO(xml_bytes.decode("utf-8-sig")))
        root = tree.getroot()

        shapes     = [el for el in root if shape_pat.match(el.tag)]
        non_shapes = [el for el in root if not shape_pat.match(el.tag)]

        # Keep non-dropout shapes only
        kept = []
        for el in shapes:
            t_el = el.find("TransferID")
            roi  = (t_el.text or "").strip() if t_el is not None else ""
            if roi not in dropout_rois:
                kept.append(el)

        # Rebuild root
        for el in list(root):
            root.remove(el)
        for el in non_shapes:
            root.append(el)

        sc = root.find("ShapeCount")
        if sc is not None:
            sc.text = str(len(kept))

        for i, el in enumerate(kept, start=1):
            el.tag = f"Shape_{i}"
            root.append(el)

        buf = io.BytesIO()
        tree.write(buf, encoding="UTF-8", xml_declaration=True)
        results.append((plate_label, buf.getvalue()))

    return results


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_process_tab():
    st.header("Post-Cutting QC")
    st.caption(
        "Mark dropout ROIs after stereomicroscope inspection. "
        "Downloads updated sample list and dropout-free XMLs for LMD re-runs on fresh sections."
    )

    t2_plates = st.session_state.get("t2_plates")
    t2_saw    = st.session_state.get("t2_saw")
    t2_stem   = st.session_state.get("t2_stem", "collection")
    t2_groups = st.session_state.get("t2_groups", {})

    # Backward compat: old single-XML key
    if t2_plates is None and st.session_state.get("t2_xml") is not None:
        t2_plates = [("Plate1", st.session_state.t2_xml)]

    uploaded_csv = st.file_uploader(
        "Upload sample list CSV (optional — overrides Tab 2 pipe)",
        type=["csv"], key="proc_upload_csv"
    )

    # Resolve sample list source
    if uploaded_csv is not None:
        raw_csv = uploaded_csv.getvalue()
        stem    = uploaded_csv.name.replace("_sample_list.csv", "").replace(".csv", "")
        st.session_state.t3_sample_list = raw_csv
        st.session_state.t3_stem        = stem
    elif st.session_state.get("t3_sample_list") is not None:
        raw_csv = st.session_state.t3_sample_list
        stem    = st.session_state.get("t3_stem", t2_stem)
        st.info("Using sample list piped from Tab 2.")
    elif t2_saw is not None:
        raw_csv = build_sample_list_from_t2(t2_saw, t2_groups, t2_stem)
        stem    = t2_stem
        st.session_state.t3_sample_list = raw_csv
        st.session_state.t3_stem        = stem
        st.info("Sample list built from Tab 2 well assignments.")
    else:
        st.info("Complete Tab 2 first, or upload a sample list CSV here.")
        st.stop()

    # Parse
    df = pd.read_csv(io.BytesIO(raw_csv))
    df.columns     = [c.strip() for c in df.columns]
    dropout_col    = "Dropout {Y/N}"
    if dropout_col not in df.columns:
        df[dropout_col] = False
    else:
        df[dropout_col] = df[dropout_col].apply(
            lambda v: True if str(v).strip().upper() == "Y" else False
        )

    n_rois = len(df)
    n_plates = df["Plate"].nunique() if "Plate" in df.columns else 1

    st.subheader("Dropout Editor")
    st.caption(f"{n_rois} ROIs across {n_plates} plate(s). Check **Dropout?** for any failed wells.")

    edited_df = st.data_editor(
        df,
        column_config={dropout_col: st.column_config.CheckboxColumn("Dropout?", default=False)},
        disabled=[c for c in df.columns if c != dropout_col],
        hide_index=True,
        use_container_width=True,
        key="proc_dropout_editor",
    )

    n_dropouts = int(edited_df[dropout_col].sum())
    if n_dropouts:
        st.warning(f"{n_dropouts} dropout(s) marked — these will be excluded from the MS queue and re-run XMLs.")
    else:
        st.success("No dropouts marked.")

    st.divider()

    # Build updated CSV (Y/N strings for export)
    export_df              = edited_df.copy()
    export_df[dropout_col] = export_df[dropout_col].apply(lambda v: "Y" if v else "")
    out_buf                = io.StringIO()
    export_df.to_csv(out_buf, index=False)
    updated_csv = out_buf.getvalue().encode("utf-8")

    # Always pipe latest to Tab 4
    st.session_state.t3_sample_list = updated_csv
    st.session_state.t3_stem        = stem

    dl1, dl2 = st.columns(2)
    dl1.download_button(
        "Download sample list (CSV)", updated_csv,
        file_name=f"{stem}_sample_list.csv", mime="text/csv",
    )

    # Dropout-free re-run XMLs
    if t2_plates is not None and n_dropouts > 0:
        dropout_rois = set(
            edited_df.loc[edited_df[dropout_col].astype(bool), "ROI"].astype(str)
        )
        rerun_xmls = build_dropout_free_xmls(t2_plates, dropout_rois)

        rerun_buf = io.BytesIO()
        with zipfile.ZipFile(rerun_buf, "w", zipfile.ZIP_DEFLATED) as z:
            for plate_label, xml_bytes in rerun_xmls:
                z.writestr(f"{stem}_{plate_label}_rerun.xml", xml_bytes)
        dl2.download_button(
            "Download re-run XMLs (zip)", rerun_buf.getvalue(),
            file_name=f"{stem}_rerun.zip", mime="application/zip",
            type="primary",
        )
    elif n_dropouts == 0:
        dl2.caption("No dropouts — no re-run XML needed.")
    else:
        dl2.caption("Complete Tab 2 to enable re-run XML generation.")

    st.caption("Tab 4 (MS Queue) reads this sample list and automatically excludes dropout ROIs.")
