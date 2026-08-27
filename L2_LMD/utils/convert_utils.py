"""
convert_utils.py
================
Tab 2 -- Convert reclassified GeoJSON to LMD XML via py-lmd (MannLabs).

Workflow:
  1. Parse GeoJSON: extract Point features (calibration candidates) + Polygon ROIs
  2. User selects 3 calibration points (from Points or manual entry)
     -> shows how many ROI centroids fall within the calibration triangle
  3. Wells assigned alphabetically across as many 96-well plates as needed
     -> plate selector dropdown for per-plate preview
  4. One Collection built per plate, saved as XML
  5. Downloads: zip of all XMLs + samples_and_wells.json {name: "Plate1_A1"}
  6. Piped to Tab 3 via session_state.t2_plates / t2_saw / t2_stem / t2_zip
"""

import io
import json
import os
import random
import tempfile
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm
import numpy as np
import streamlit as st

try:
    from lmd.lib import Collection
    HAS_LMD = True
except ImportError:
    HAS_LMD = False

try:
    from shapely.geometry import Point, Polygon as ShapelyPolygon
    HAS_SHAPELY = True
except ImportError:
    HAS_SHAPELY = False


ROWS      = list("ABCDEFGH")
COLS      = list(range(1, 13))
ALL_WELLS = [f"{r}{c}" for r in ROWS for c in COLS]   # A1 .. H12


# ============================================================
# CORE LOGIC
# ============================================================
def parse_geojson(geojson_bytes: bytes):
    gj       = json.loads(geojson_bytes.decode("utf-8"))
    points   = []
    polygons = []
    for f in gj.get("features", []):
        geom  = f.get("geometry", {})
        gtype = geom.get("type", "")
        name  = f.get("properties", {}).get("name", "").strip()
        if gtype == "Point":
            points.append({"name": name, "coords": geom["coordinates"], "feature": f})
        elif gtype in ("Polygon", "MultiPolygon"):
            polygons.append({"name": name, "geom": geom, "feature": f})
    return points, polygons


def polygon_exterior(geom: dict) -> np.ndarray:
    if geom["type"] == "Polygon":
        ring = geom["coordinates"][0]
    else:
        ring = geom["coordinates"][0][0]
    return np.array([[c[0], c[1]] for c in ring], dtype=float)


def polygon_centroid(geom: dict) -> np.ndarray:
    return polygon_exterior(geom).mean(axis=0)


def count_inside_triangle(calib_pts: np.ndarray, polygons: list):
    """Count how many polygon centroids lie inside the calibration triangle."""
    if not HAS_SHAPELY or calib_pts is None:
        return None
    tri = ShapelyPolygon(calib_pts)
    if not tri.is_valid or tri.area == 0:
        return None
    return sum(
        1 for p in polygons
        if tri.contains(Point(polygon_centroid(p["geom"])))
    )


def assign_wells(polygons: list, randomize: bool = False, seed: int = 42) -> dict:
    """Returns {name: 'Plate1_A1'} across as many 96-well plates as needed."""
    names = sorted(set(p["name"] for p in polygons))
    if randomize:
        random.Random(seed).shuffle(names)
    result = {}
    for i, name in enumerate(names):
        plate_num    = i // 96 + 1
        well         = ALL_WELLS[i % 96]
        result[name] = f"Plate{plate_num}_{well}"
    return result


def get_plate_labels(well_map: dict) -> list:
    """Sorted list of unique plate labels e.g. ['Plate1', 'Plate2']."""
    return sorted({v.split("_")[0] for v in well_map.values()})


def filter_plate_well_map(well_map: dict, plate_label: str) -> dict:
    """Returns {name: 'A1'} for one plate, stripping the Plate prefix."""
    prefix = f"{plate_label}_"
    return {name: w[len(prefix):] for name, w in well_map.items() if w.startswith(prefix)}


def _build_single_xml(calib_pts: np.ndarray, polygons: list, pwm: dict) -> bytes:
    """Build one XML for one plate. pwm: {name: 'A1'}"""
    import xml.etree.ElementTree as ET

    well_to_name = {v: k for k, v in pwm.items()}
    col          = Collection(calibration_points=calib_pts)
    for p in polygons:
        coords = polygon_exterior(p["geom"])
        col.new_shape(coords, well=pwm.get(p["name"], "A1"), name=p["name"])

    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as f:
        tmppath = f.name
    try:
        col.save(tmppath)
        tree = ET.parse(tmppath)
        root = tree.getroot()
        for el in root:
            if not el.tag.startswith("Shape_"):
                continue
            cap_el = el.find("CapID")
            if cap_el is None:
                continue
            sample_name = well_to_name.get(cap_el.text.strip() if cap_el.text else "", "")
            if el.find("TransferID") is None and sample_name:
                children = list(el)
                cap_pos  = children.index(cap_el)
                tid      = ET.Element("TransferID")
                tid.text = sample_name
                el.insert(cap_pos + 1, tid)
        tree.write(tmppath, encoding="UTF-8", xml_declaration=True)
        with open(tmppath, "rb") as f:
            return f.read()
    finally:
        os.unlink(tmppath)


def build_xml_plates(calib_pts: np.ndarray, polygons: list, well_map: dict) -> list:
    """Returns [(plate_label, xml_bytes), ...] one entry per plate."""
    results = []
    for plate_label in get_plate_labels(well_map):
        pwm         = filter_plate_well_map(well_map, plate_label)
        plate_polys = [p for p in polygons if p["name"] in pwm]
        xml_bytes   = _build_single_xml(calib_pts, plate_polys, pwm)
        results.append((plate_label, xml_bytes))
    return results


def build_xml_zip(plate_xml_list: list, well_map: dict, stem: str) -> bytes:
    """Zip all plate XMLs + samples_and_wells.json."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for plate_label, xml_bytes in plate_xml_list:
            z.writestr(f"{stem}_{plate_label}.xml", xml_bytes)
        z.writestr("samples_and_wells.json",
                   json.dumps(well_map, indent=2))
    return buf.getvalue()


# ============================================================
# 96-WELL PLATE PREVIEW
# ============================================================
def plot_well_preview(well_map: dict, title: str) -> bytes:
    """well_map: {name: 'A1'} for a single plate."""
    groups       = sorted({n.split("_")[0] for n in well_map})
    palette      = cm.tab20
    group_color  = {g: palette(i / max(len(groups), 1)) for i, g in enumerate(groups)}
    well_to_name = {v: k for k, v in well_map.items()}

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(-0.8, 12.5)
    ax.set_ylim(-0.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)

    for r_idx, r in enumerate(ROWS):
        for c_idx, c in enumerate(COLS):
            x    = c_idx
            y    = 7 - r_idx
            well = f"{r}{c}"
            name = well_to_name.get(well, "")
            grp  = name.split("_")[0] if name else ""
            color = group_color.get(grp, "whitesmoke") if name else "whitesmoke"
            edge  = "#999999" if not name else "#333333"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.7, zorder=2))
            if name:
                ax.text(x, y, name.replace("_", "\n"), ha="center", va="center",
                        fontsize=3.5, zorder=3, color="black")

    for r_idx, r in enumerate(ROWS):
        ax.text(-0.65, 7 - r_idx, r, ha="right", va="center", fontsize=8, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=8, fontweight="bold")

    patches = [mpatches.Patch(color=group_color[g], label=g) for g in groups]
    if patches:
        ax.legend(handles=patches, bbox_to_anchor=(1.01, 1), loc="upper left",
                  fontsize=6, title="Group", title_fontsize=7)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_convert_tab():
    st.header("Convert GeoJSON to LMD XML")
    st.caption(
        "Converts reclassified GeoJSON to LMD XML using "
        "[py-lmd](https://github.com/MannLabs/py-lmd) (MannLabs, Apache-2.0). "
        "Select 3 calibration points, then assign wells across one or more 96-well plates."
    )

    if not HAS_LMD:
        st.error("`py-lmd` not installed. Run `pip install py-lmd`.")
        return

    pipe     = st.session_state.get("t1_geojson")
    uploaded = st.file_uploader(
        "Upload reclassified GeoJSON", type=["geojson", "json"], key="conv_upload"
    )

    if uploaded:
        raw  = uploaded.read()
        stem = uploaded.name.replace(".geojson", "").replace(".json", "")
        if st.session_state.get("conv_last") != uploaded.name:
            st.session_state.t2_plates = None
            st.session_state.conv_last = uploaded.name
    elif pipe is not None:
        raw  = pipe
        stem = st.session_state.get("t1_stem", "geojson")
        st.info("Using GeoJSON piped from Tab 1.")
    else:
        st.stop()

    points, polygons = parse_geojson(raw)

    n_unique = len(set(p["name"] for p in polygons))
    n_plates = max(1, (n_unique + 95) // 96)
    c1, c2, c3 = st.columns(3)
    c1.metric("Calibration point candidates", len(points))
    c2.metric("Polygon ROIs", len(polygons))
    c3.metric("Plates needed", n_plates)

    if not polygons:
        st.error("No polygon features found in GeoJSON.")
        return

    dup_names = [n for n in (p["name"] for p in polygons)
                 if sum(1 for p2 in polygons if p2["name"] == n) > 1]
    if dup_names:
        st.warning(f"Duplicate annotation names: {sorted(set(dup_names))}. "
                   "Only the first occurrence of each name gets a unique well.")

    st.divider()

    # Calibration points
    st.subheader("Calibration Points")
    calib_pts = None
    calib_ok  = False

    if len(points) >= 3:
        pt_labels = [
            f"{p['name']}  ({p['coords'][0]:.1f}, {p['coords'][1]:.1f})"
            for p in points
        ]
        st.caption("Select exactly 3 different Point annotations as calibration markers.")
        sel_cols = st.columns(3)
        selected = [
            sel_cols[i].selectbox(
                f"Calibration {i+1}", options=range(len(points)),
                format_func=lambda x, L=pt_labels: L[x],
                key=f"calib_{i}"
            )
            for i in range(3)
        ]
        if len(set(selected)) < 3:
            st.warning("Select 3 different points.")
        else:
            calib_ok  = True
            calib_pts = np.array([points[i]["coords"][:2] for i in selected], dtype=float)
    else:
        st.caption(
            "No Point annotations found. "
            "Enter calibration coordinates manually (pixel x, y from QuPath)."
        )
        calib_pts_raw = []
        for i in range(3):
            cx, cy = st.columns(2)
            x = cx.number_input(f"Calibration {i+1} X", value=0.0, key=f"cx_{i}")
            y = cy.number_input(f"Calibration {i+1} Y", value=0.0, key=f"cy_{i}")
            calib_pts_raw.append([x, y])
        calib_pts = np.array(calib_pts_raw, dtype=float)
        calib_ok  = True

    # Triangle coverage
    if calib_ok and calib_pts is not None:
        n_inside = count_inside_triangle(calib_pts, polygons)
        if n_inside is not None:
            pct = 100 * n_inside / len(polygons) if polygons else 0
            if n_inside < len(polygons):
                st.warning(
                    f"{n_inside}/{len(polygons)} ROI centroids ({pct:.0f}%) inside the "
                    f"calibration triangle. ROIs outside may be cut inaccurately."
                )
            else:
                st.success(f"All {n_inside} ROI centroids lie inside the calibration triangle.")
        elif not HAS_SHAPELY:
            st.caption("Install `shapely` for calibration triangle coverage check.")

    st.divider()

    # Well assignment
    st.subheader("Well Assignment")
    randomize = st.checkbox("Randomize well order", value=False)
    seed      = 42
    if randomize:
        seed = int(st.number_input("Random seed", value=42, step=1, key="rand_seed"))

    well_map     = assign_wells(polygons, randomize=randomize, seed=seed)
    plate_labels = get_plate_labels(well_map)

    # Show all plates in tabs
    plate_tabs = st.tabs(plate_labels)
    for tab, pl in zip(plate_tabs, plate_labels):
        with tab:
            pwm        = filter_plate_well_map(well_map, pl)
            n_in_plate = len(pwm)
            st.caption(f"{pl}: {n_in_plate} ROIs")
            col_tbl, col_plate = st.columns([1, 2])
            with col_tbl:
                with st.expander("Table", expanded=False):
                    st.dataframe(
                        [{"Annotation": k, "Well": v} for k, v in sorted(pwm.items())],
                        use_container_width=True, hide_index=True
                    )
            with col_plate:
                preview_png = plot_well_preview(pwm, f"{pl} — {n_in_plate} ROIs")
                st.image(preview_png, use_container_width=True)

    st.divider()

    if st.button("Convert to XML", type="primary", disabled=not calib_ok):
        with st.spinner(f"Building {len(plate_labels)} plate(s)..."):
            plate_xml_list = build_xml_plates(calib_pts, polygons, well_map)
            zip_bytes      = build_xml_zip(plate_xml_list, well_map, stem)

        st.session_state.t2_plates = plate_xml_list
        st.session_state.t2_saw    = well_map
        st.session_state.t2_stem   = stem
        st.session_state.t2_zip    = zip_bytes
        st.success(f"Converted {len(polygons)} ROIs across {len(plate_labels)} plate(s).")

    if st.session_state.get("t2_plates"):
        stem_out      = st.session_state.get("t2_stem") or stem
        zip_bytes_dl  = st.session_state.t2_zip
        saw_bytes     = json.dumps(st.session_state.t2_saw, indent=2).encode("utf-8")
        plates_out    = st.session_state.t2_plates

        # Download buttons: zip + individual plates + JSON
        n_btns   = len(plates_out) + 2
        dl_cols  = st.columns(n_btns)
        dl_cols[0].download_button(
            "All plates (zip)", zip_bytes_dl,
            file_name=f"{stem_out}_lmd.zip", mime="application/zip", type="primary"
        )
        for i, (plate_label, xml_bytes) in enumerate(plates_out):
            dl_cols[i + 1].download_button(
                f"{plate_label}.xml", xml_bytes,
                file_name=f"{stem_out}_{plate_label}.xml", mime="application/xml"
            )
        dl_cols[-1].download_button(
            "samples_and_wells.json", saw_bytes,
            file_name=f"{stem_out}_samples_and_wells.json", mime="application/json"
        )
        st.caption("Outputs piped to Tab 3 (Process) without re-uploading.")
