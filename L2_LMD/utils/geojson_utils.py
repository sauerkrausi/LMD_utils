"""
geojson_utils.py
================
Tab 1 — Reclassify GeoJSON.

Copies properties.name -> properties.classification.name for each feature.
Required when QuPath annotations were drawn without a classification assigned.
"""

import json
import streamlit as st


# ============================================================
# CORE LOGIC
# ============================================================
def reclassify_geojson(geojson_bytes: bytes) -> tuple:
    """
    For each feature: copy properties.name into properties.classification.name.
    Returns (corrected_bytes, n_fixed, names_list, rejected_list).
    """
    geojson  = json.loads(geojson_bytes.decode("utf-8"))
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

    corrected = json.dumps(geojson, indent=2).encode("utf-8")
    return corrected, n_fixed, names, rejected


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_reclassify_tab():
    st.header("Reclassify GeoJSON")
    st.caption(
        "Copies `properties.name` into `properties.classification.name` for each annotation. "
        "Required before converting if QuPath annotations lack a classification."
    )

    # Accept piped output from a previous session or fresh upload
    pipe = st.session_state.get("t1_geojson")

    uploaded = st.file_uploader(
        "Upload QuPath GeoJSON export", type=["geojson", "json"], key="geo_upload"
    )

    # Determine active source
    if uploaded:
        raw   = uploaded.read()
        stem  = uploaded.name.replace(".geojson", "").replace(".json", "")
        # clear any previous result when a new file arrives
        if st.session_state.get("geo_last") != uploaded.name:
            st.session_state.geo_result = None
            st.session_state.geo_last   = uploaded.name
    elif pipe is not None:
        raw  = pipe
        stem = st.session_state.get("t1_stem", "geojson")
        st.info("Using GeoJSON piped from a previous step.")
    else:
        raw  = None
        stem = None

    if raw is None:
        st.stop()

    # Parse and preview
    try:
        geojson     = json.loads(raw.decode("utf-8"))
        features    = geojson.get("features", [])
        all_names   = [f.get("properties", {}).get("name", "") for f in features]
        n_total     = len(features)
        n_named     = sum(1 for n in all_names if n.strip())
        n_points    = sum(1 for f in features
                         if f.get("geometry", {}).get("type") == "Point")
        n_polygons  = n_total - n_points
    except Exception as e:
        st.error(f"Could not parse GeoJSON: {e}")
        st.stop()

    col1, col2, col3 = st.columns(3)
    col1.metric("Total features", n_total)
    col2.metric("Polygons / annotations", n_polygons)
    col3.metric("Named features", n_named)

    if n_named < n_total:
        st.warning(f"{n_total - n_named} features have no name and will be skipped.")

    with st.expander("Preview annotation names", expanded=False):
        named = [n for n in all_names if n.strip()]
        st.write(sorted(set(named)))

    st.divider()

    if st.button("Reclassify", type="primary"):
        corrected, n_fixed, names, rejected = reclassify_geojson(raw)
        st.session_state.geo_result = {
            "bytes":    corrected,
            "stem":     stem,
            "n_fixed":  n_fixed,
            "names":    names,
            "rejected": rejected,
        }

    r = st.session_state.get("geo_result")
    if r:
        out_name = f"{r['stem']}_reclassified.geojson"

        if r["n_fixed"] > 0:
            st.success(f"{r['n_fixed']} features updated.")
        else:
            st.info("All features already had correct classification — file passed through unchanged.")

        if r["rejected"]:
            st.warning(f"{len(r['rejected'])} unnamed features skipped.")

        # Pass to Tab 2 via session state
        st.session_state.t1_geojson = r["bytes"]
        st.session_state.t1_stem    = r["stem"]

        st.download_button(
            label    = f"Download {out_name}",
            data     = r["bytes"],
            file_name= out_name,
            mime     = "application/json",
            type     = "primary",
        )
        st.caption("Output is available in Tab 2 (Convert) without re-uploading.")
