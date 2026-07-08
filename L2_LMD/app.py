"""
L2-LMD — Lange Lab LMD Platform
=================================
Integrated Streamlit app for the full LMD spatial proteomics workflow.

Tabs:
  1. Reclassify  — fix QuPath GeoJSON classification names
  2. Convert     — GeoJSON -> LMD XML via py-lmd (MannLabs)
  3. Process     — sort XML, generate plate map + sample list
  4. MS Queue    — generate Bruker timsTOF instrument queue

Run:
  streamlit run app.py

Attribution:
  XML generation: py-lmd (MannLabs) — https://github.com/MannLabs/py-lmd — Apache-2.0
  Workflow reference: Qupath_to_LMD (Coscia Lab) — https://github.com/CosciaLab/Qupath_to_LMD — GPL-3.0
"""

import streamlit as st

st.set_page_config(page_title="L2-LMD", layout="wide", page_icon="🔬")

st.markdown(
    "<div style='position:fixed;bottom:10px;right:16px;font-size:11px;color:#aaa;'>"
    "py-lmd: <a href='https://github.com/MannLabs/py-lmd' target='_blank'>MannLabs</a> | "
    "Workflow ref: <a href='https://github.com/CosciaLab/Qupath_to_LMD' target='_blank'>Coscia Lab</a> | "
    "Source: <a href='https://github.com/sauerkrausi/LMD_utils' target='_blank'>LMD_utils</a>"
    "</div>",
    unsafe_allow_html=True,
)

# ============================================================
# SESSION STATE — cross-tab data transfer
# ============================================================
# Tab 1 -> Tab 2: reclassified GeoJSON bytes
if "t1_geojson" not in st.session_state:
    st.session_state.t1_geojson = None      # bytes
if "t1_stem" not in st.session_state:
    st.session_state.t1_stem = None         # str, filename stem
if "geo_result" not in st.session_state:
    st.session_state.geo_result = None
if "geo_last" not in st.session_state:
    st.session_state.geo_last = None

# Tab 2 internal
if "conv_last" not in st.session_state:
    st.session_state.conv_last = None

# Tab 2 -> Tab 3: XML + samples_and_wells.json bytes
if "t2_xml" not in st.session_state:
    st.session_state.t2_xml = None          # bytes
if "t2_saw" not in st.session_state:
    st.session_state.t2_saw = None          # dict {sample_name: well_id}
if "t2_stem" not in st.session_state:
    st.session_state.t2_stem = None         # str

# Tab 3 internal
if "proc_result" not in st.session_state:
    st.session_state.proc_result = None
if "proc_png" not in st.session_state:
    st.session_state.proc_png = None
if "proc_zip" not in st.session_state:
    st.session_state.proc_zip = None
if "proc_last_key" not in st.session_state:
    st.session_state.proc_last_key = None

# Tab 3 -> Tab 4: sample list CSV bytes
if "t3_sample_list" not in st.session_state:
    st.session_state.t3_sample_list = None  # bytes
if "t3_stem" not in st.session_state:
    st.session_state.t3_stem = None         # str

# Tab 4 internal
if "msq_results" not in st.session_state:
    st.session_state.msq_results = None
if "msq_zip" not in st.session_state:
    st.session_state.msq_zip = None
if "msq_last" not in st.session_state:
    st.session_state.msq_last = None

# ============================================================
# TABS
# ============================================================
tab1, tab2, tab3, tab4 = st.tabs([
    "1. Reclassify",
    "2. Convert",
    "3. Process",
    "4. MS Queue",
])

# ============================================================
# TAB 1 — Reclassify GeoJSON
# ============================================================
with tab1:
    from utils.geojson_utils import render_reclassify_tab
    render_reclassify_tab()

# ============================================================
# TAB 2 — Convert GeoJSON -> LMD XML
# ============================================================
with tab2:
    from utils.convert_utils import render_convert_tab
    render_convert_tab()

# ============================================================
# TAB 3 — Process LMD Collection
# ============================================================
with tab3:
    from utils.process_utils import render_process_tab
    render_process_tab()

# ============================================================
# TAB 4 — MS Sample Queue
# ============================================================
with tab4:
    from utils.ms_queue_utils import render_ms_queue_tab
    render_ms_queue_tab()
