"""
L2 msQ
======
Standalone MS queue generator for the timsTOF Ultra (HyStar XLSX) and the
Astral Zoom with Evosep One or Vanquish Neo (Xcalibur CSV).

Queue code in ms_queue_utils.py next to this file, its own copy, independent
of L2_LMD. Bulk and SCP run modes, 96 or 384 well plates on the Ultra,
CellenONE input supported.

Run:
    cd L2_msQ
    pip install -r requirements.txt
    streamlit run app.py
"""

import sys
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="L2 msQ", layout="wide")

# Queue code lives next to this file
sys.path.insert(0, str(Path(__file__).resolve().parent))
from ms_queue_utils import render_ms_queue_tab  # noqa: E402

st.title("L2 msQ")
st.caption("Standalone MS queue generator for the timsTOF Ultra and the Astral Zoom.")
render_ms_queue_tab(standalone=True)
