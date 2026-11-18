"""
L2 msQ
======
Standalone MS queue generator for the timsTOF Ultra (HyStar XLSX) and the
Astral Zoom with Evosep One or Vanquish Neo (Xcalibur CSV).

Runs the same queue code as Tab 4 of L2_LMD (L2_LMD/utils/ms_queue_utils.py),
so both apps stay identical. Bulk and SCP run modes, 96 or 384 well plates on
the Ultra, CellenONE input supported.

Run:
    cd L2_msQ
    pip install -r requirements.txt
    streamlit run app.py
"""

import sys
from pathlib import Path

import streamlit as st

st.set_page_config(page_title="L2 msQ", layout="wide")

# Shared queue code lives in the L2_LMD app of the same repo
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "L2_LMD"))
from utils.ms_queue_utils import render_ms_queue_tab  # noqa: E402

st.title("L2 msQ")
st.caption("Standalone MS queue generator, same queue logic as Tab 4 of L2-LMD.")
render_ms_queue_tab(standalone=True)
