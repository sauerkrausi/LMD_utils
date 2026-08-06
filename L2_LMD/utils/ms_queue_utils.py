"""
ms_queue_utils.py
=================
Tab 4 -- MS Sample Queue Generator.

Ported from create_ms_queue_streamlit.py.
Accepts sample list CSV piped from Tab 3 or via direct upload.

Outputs: queue XLSX + slot plate CSVs + plate map PNGs (zip).
"""

import csv
import datetime
import io
import json
import math
import re
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st

# ============================================================
# METHOD LOOKUP
# ============================================================
LC_METHODS = {
    "WhisperZOOM40": (
        r"D:\Methods\LC_Methods\Evosep"
        r"\WhisperZOOM_40_SPD_32p5min.m?HyStar_LC"
    ),
}

MS_METHODS = {
    "diaPASEF": (
        r"D:\Methods\MS_Methods\DIA\Farah\TimsControl methods"
        r"\DIA_PASEF_Var_windows_test4_pydiAID_300to1200_80PASEF_scans_-05shift"
        r".proteoscape.m?OtofImpacTEMControl",
        r"D:\Methods\MS_Methods\DIA\Farah\TimsControl methods"
        r"\DIA_PASEF_Var_windows_test4_pydiAID_300to1200_80PASEF_scans_-05shift"
        r".proteoscape.m?DataAnalysis",
    ),
}

LC_OPTIONS = list(LC_METHODS.keys()) + ["Custom"]
MS_OPTIONS = list(MS_METHODS.keys()) + ["Custom"]

ROWS      = list("ABCDEFGH")
COLS      = list(range(1, 13))
GROUP_SIZE = 6

QUEUE_COLS = [
    "Vial", "Sample ID", "Method Set", "Separation Method",
    "Injection Method", "MS Method", "Processing Method",
    "Sample Type", "Volume [µl]", "Data Path", "Run Automated Processing",
]

CTRL_COLORS       = {"K562": "#4C9BE8", "Supermix": "#F4A261", "Blank": "#B7E4C7"}
CTRL_COLORS_SPARE = {"K562": "#C5DDF7", "Supermix": "#FAE0C8", "Blank": "#E6F7EC"}
GREY_REPLACEMENTS = {14: "#777777", 15: "#555555"}


# ============================================================
# HELPERS
# ============================================================
def well_to_slot1(well_id):
    row = ord(well_id[0].upper()) - ord('A')
    col = int(well_id[1:])
    return f"Slot1:{row * 12 + col}"


def index_to_well(index):
    i = index - 1
    return ROWS[i // 12], (i % 12) + 1


def split_groups(samples, max_size=GROUP_SIZE):
    n = len(samples)
    if n == 0:
        return []
    num_groups = max(1, math.ceil(n / max_size))
    base, extra = divmod(n, num_groups)
    groups, start = [], 0
    for i in range(num_groups):
        size = base + (1 if i < extra else 0)
        groups.append(samples[start:start + size])
        start += size
    return groups


def make_row(vial, sample_id, data_path, sep_method, inj_method, ms_method, proc_method):
    return {
        "Vial": vial,
        "Sample ID": sample_id,
        "Method Set": "",
        "Separation Method": sep_method,
        "Injection Method": inj_method,
        "MS Method": ms_method,
        "Processing Method": proc_method,
        "Sample Type": "Sample",
        "Volume [µl]": 1,
        "Data Path": data_path,
        "Run Automated Processing": "False",
    }


# ============================================================
# PLATE VISUALIZATION
# ============================================================
def plot_plate(grid, color_map, title, label_map=None, legend_group_map=None,
               label_color_map=None) -> bytes:
    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(-0.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    for r_idx, r in enumerate(ROWS):
        for c_idx, c in enumerate(COLS):
            x, y  = c_idx, 7 - r_idx
            label = grid[r][c]
            color = color_map.get(label, "white") if label else "whitesmoke"
            edge  = "#aaaaaa" if not label else "#444444"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.8, zorder=2))
            display    = (label_map or {}).get(label, label)
            font_color = (label_color_map or {}).get(label, "black")
            if display:
                fsize = 4.5 if len(display) > 10 else 5.5
                ax.text(x, y, display, ha="center", va="center",
                        fontsize=fsize, zorder=3, color=font_color, clip_on=True)

    for r_idx, r in enumerate(ROWS):
        ax.text(-0.55, 7 - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")

    seen = {}
    for r in ROWS:
        for c in COLS:
            lbl = grid[r][c]
            if lbl and lbl not in seen:
                seen[lbl] = color_map.get(lbl, "white")

    if legend_group_map:
        group_seen = {}
        for lbl in sorted(seen):
            grp = legend_group_map.get(lbl, lbl)
            if grp not in group_seen:
                group_seen[grp] = seen[lbl]
        patches = [mpatches.Patch(color=col, label=grp) for grp, col in group_seen.items()]
    else:
        patches = [mpatches.Patch(color=seen[l], label=(label_map or {}).get(l, l))
                   for l in sorted(seen)]

    if patches:
        fig.legend(handles=patches, loc="lower center",
                   bbox_to_anchor=(0.5, -0.02), ncol=min(len(patches), 6),
                   fontsize=7, framealpha=0.9, title="Legend", title_fontsize=8)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# ============================================================
# GROUPING
# ============================================================
def suggest_group(name: str) -> str:
    """Strip trailing numbers/spaces to get a run-group key."""
    key = re.sub(r'[\s_-]*\d+\s*$', '', name.strip()).strip()
    return key or name.strip()


# ============================================================
# CONTROL COUNTING (pre-pass)
# ============================================================
def count_controls(groups_seen, group_map, use_k562, use_supermix):
    """Pre-count controls needed to compute row-band offsets."""
    k562_used     = len(groups_seen) if use_k562     else 0
    supermix_used = len(groups_seen) if use_supermix else 0
    blank_used    = sum(1 + len(split_groups(group_map[g])) for g in groups_seen)

    k562_spares     = max(3, math.ceil(k562_used     * 0.10)) if use_k562     else 0
    supermix_spares = max(3, math.ceil(supermix_used * 0.10)) if use_supermix else 0
    blank_spares    = max(3, math.ceil(blank_used    * 0.10))

    return {
        "k562_total":     k562_used + k562_spares,
        "supermix_total": supermix_used + supermix_spares,
        "blank_total":    blank_used + blank_spares,
        "k562_used":      k562_used,
        "supermix_used":  supermix_used,
        "blank_used":     blank_used,
        "k562_spares":    k562_spares,
        "supermix_spares":supermix_spares,
        "blank_spares":   blank_spares,
    }


# ============================================================
# CONTROL ALLOCATOR — row-banded, dynamic multi-slot
# ============================================================
class ControlAllocator:
    """
    Assigns K562 / Supermix / Blank to separate row bands so each
    type occupies contiguous rows (no empty-row gaps between types).
    Overflows to Slot3, Slot4, ... as needed.
    """
    def __init__(self, ctrl_counts: dict, use_k562: bool, use_supermix: bool,
                 start_slot: int = 2):
        k562_rows     = math.ceil(ctrl_counts["k562_total"]     / 12) if use_k562     else 0
        supermix_rows = math.ceil(ctrl_counts["supermix_total"] / 12) if use_supermix else 0

        # 1-indexed absolute position offsets per type
        self._offsets = {
            "K562":     0,
            "Supermix": k562_rows * 12,
            "Blank":    (k562_rows + supermix_rows) * 12,
        }
        self._counts     = {"K562": 0, "Supermix": 0, "Blank": 0}
        self._start_slot = start_slot
        self.entries     = []   # (slot, pos, ctype, sid, in_queue)

    def add(self, ctype: str, sid: str, in_queue: bool = True) -> str:
        self._counts[ctype] += 1
        abs_pos  = self._offsets[ctype] + self._counts[ctype]
        slot     = self._start_slot + (abs_pos - 1) // 96
        slot_pos = ((abs_pos - 1) % 96) + 1
        self.entries.append((slot, slot_pos, ctype, sid, in_queue))
        return f"Slot{slot}:{slot_pos}"

    def plates(self) -> dict:
        """Return {slot_num: {row: {col: sid}}} for all used ctrl slots."""
        grids = {}
        for slot, pos, ctype, sid, _ in self.entries:
            if slot not in grids:
                grids[slot] = {r: {c: "" for c in COLS} for r in ROWS}
            r, c = index_to_well(pos)
            grids[slot][r][c] = sid
        return grids


# ============================================================
# QUEUE BUILDER
# ============================================================
def build_queue_core(csv_bytes: bytes, group_assignments: dict, p: dict) -> dict:
    """
    group_assignments: {roi_name: group_label}
    Controls overflow dynamically into Slot3, Slot4 as needed.
    """
    date          = p["date"]
    initials      = p["initials"]
    lc_short      = p["lc_short"]
    ms_short      = p["ms_short"]
    sample_load   = p["sample_load"]
    k562_load     = p.get("k562_load", "")
    supermix_load = p.get("supermix_load", "")
    use_k562      = p["use_k562"]
    use_supermix  = p["use_supermix"]
    sep_method    = p["sep_method"]
    inj_method    = p["inj_method"]
    ms_method     = p["ms_method"]
    proc_method   = p["proc_method"]
    sample_path   = p["sample_path"]
    blank_path    = p["blank_path"]
    stem          = p["stem"]

    text   = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    reader.fieldnames = [h.strip() for h in reader.fieldnames]
    all_rows = list(reader)

    dropout_wells = {r["Well_ID"].strip() for r in all_rows
                     if r.get("Dropout {Y/N}", "").strip().upper() == "Y"}
    samples = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    # Build groups preserving CSV order
    groups_seen, group_map = [], {}
    for row in samples:
        roi = row["ROI"].strip()
        grp = group_assignments.get(roi, suggest_group(roi))
        if grp not in group_map:
            group_map[grp] = []
            groups_seen.append(grp)
        group_map[grp].append(row)

    # Pre-count controls to set row-band offsets
    ctrl_counts = count_controls(groups_seen, group_map, use_k562, use_supermix)
    alloc       = ControlAllocator(ctrl_counts, use_k562, use_supermix, start_slot=2)
    counts      = {"K562": 0, "Supermix": 0, "Blank": 0}
    queue       = []

    def add_k562():
        counts["K562"] += 1
        sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{counts['K562']}"
        vial = alloc.add("K562", sid)
        queue.append(make_row(vial, sid, sample_path, sep_method, inj_method, ms_method, proc_method))

    def add_supermix():
        counts["Supermix"] += 1
        sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{counts['Supermix']}"
        vial = alloc.add("Supermix", sid)
        queue.append(make_row(vial, sid, sample_path, sep_method, inj_method, ms_method, proc_method))

    def add_blank():
        counts["Blank"] += 1
        sid  = f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{counts['Blank']}"
        vial = alloc.add("Blank", sid)
        queue.append(make_row(vial, sid, blank_path, sep_method, inj_method, ms_method, proc_method))

    for grp in groups_seen:
        if use_k562:     add_k562()
        if use_supermix: add_supermix()
        add_blank()
        for batch in split_groups(group_map[grp]):
            for row in batch:
                roi  = row["ROI"].strip()
                vial = well_to_slot1(row["Well_ID"].strip())
                sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{sample_load}_{roi}"
                queue.append(make_row(vial, sid, sample_path,
                                      sep_method, inj_method, ms_method, proc_method))
            add_blank()

    # Spares (pre-counted — added to ctrl slots, not in queue)
    k562_spares     = ctrl_counts["k562_spares"]
    supermix_spares = ctrl_counts["supermix_spares"]
    blank_spares    = ctrl_counts["blank_spares"]

    for i in range(1, k562_spares + 1):
        n   = counts["K562"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{n}_spare"
        alloc.add("K562", sid, in_queue=False)
    for i in range(1, supermix_spares + 1):
        n   = counts["Supermix"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{n}_spare"
        alloc.add("Supermix", sid, in_queue=False)
    for i in range(1, blank_spares + 1):
        n   = counts["Blank"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{n}_spare"
        alloc.add("Blank", sid, in_queue=False)

    # Queue XLSX
    try:
        import openpyxl
        wb = openpyxl.Workbook()
        ws = wb.active
        ws.append(QUEUE_COLS)
        for row in queue:
            ws.append([row.get(col, "") for col in QUEUE_COLS])
        xlsx_buf = io.BytesIO()
        wb.save(xlsx_buf)
        queue_xlsx = xlsx_buf.getvalue()
    except ImportError:
        queue_xlsx = None

    # Slot1 - samples
    slot1_grid    = {r: {c: "" for c in COLS} for r in ROWS}
    well_to_group = {}
    for row in all_rows:
        w_id = row["Well_ID"].strip()
        if w_id:
            roi = row["ROI"].strip()
            slot1_grid[w_id[0]][int(w_id[1:])] = roi
            well_to_group[roi] = group_assignments.get(roi, suggest_group(roi))

    n_g = max(len(groups_seen), 1)
    group_colors = {}
    for i, g in enumerate(groups_seen):
        idx = round(i / n_g * 20)
        group_colors[g] = GREY_REPLACEMENTS.get(idx, cm.tab20(i / n_g))

    slot1_color_map   = {}
    slot1_label_color = {}
    for r in ROWS:
        for c in COLS:
            roi = slot1_grid[r][c]
            if roi:
                grp    = well_to_group.get(roi, "")
                wid    = f"{r}{c}"
                if wid in dropout_wells:
                    slot1_color_map[roi]   = "#dddddd"
                    slot1_label_color[roi] = "red"
                else:
                    slot1_color_map[roi]   = group_colors.get(grp, "white")
                    slot1_label_color[roi] = "black"

    roi_to_group = {roi: well_to_group.get(roi, roi) for roi in slot1_color_map}
    slot1_png    = plot_plate(slot1_grid, slot1_color_map, f"Slot1 - Samples ({stem})",
                              legend_group_map=roi_to_group,
                              label_color_map=slot1_label_color)

    slot1_buf = io.StringIO()
    wtr = csv.writer(slot1_buf)
    wtr.writerow([""] + COLS)
    for r in ROWS:
        wtr.writerow([r] + [slot1_grid[r][c] for c in COLS])

    # Control slots
    ctrl_color_map = {}
    ctrl_label_map = {}
    for _, _, ctype, sid, in_queue in alloc.entries:
        ctrl_color_map[sid] = (CTRL_COLORS if in_queue else CTRL_COLORS_SPARE).get(ctype, "white")
        num = sid.split("_")[-1]
        ctrl_label_map[sid] = f"{ctype}\n{num}"
    sid_to_type = {sid: ctype for _, _, ctype, sid, _ in alloc.entries}

    ctrl_slot_outputs = {}
    for slot_num, grid in alloc.plates().items():
        png = plot_plate(grid, ctrl_color_map, f"Slot{slot_num} - Controls ({stem})",
                         label_map=ctrl_label_map, legend_group_map=sid_to_type)
        buf = io.StringIO()
        wtr = csv.writer(buf)
        wtr.writerow([""] + COLS)
        for r in ROWS:
            wtr.writerow([r] + [grid[r][c] for c in COLS])
        ctrl_slot_outputs[slot_num] = {"png": png, "csv": buf.getvalue().encode("utf-8")}

    return {
        "queue_xlsx":      queue_xlsx,
        "slot1_csv":       slot1_buf.getvalue().encode("utf-8"),
        "slot1_png":       slot1_png,
        "ctrl_slots":      ctrl_slot_outputs,
        "n_queue":         len(queue),
        "counts":          counts,
        "k562_spares":     k562_spares,
        "supermix_spares": supermix_spares,
        "blank_spares":    blank_spares,
        "stem":            stem,
        "groups":          groups_seen,
    }


def build_zip(res: dict) -> bytes:
    buf  = io.BytesIO()
    stem = res["stem"]
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        if res["queue_xlsx"]:
            z.writestr(f"{stem}_queue.xlsx", res["queue_xlsx"])
        z.writestr(f"{stem}_slot1.csv", res["slot1_csv"])
        z.writestr(f"{stem}_slot1.png", res["slot1_png"])
        for slot_num, data in res["ctrl_slots"].items():
            z.writestr(f"{stem}_slot{slot_num}.csv", data["csv"])
            z.writestr(f"{stem}_slot{slot_num}.png", data["png"])
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_ms_queue_tab():
    st.header("MS Sample Queue")
    st.caption("Generates Bruker timsTOF queue (XLSX + plate maps) from sample list CSV.")

    for key in ("msq_results", "msq_zip", "msq_last", "msq_group_assignments", "msq_csv_hash"):
        if key not in st.session_state:
            st.session_state[key] = None

    date = datetime.date.today().strftime("%Y%m%d")

    # Parameters
    c1, c2, c3 = st.columns([1, 2, 2])
    initials = c1.text_input("Initials", value="FK", key="msq_initials")

    lc_sel   = c2.selectbox("LC method", LC_OPTIONS, key="msq_lc_sel")
    lc_short = c2.text_input("Custom LC short name", value="", key="msq_lc_custom") \
               if lc_sel == "Custom" else lc_sel

    ms_sel   = c3.selectbox("MS method", MS_OPTIONS, key="msq_ms_sel")
    ms_short = c3.text_input("Custom MS short name", value="", key="msq_ms_custom") \
               if ms_sel == "Custom" else ms_sel

    c4, c5, c6, c7, c8 = st.columns(5)
    sample_load   = c4.text_input("Sample load",   value="20ul", key="msq_sload")
    use_k562      = c5.checkbox("K562",             value=True,   key="msq_k562")
    k562_load     = c6.text_input("K562 load",      value="1ng",  key="msq_k562l") if use_k562     else ""
    use_supermix  = c7.checkbox("Supermix",         value=True,   key="msq_smix")
    supermix_load = c8.text_input("Supermix load",  value="20ng", key="msq_smixl") if use_supermix else ""

    with st.expander("Instrument method paths", expanded=False):
        default_sep              = LC_METHODS.get(lc_short, "")
        default_ms, default_proc = MS_METHODS.get(ms_short, ("", ""))
        sep_method  = st.text_area("Separation Method", value=default_sep,   height=60, key="msq_sep")
        inj_method  = st.text_input("Injection Method",  value="Standard",               key="msq_inj")
        ms_method   = st.text_area("MS Method",          value=default_ms,   height=60, key="msq_ms")
        proc_method = st.text_area("Processing Method",  value=default_proc, height=60, key="msq_proc")

    _year      = date[:4]
    _month_str = datetime.datetime.strptime(date, "%Y%m%d").strftime("%m %B")
    with st.expander("Data paths", expanded=True):
        cp1, cp2    = st.columns(2)
        sample_path = cp1.text_input("Sample / K562 / Supermix path",
                                     value=rf"D:\Data\{_year}\{_month_str}\Sample\{initials}",
                                     key="msq_spath")
        blank_path  = cp2.text_input("Blank path",
                                     value=rf"D:\Data\{_year}\{_month_str}\Blank",
                                     key="msq_bpath")

    st.divider()

    # Source
    pipe_csv  = st.session_state.get("t3_sample_list")
    pipe_stem = st.session_state.get("t3_stem")
    uploaded  = st.file_uploader("Upload sample list CSV", type=["csv"], key="msq_upload")

    if uploaded:
        csv_bytes = uploaded.getvalue()
        stem      = re.sub(r'^\d{8}', date, uploaded.name.replace(".csv", ""))
        if st.session_state.msq_last != uploaded.name:
            st.session_state.msq_results          = None
            st.session_state.msq_last             = uploaded.name
            st.session_state.msq_group_assignments = None
    elif pipe_csv is not None:
        csv_bytes = pipe_csv
        stem      = pipe_stem or "sample_list"
        st.info("Using sample list piped from Tab 3.")
    else:
        st.stop()

    # Reset groups if CSV content changed
    csv_hash = hash(csv_bytes)
    if st.session_state.msq_csv_hash != csv_hash:
        st.session_state.msq_group_assignments = None
        st.session_state.msq_results           = None
        st.session_state.msq_csv_hash          = csv_hash

    st.caption(f"Output stem: `{stem}`")

    # Parse samples for grouping editor
    text   = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    reader.fieldnames = [h.strip() for h in reader.fieldnames]
    all_rows = list(reader)
    samples  = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    if not samples:
        st.error("No active samples in CSV (all marked as dropouts).")
        return

    st.divider()
    st.subheader("Run Grouping")
    st.caption(
        "Groups are auto-suggested by stripping trailing numbers from sample names. "
        "Edit **Group** to reassign samples. Each group gets K562/Supermix/Blank at the start, "
        "then a Blank after every 6 samples."
    )

    confirmed_assignments = st.session_state.msq_group_assignments or {}

    init_data = [
        {
            "ROI":   r["ROI"].strip(),
            "Well":  r.get("Well_ID", "").strip(),
            "Group": confirmed_assignments.get(r["ROI"].strip(), suggest_group(r["ROI"].strip())),
        }
        for r in samples
    ]
    group_df = pd.DataFrame(init_data)

    edited_df = st.data_editor(
        group_df,
        column_config={
            "ROI":   st.column_config.TextColumn("ROI",   disabled=True),
            "Well":  st.column_config.TextColumn("Well",  disabled=True),
            "Group": st.column_config.TextColumn("Group", help="Edit to change group assignment"),
        },
        hide_index=True,
        use_container_width=True,
        key="msq_group_editor",
    )

    new_assignments = dict(zip(edited_df["ROI"], edited_df["Group"]))

    if st.button("Confirm grouping", type="primary", key="msq_confirm"):
        changed = new_assignments != st.session_state.msq_group_assignments
        st.session_state.msq_group_assignments = new_assignments
        if changed:
            st.session_state.msq_results = None
        st.rerun()

    if not st.session_state.msq_group_assignments:
        st.info("Confirm grouping above to continue.")
        st.stop()

    group_assignments = st.session_state.msq_group_assignments

    # Group summary
    summary = {}
    for roi, grp in group_assignments.items():
        summary.setdefault(grp, 0)
        summary[grp] += 1
    st.success("Groups: " + "  |  ".join(f"**{g}** ({n})" for g, n in summary.items()))

    st.divider()

    if st.button("Generate queue", type="primary", key="msq_gen"):
        params = dict(
            date=date, initials=initials, lc_short=lc_short, ms_short=ms_short,
            sample_load=sample_load, k562_load=k562_load, supermix_load=supermix_load,
            use_k562=use_k562, use_supermix=use_supermix,
            sep_method=sep_method, inj_method=inj_method,
            ms_method=ms_method, proc_method=proc_method,
            sample_path=sample_path, blank_path=blank_path,
            stem=stem,
        )
        with st.spinner("Generating..."):
            res = build_queue_core(csv_bytes, group_assignments, params)
            st.session_state.msq_results = res
            st.session_state.msq_zip     = build_zip(res)

    res = st.session_state.msq_results
    if not res:
        return

    c        = res["counts"]
    n_cslots = len(res["ctrl_slots"])
    st.success(
        f"{res['n_queue']} queue rows | "
        f"K562: {c['K562']} (+{res['k562_spares']} spare) | "
        f"Supermix: {c['Supermix']} (+{res['supermix_spares']} spare) | "
        f"Blank: {c['Blank']} (+{res['blank_spares']} spare) | "
        f"Control slots used: {n_cslots}"
    )

    # Plate maps — Slot1 + ctrl slots side by side (max 3 per row)
    all_plate_items = [("Slot1 - Samples", res["slot1_png"])] + [
        (f"Slot{sn} - Controls", data["png"])
        for sn, data in res["ctrl_slots"].items()
    ]
    for i in range(0, len(all_plate_items), 2):
        chunk = all_plate_items[i:i + 2]
        cols  = st.columns(len(chunk))
        for col, (title, png) in zip(cols, chunk):
            with col:
                st.subheader(title)
                st.image(png, use_container_width=True)

    st.subheader("Downloads")
    stem_out = res["stem"]

    def build_all_zip(res, msq_zip_bytes):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            if st.session_state.get("t1_geojson") and st.session_state.get("t1_stem"):
                z.writestr(f"{st.session_state.t1_stem}_reclassified.geojson",
                           st.session_state.t1_geojson)
            if st.session_state.get("t2_xml") and st.session_state.get("t2_stem"):
                z.writestr(f"{st.session_state.t2_stem}.xml", st.session_state.t2_xml)
            if st.session_state.get("t2_saw"):
                z.writestr("samples_and_wells.json",
                           json.dumps(st.session_state.t2_saw, indent=2).encode("utf-8"))
            if st.session_state.get("proc_result"):
                r3 = st.session_state.proc_result
                s3 = r3["stem"]
                z.writestr(f"{s3}_sorted.xml",      r3["sorted_xml"])
                z.writestr(f"{s3}_96wellplate.csv", r3["wellplate_csv"])
                z.writestr(f"{s3}_platemap.png",    st.session_state.proc_png)
            if st.session_state.get("t3_sample_list"):
                z.writestr(f"{stem_out}_sample_list.csv", st.session_state.t3_sample_list)
            with zipfile.ZipFile(io.BytesIO(msq_zip_bytes)) as msq_z:
                for name in msq_z.namelist():
                    z.writestr(name, msq_z.read(name))
        return buf.getvalue()

    dl_all, dl_ms = st.columns(2)
    dl_all.download_button(
        "Download all steps (zip)", build_all_zip(res, st.session_state.msq_zip),
        file_name=f"{stem_out}_all_steps.zip", mime="application/zip", type="primary"
    )
    dl_ms.download_button(
        "Download MS queue only (zip)", st.session_state.msq_zip,
        file_name=f"{stem_out}_ms_queue.zip", mime="application/zip"
    )

    st.markdown("---")

    # Individual file downloads
    dl_items = []
    if res["queue_xlsx"]:
        dl_items.append(("Queue XLSX", res["queue_xlsx"],
                         f"{stem_out}_queue.xlsx",
                         "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"))
    dl_items.append(("Slot1 CSV", res["slot1_csv"], f"{stem_out}_slot1.csv", "text/csv"))
    dl_items.append(("Slot1 PNG", res["slot1_png"], f"{stem_out}_slot1.png", "image/png"))
    for sn, data in res["ctrl_slots"].items():
        dl_items.append((f"Slot{sn} CSV", data["csv"], f"{stem_out}_slot{sn}.csv", "text/csv"))
        dl_items.append((f"Slot{sn} PNG", data["png"], f"{stem_out}_slot{sn}.png", "image/png"))

    cols = st.columns(min(len(dl_items), 5))
    for i, (label, data, fname, mime) in enumerate(dl_items):
        cols[i % 5].download_button(label, data, file_name=fname, mime=mime)
