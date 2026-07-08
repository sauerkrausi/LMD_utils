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
import streamlit as st

# ============================================================
# METHOD LOOKUP — edit these to add your own methods
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

ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))
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


def core_from_roi(name):
    parts = name.split("_")
    try:
        int(parts[-1]); int(parts[-2])
        return "_".join(parts[:-2]) or name
    except (ValueError, IndexError):
        return name


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
# FEASIBILITY CHECK
# ============================================================
def check_one_slot_feasible(csv_bytes: bytes, use_k562: bool, use_supermix: bool):
    text   = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    reader.fieldnames = [h.strip() for h in reader.fieldnames]
    all_rows = list(reader)
    samples  = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    cores_seen, core_map = [], {}
    for row in samples:
        core = core_from_roi(row["ROI"].strip())
        if core not in core_map:
            core_map[core] = []
            cores_seen.append(core)
        core_map[core].append(row)

    n_cores         = len(cores_seen)
    k562_needed     = n_cores if use_k562 else 0
    supermix_needed = n_cores if use_supermix else 0
    blank_needed    = n_cores + sum(len(split_groups(core_map[c])) for c in cores_seen)

    if k562_needed > 12:
        return False, f"K562 needs {k562_needed} vials (max 12 per row)"
    if supermix_needed > 12:
        return False, f"Supermix needs {supermix_needed} vials (max 12 per row)"
    if blank_needed > 12:
        return False, f"Blank needs {blank_needed} vials (max 12 per row)"

    n_ctrl_rows = 1 + (1 if use_supermix else 0) + (1 if use_k562 else 0)
    reserved    = set(ROWS[8 - n_ctrl_rows:])

    for row in samples:
        w = row.get("Well_ID", "").strip()
        if w and w[0].upper() in reserved:
            return False, f"Sample in well {w} conflicts with control row {w[0].upper()}"

    return True, (f"K562={k562_needed}, Supermix={supermix_needed}, Blank={blank_needed} "
                  f"| control rows: {', '.join(sorted(reserved))}")


# ============================================================
# QUEUE BUILDER
# ============================================================
def build_queue_core(csv_bytes: bytes, p: dict) -> dict:
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
    one_slot      = p.get("one_slot", False)

    text   = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    reader.fieldnames = [h.strip() for h in reader.fieldnames]
    all_rows = list(reader)

    dropout_wells = {r["Well_ID"].strip() for r in all_rows
                     if r.get("Dropout {Y/N}", "").strip().upper() == "Y"}
    samples = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    cores_seen, core_map = [], {}
    for row in samples:
        core = core_from_roi(row["ROI"].strip())
        if core not in core_map:
            core_map[core] = []
            cores_seen.append(core)
        core_map[core].append(row)

    n_cores         = len(cores_seen)
    k562_needed     = n_cores if use_k562 else 0
    supermix_needed = n_cores if use_supermix else 0
    blank_needed    = n_cores + sum(len(split_groups(core_map[c])) for c in cores_seen)

    k562_spares     = max(3, math.ceil(k562_needed     * 0.10)) if use_k562     else 0
    supermix_spares = max(3, math.ceil(supermix_needed * 0.10)) if use_supermix else 0
    blank_spares    = max(3, math.ceil(blank_needed    * 0.10))

    if one_slot:
        ctrl_types = []
        if use_k562:     ctrl_types.append("K562")
        if use_supermix: ctrl_types.append("Supermix")
        ctrl_types.append("Blank")
        ctrl_offsets = {ct: (7 - i) * 12 for i, ct in enumerate(reversed(ctrl_types))}
        k562_offset     = ctrl_offsets.get("K562", 0)
        supermix_offset = ctrl_offsets.get("Supermix", 0)
        blank_offset    = ctrl_offsets.get("Blank", 84)
        ctrl_slot       = "Slot1"
    else:
        k562_offset     = 0
        supermix_offset = 24
        blank_offset    = 48
        ctrl_slot       = "Slot2"

    counts        = {"K562": 0, "Supermix": 0, "Blank": 0}
    queue         = []
    slot2_entries = []

    def add_k562():
        counts["K562"] += 1
        pos = k562_offset + counts["K562"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{counts['K562']}"
        slot2_entries.append((pos, "K562", sid, True))
        queue.append(make_row(f"{ctrl_slot}:{pos}", sid, sample_path,
                              sep_method, inj_method, ms_method, proc_method))

    def add_supermix():
        counts["Supermix"] += 1
        pos = supermix_offset + counts["Supermix"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{counts['Supermix']}"
        slot2_entries.append((pos, "Supermix", sid, True))
        queue.append(make_row(f"{ctrl_slot}:{pos}", sid, sample_path,
                              sep_method, inj_method, ms_method, proc_method))

    def add_blank():
        counts["Blank"] += 1
        pos = blank_offset + counts["Blank"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{counts['Blank']}"
        slot2_entries.append((pos, "Blank", sid, True))
        queue.append(make_row(f"{ctrl_slot}:{pos}", sid, blank_path,
                              sep_method, inj_method, ms_method, proc_method))

    for core in cores_seen:
        if use_k562:     add_k562()
        if use_supermix: add_supermix()
        add_blank()
        for group in split_groups(core_map[core]):
            for row in group:
                roi  = row["ROI"].strip()
                vial = well_to_slot1(row["Well_ID"].strip())
                sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{sample_load}_{roi}"
                queue.append(make_row(vial, sid, sample_path,
                                      sep_method, inj_method, ms_method, proc_method))
            add_blank()

    for i in range(1, k562_spares + 1):
        n = counts["K562"] + i
        slot2_entries.append((k562_offset + n, "K562",
                              f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{n}_spare", False))
    for i in range(1, supermix_spares + 1):
        n = counts["Supermix"] + i
        slot2_entries.append((supermix_offset + n, "Supermix",
                              f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{n}_spare", False))
    for i in range(1, blank_spares + 1):
        n = counts["Blank"] + i
        slot2_entries.append((blank_offset + n, "Blank",
                              f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{n}_spare", False))

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

    # Slot1 grid + PNG
    slot1_grid   = {r: {c: "" for c in COLS} for r in ROWS}
    well_to_core = {}
    for row in all_rows:
        w_id = row["Well_ID"].strip()
        if w_id:
            roi = row["ROI"].strip()
            slot1_grid[w_id[0]][int(w_id[1:])] = roi
            well_to_core[roi] = core_from_roi(roi)

    slot1_buf = io.StringIO()
    w = csv.writer(slot1_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [slot1_grid[r][c] for c in COLS])

    n_c = max(len(cores_seen), 1)
    core_colors = {}
    for i, c in enumerate(cores_seen):
        idx = round(i / n_c * 20)
        core_colors[c] = GREY_REPLACEMENTS.get(idx, cm.tab20(i / n_c))

    slot1_color_map   = {}
    slot1_label_color = {}
    for r in ROWS:
        for c in COLS:
            roi = slot1_grid[r][c]
            if roi:
                core    = well_to_core.get(roi, "")
                well_id = f"{r}{c}"
                if well_id in dropout_wells:
                    slot1_color_map[roi]   = "#dddddd"
                    slot1_label_color[roi] = "red"
                else:
                    slot1_color_map[roi]   = core_colors.get(core, "white")
                    slot1_label_color[roi] = "black"

    roi_to_core = {roi: well_to_core.get(roi, roi) for roi in slot1_color_map}
    slot1_png = plot_plate(slot1_grid, slot1_color_map, f"Slot1 - Samples ({stem})",
                           legend_group_map=roi_to_core,
                           label_color_map=slot1_label_color)

    ctrl_color_map = {}
    ctrl_label_map = {}
    for pos, ctype, sid, in_queue in slot2_entries:
        ctrl_color_map[sid] = (CTRL_COLORS if in_queue else CTRL_COLORS_SPARE).get(ctype, "white")
        num = sid.split("_")[-1]
        ctrl_label_map[sid] = f"{ctype}\n{num}"
    sid_to_type = {sid: ctype for _, ctype, sid, _ in slot2_entries}

    if one_slot:
        combined_grid       = {r: {c: slot1_grid[r][c] for c in COLS} for r in ROWS}
        combined_color_map  = dict(slot1_color_map)
        combined_label_map  = {}
        combined_lcolor_map = dict(slot1_label_color)
        for pos, ctype, sid, _ in slot2_entries:
            r, c = index_to_well(pos)
            combined_grid[r][c]      = sid
            combined_color_map[sid]  = ctrl_color_map[sid]
            combined_label_map[sid]  = ctrl_label_map[sid]
            combined_lcolor_map[sid] = "black"
        combined_legend = {**roi_to_core, **sid_to_type}
        combined_png = plot_plate(combined_grid, combined_color_map,
                                  f"Combined plate ({stem})",
                                  label_map=combined_label_map,
                                  legend_group_map=combined_legend,
                                  label_color_map=combined_lcolor_map)
        combined_buf = io.StringIO()
        w = csv.writer(combined_buf)
        w.writerow([""] + COLS)
        for r in ROWS:
            w.writerow([r] + [combined_grid[r][c] for c in COLS])

        return {
            "queue_xlsx":      queue_xlsx,
            "slot1_csv":       combined_buf.getvalue().encode("utf-8"),
            "slot1_png":       combined_png,
            "slot2_csv":       None,
            "slot2_png":       None,
            "n_queue":         len(queue),
            "counts":          counts,
            "k562_spares":     k562_spares,
            "supermix_spares": supermix_spares,
            "blank_spares":    blank_spares,
            "stem":            stem,
            "one_slot":        True,
        }

    slot2_grid = {r: {c: "" for c in COLS} for r in ROWS}
    for pos, ctype, sid, _ in slot2_entries:
        r, c = index_to_well(pos)
        slot2_grid[r][c] = sid

    slot2_buf = io.StringIO()
    w = csv.writer(slot2_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [slot2_grid[r][c] for c in COLS])

    slot2_png = plot_plate(slot2_grid, ctrl_color_map, f"Slot2 - Controls ({stem})",
                           label_map=ctrl_label_map, legend_group_map=sid_to_type)

    return {
        "queue_xlsx":      queue_xlsx,
        "slot1_csv":       slot1_buf.getvalue().encode("utf-8"),
        "slot2_csv":       slot2_buf.getvalue().encode("utf-8"),
        "slot1_png":       slot1_png,
        "slot2_png":       slot2_png,
        "n_queue":         len(queue),
        "counts":          counts,
        "k562_spares":     k562_spares,
        "supermix_spares": supermix_spares,
        "blank_spares":    blank_spares,
        "stem":            stem,
        "one_slot":        False,
    }


def build_zip(res: dict) -> bytes:
    buf  = io.BytesIO()
    stem = res["stem"]
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        if res["queue_xlsx"]:
            z.writestr(f"{stem}_queue.xlsx", res["queue_xlsx"])
        if res.get("one_slot"):
            z.writestr(f"{stem}_combined.csv", res["slot1_csv"])
            z.writestr(f"{stem}_combined.png", res["slot1_png"])
        else:
            z.writestr(f"{stem}_slot1.csv", res["slot1_csv"])
            z.writestr(f"{stem}_slot2.csv", res["slot2_csv"])
            z.writestr(f"{stem}_slot1.png", res["slot1_png"])
            z.writestr(f"{stem}_slot2.png", res["slot2_png"])
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
def render_ms_queue_tab():
    st.header("MS Sample Queue")
    st.caption("Generates Bruker timsTOF queue (XLSX + plate maps) from sample list CSV.")

    # Session state
    for key in ("msq_results", "msq_zip", "msq_last"):
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
    sample_load   = c4.text_input("Sample load",   value="20ul",  key="msq_sload")
    use_k562      = c5.checkbox("K562",             value=True,    key="msq_k562")
    k562_load     = c6.text_input("K562 load",      value="1ng",   key="msq_k562l") if use_k562     else ""
    use_supermix  = c7.checkbox("Supermix",         value=True,    key="msq_smix")
    supermix_load = c8.text_input("Supermix load",  value="20ng",  key="msq_smixl") if use_supermix else ""

    with st.expander("Instrument method paths", expanded=False):
        default_sep            = LC_METHODS.get(lc_short, "")
        default_ms, default_proc = MS_METHODS.get(ms_short, ("", ""))
        sep_method  = st.text_area("Separation Method", value=default_sep,  height=60, key="msq_sep")
        inj_method  = st.text_input("Injection Method",  value="Standard",              key="msq_inj")
        ms_method   = st.text_area("MS Method",          value=default_ms,  height=60, key="msq_ms")
        proc_method = st.text_area("Processing Method",  value=default_proc,height=60, key="msq_proc")

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

    # Source: piped from Tab 3 or upload
    pipe_csv  = st.session_state.get("t3_sample_list")
    pipe_stem = st.session_state.get("t3_stem")
    uploaded  = st.file_uploader("Upload sample list CSV", type=["csv"], key="msq_upload")

    if uploaded:
        csv_bytes = uploaded.getvalue()
        stem      = re.sub(r'^\d{8}', date, uploaded.name.replace(".csv", ""))
        if st.session_state.msq_last != uploaded.name:
            st.session_state.msq_results = None
            st.session_state.msq_last    = uploaded.name
    elif pipe_csv is not None:
        csv_bytes = pipe_csv
        stem      = pipe_stem or "sample_list"
        st.info("Using sample list piped from Tab 3.")
    else:
        st.stop()

    st.caption(f"Output stem: `{stem}`")

    feasible, reason = check_one_slot_feasible(csv_bytes, use_k562, use_supermix)

    btn_cols = st.columns([2, 2, 3])
    gen_clicked      = btn_cols[0].button("Generate queue",    type="primary",    key="msq_gen")
    one_slot_clicked = btn_cols[1].button("Fit into one slot", type="secondary",  key="msq_1slot",
                                          disabled=not feasible)
    if feasible:
        btn_cols[2].success(f"One-slot possible: {reason}")
    else:
        btn_cols[2].caption(f"One-slot not possible: {reason}")

    def _run(one_slot: bool):
        params = dict(
            date=date, initials=initials, lc_short=lc_short, ms_short=ms_short,
            sample_load=sample_load, k562_load=k562_load, supermix_load=supermix_load,
            use_k562=use_k562, use_supermix=use_supermix,
            sep_method=sep_method, inj_method=inj_method,
            ms_method=ms_method, proc_method=proc_method,
            sample_path=sample_path, blank_path=blank_path,
            stem=stem, one_slot=one_slot,
        )
        with st.spinner("Generating..."):
            res = build_queue_core(csv_bytes, params)
            st.session_state.msq_results = res
            st.session_state.msq_zip     = build_zip(res)

    if gen_clicked:
        _run(False)
    if one_slot_clicked:
        _run(True)

    # Results
    res = st.session_state.msq_results
    if not res:
        return

    c = res["counts"]
    st.success(
        f"{res['n_queue']} queue rows | "
        f"K562: {c['K562']} (+{res['k562_spares']} spare) | "
        f"Supermix: {c['Supermix']} (+{res['supermix_spares']} spare) | "
        f"Blank: {c['Blank']} (+{res['blank_spares']} spare)"
    )

    if res.get("one_slot"):
        st.subheader("Combined Plate — Samples + Controls")
        st.image(res["slot1_png"], use_container_width=True)
    else:
        col1, col2 = st.columns(2)
        with col1:
            st.subheader("Slot 1 — Samples")
            st.image(res["slot1_png"], use_container_width=True)
        with col2:
            st.subheader("Slot 2 — Controls")
            st.image(res["slot2_png"], use_container_width=True)

    st.subheader("Downloads")
    stem_out = res["stem"]

    # All-steps zip: bundle everything available from tabs 1-4
    def build_all_zip(res: dict, msq_zip_bytes: bytes) -> bytes:
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            # Tab 1
            if st.session_state.get("t1_geojson") and st.session_state.get("t1_stem"):
                z.writestr(f"{st.session_state.t1_stem}_reclassified.geojson",
                           st.session_state.t1_geojson)
            # Tab 2
            if st.session_state.get("t2_xml") and st.session_state.get("t2_stem"):
                z.writestr(f"{st.session_state.t2_stem}.xml", st.session_state.t2_xml)
            if st.session_state.get("t2_saw"):
                z.writestr("samples_and_wells.json",
                           json.dumps(st.session_state.t2_saw, indent=2).encode("utf-8"))
            # Tab 3
            if st.session_state.get("proc_result"):
                r3   = st.session_state.proc_result
                s3   = r3["stem"]
                z.writestr(f"{s3}_sorted.xml",        r3["sorted_xml"])
                z.writestr(f"{s3}_96wellplate.csv",   r3["wellplate_csv"])
                z.writestr(f"{s3}_platemap.png",      st.session_state.proc_png)
            # Always include the current (possibly dropout-updated) sample list
            if st.session_state.get("t3_sample_list"):
                z.writestr(f"{stem_out}_sample_list.csv", st.session_state.t3_sample_list)
            # Tab 4
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

    if res.get("one_slot"):
        cols = st.columns(3)
        if res["queue_xlsx"]:
            cols[0].download_button("Queue XLSX", res["queue_xlsx"],
                                    file_name=f"{stem_out}_queue.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        cols[1].download_button("Plate CSV", res["slot1_csv"],
                                file_name=f"{stem_out}_combined.csv", mime="text/csv")
        cols[2].download_button("Plate PNG", res["slot1_png"],
                                file_name=f"{stem_out}_combined.png", mime="image/png")
    else:
        cols = st.columns(5)
        if res["queue_xlsx"]:
            cols[0].download_button("Queue XLSX", res["queue_xlsx"],
                                    file_name=f"{stem_out}_queue.xlsx",
                                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        cols[1].download_button("Slot1 CSV", res["slot1_csv"],
                                file_name=f"{stem_out}_slot1.csv", mime="text/csv")
        cols[2].download_button("Slot2 CSV", res["slot2_csv"],
                                file_name=f"{stem_out}_slot2.csv", mime="text/csv")
        cols[3].download_button("Slot1 PNG", res["slot1_png"],
                                file_name=f"{stem_out}_slot1.png", mime="image/png")
        cols[4].download_button("Slot2 PNG", res["slot2_png"],
                                file_name=f"{stem_out}_slot2.png", mime="image/png")
