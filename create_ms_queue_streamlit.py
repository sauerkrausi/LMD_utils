"""
create_ms_queue_streamlit.py
============================
Streamlit UI for the MS sample queue generator.

OVERVIEW
--------
1. User configures run parameters in the sidebar
2. Uploads sample list CSV (Core, ROI, Well_ID, Dropout columns)
3. Computed data paths shown for verification before processing
4. Generates queue CSV/XLSX + Slot1/Slot2 CSV + PNG plate maps
5. All outputs downloadable individually or as a zip

Deploy locally:  streamlit run create_ms_queue_streamlit.py
"""

import io
import csv
import math
import re
import zipfile
import datetime

import streamlit as st
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.cm as cm

# ============================================================
# PAGE CONFIG
# ============================================================
st.set_page_config(page_title="MS Queue Generator", layout="wide")

st.markdown(
    "<div style='position:fixed;bottom:10px;right:16px;font-size:12px;color:#888;'>"
    "Questions / issues: <a href='https://github.com/sauerkrausi/LMD_utils' target='_blank'>"
    "github.com/sauerkrausi/LMD_utils</a></div>",
    unsafe_allow_html=True,
)

# ============================================================
# CONSTANTS
# ============================================================
ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))
GROUP_SIZE = 6

QUEUE_COLS = [
    "Vial", "Sample ID", "Method Set", "Separation Method",
    "Injection Method", "MS Method", "Processing Method",
    "Sample Type", "Volume [µl]", "Data Path", "Run Automated Processing",
]

LC_OPTIONS    = ["WhisperZOOM40", "WhisperZOOM20", "WhisperZOOM100",
                 "30SPD", "60SPD", "100SPD", "Custom"]
MS_OPTIONS    = ["diaPASEF", "ddaPASEF", "Custom"]

# ============================================================
# HELPERS
# ============================================================
def well_to_slot1(well_id):
    row = ord(well_id[0].upper()) - ord('A')
    col = int(well_id[1:])
    return f"Slot1.{row * 12 + col}"

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
            display = (label_map or {}).get(label, label)
            if display:
                fsize      = 4.5 if len(display) > 10 else 5.5
                font_color = (label_color_map or {}).get(label, "black")
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
# CORE QUEUE BUILDER  (in-memory, no disk I/O)
# ============================================================
def build_queue_core(csv_bytes: bytes, p: dict) -> dict:
    """
    p: dict of all run parameters
    Returns dict of output bytes + metadata.
    """
    date         = p["date"]
    initials     = p["initials"]
    lc_short     = p["lc_short"]
    ms_short     = p["ms_short"]
    sample_load  = p["sample_load"]
    k562_load    = p.get("k562_load", "")
    supermix_load = p.get("supermix_load", "")
    use_k562     = p["use_k562"]
    use_supermix = p["use_supermix"]
    sep_method   = p["sep_method"]
    inj_method   = p["inj_method"]
    ms_method    = p["ms_method"]
    proc_method  = p["proc_method"]
    sample_path  = p["sample_path"]
    blank_path   = p["blank_path"]
    stem         = p["stem"]

    # Read CSV
    text = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    reader.fieldnames = [h.strip() for h in reader.fieldnames]
    all_rows = list(reader)

    dropout_wells = {r["Well_ID"].strip() for r in all_rows
                     if r.get("Dropout {Y/N}", "").strip().upper() == "Y"}
    samples = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    cores_seen, core_map = [], {}
    for row in samples:
        core = row["Core"].strip()
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

    k562_offset     = 0
    supermix_offset = 24
    blank_offset    = 48

    counts        = {"K562": 0, "Supermix": 0, "Blank": 0}
    queue         = []
    slot2_entries = []  # (pos, type, sid, in_queue)

    def add_k562():
        counts["K562"] += 1
        pos = k562_offset + counts["K562"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{counts['K562']}"
        slot2_entries.append((pos, "K562", sid, True))
        queue.append(make_row(f"Slot2.{pos}", sid, sample_path,
                              sep_method, inj_method, ms_method, proc_method))

    def add_supermix():
        counts["Supermix"] += 1
        pos = supermix_offset + counts["Supermix"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{counts['Supermix']}"
        slot2_entries.append((pos, "Supermix", sid, True))
        queue.append(make_row(f"Slot2.{pos}", sid, sample_path,
                              sep_method, inj_method, ms_method, proc_method))

    def add_blank():
        counts["Blank"] += 1
        pos = blank_offset + counts["Blank"]
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{counts['Blank']}"
        slot2_entries.append((pos, "Blank", sid, True))
        queue.append(make_row(f"Slot2.{pos}", sid, blank_path,
                              sep_method, inj_method, ms_method, proc_method))

    for core in cores_seen:
        if use_k562:
            add_k562()
        if use_supermix:
            add_supermix()
        add_blank()
        for group in split_groups(core_map[core]):
            for row in group:
                roi  = row["ROI"].strip()
                vial = well_to_slot1(row["Well_ID"].strip())
                sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{sample_load}_{roi}"
                queue.append(make_row(vial, sid, sample_path,
                                      sep_method, inj_method, ms_method, proc_method))
            add_blank()

    # Spare vials (plate only)
    for i in range(1, k562_spares + 1):
        n = counts["K562"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{n}_spare"
        slot2_entries.append((k562_offset + n, "K562", sid, False))
    for i in range(1, supermix_spares + 1):
        n = counts["Supermix"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{n}_spare"
        slot2_entries.append((supermix_offset + n, "Supermix", sid, False))
    for i in range(1, blank_spares + 1):
        n = counts["Blank"] + i
        sid = f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{n}_spare"
        slot2_entries.append((blank_offset + n, "Blank", sid, False))

    # Queue CSV
    queue_buf = io.StringIO()
    w = csv.DictWriter(queue_buf, fieldnames=QUEUE_COLS)
    w.writeheader()
    w.writerows(queue)

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

    # Slot1 CSV + PNG
    slot1_grid    = {r: {c: "" for c in COLS} for r in ROWS}
    well_to_core  = {}
    for row in all_rows:
        w_id = row["Well_ID"].strip()
        if w_id:
            slot1_grid[w_id[0]][int(w_id[1:])] = row["ROI"].strip()
            well_to_core[row["ROI"].strip()] = row["Core"].strip()

    slot1_buf = io.StringIO()
    w = csv.writer(slot1_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [slot1_grid[r][c] for c in COLS])

    GREY_REPLACEMENTS = {14: "#777777", 15: "#555555"}
    n_c = max(len(cores_seen), 1)
    core_colors = {}
    for i, c in enumerate(cores_seen):
        idx = round(i / n_c * 20)
        core_colors[c] = GREY_REPLACEMENTS.get(idx, cm.tab20(i / n_c))

    slot1_color_map  = {}
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

    # Slot2 CSV + PNG
    slot2_grid      = {r: {c: "" for c in COLS} for r in ROWS}
    slot2_type_grid = {r: {c: "" for c in COLS} for r in ROWS}
    for pos, ctype, sid, _ in slot2_entries:
        r, c = index_to_well(pos)
        slot2_grid[r][c]      = sid
        slot2_type_grid[r][c] = ctype

    slot2_buf = io.StringIO()
    w = csv.writer(slot2_buf)
    w.writerow([""] + COLS)
    for r in ROWS:
        w.writerow([r] + [slot2_grid[r][c] for c in COLS])

    CTRL_COLORS       = {"K562": "#4C9BE8", "Supermix": "#F4A261", "Blank": "#B7E4C7"}
    CTRL_COLORS_SPARE = {"K562": "#C5DDF7", "Supermix": "#FAE0C8", "Blank": "#E6F7EC"}
    slot2_color_map = {}
    slot2_label_map = {}
    for pos, ctype, sid, in_queue in slot2_entries:
        slot2_color_map[sid] = (CTRL_COLORS if in_queue else CTRL_COLORS_SPARE).get(ctype, "white")
        num = sid.split("_")[-1]
        slot2_label_map[sid] = f"{ctype}\n{num}"
    sid_to_type = {sid: ctype for _, ctype, sid, _ in slot2_entries}
    slot2_png = plot_plate(slot2_grid, slot2_color_map, f"Slot2 - Controls ({stem})",
                           label_map=slot2_label_map, legend_group_map=sid_to_type)

    return {
        "queue_csv":    queue_buf.getvalue().encode("utf-8"),
        "queue_xlsx":   queue_xlsx,
        "slot1_csv":    slot1_buf.getvalue().encode("utf-8"),
        "slot2_csv":    slot2_buf.getvalue().encode("utf-8"),
        "slot1_png":    slot1_png,
        "slot2_png":    slot2_png,
        "n_queue":      len(queue),
        "counts":       counts,
        "k562_spares":  k562_spares,
        "supermix_spares": supermix_spares,
        "blank_spares": blank_spares,
        "stem":         stem,
    }

def build_zip(res: dict) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{res['stem']}_queue.csv",   res["queue_csv"])
        if res["queue_xlsx"]:
            z.writestr(f"{res['stem']}_queue.xlsx", res["queue_xlsx"])
        z.writestr(f"{res['stem']}_slot1.csv",   res["slot1_csv"])
        z.writestr(f"{res['stem']}_slot2.csv",   res["slot2_csv"])
        z.writestr(f"{res['stem']}_slot1.png",   res["slot1_png"])
        z.writestr(f"{res['stem']}_slot2.png",   res["slot2_png"])
    return buf.getvalue()

# ============================================================
# SESSION STATE INIT
# ============================================================
for key in ("results", "zip_bytes", "last_file"):
    if key not in st.session_state:
        st.session_state[key] = None

# ============================================================
# MAIN UI — inputs at top
# ============================================================
st.title("MS Queue Generator")

date = datetime.date.today().strftime("%Y%m%d")

# Row 1: Initials | LC_SHORT | MS_SHORT
c1, c2, c3 = st.columns([1, 2, 2])
initials = c1.text_input("Initials", value="FK")

lc_sel   = c2.selectbox("LC method (LC_SHORT)", LC_OPTIONS)
lc_short = c2.text_input("Custom LC_SHORT", value="") if lc_sel == "Custom" else lc_sel

ms_sel   = c3.selectbox("MS method (MS_SHORT)", MS_OPTIONS)
ms_short = c3.text_input("Custom MS_SHORT", value="") if ms_sel == "Custom" else ms_sel

# Row 2: loads + controls
c4, c5, c6, c7, c8 = st.columns([1, 1, 1, 1, 1])
sample_load  = c4.text_input("Sample load", value="1ng")
use_k562     = c5.checkbox("Include K562", value=True)
k562_load    = c6.text_input("K562 load",    value="1ng")  if use_k562    else ""
use_supermix = c7.checkbox("Include Supermix", value=True)
supermix_load = c8.text_input("Supermix load", value="20ng") if use_supermix else ""

st.divider()

# Method paths — collapsible
with st.expander("Instrument method paths", expanded=False):
    sep_method  = st.text_area("Separation Method", height=60,
        value=r"D:\Methods\LC_Methods\Evosep\WhisperZOOM_40_SPD_32p5min.m?HyStar_LC")
    inj_method  = st.text_input("Injection Method", value="Standard")
    ms_method   = st.text_area("MS Method", height=60,
        value=(r"D:\Methods\MS_Methods\DIA\Farah\TimsControl methods"
               r"\DIA_PASEF_Var_windows_test4_pydiAID_300to1200_80PASEF_scans_-05shift"
               r".proteoscape.m?OtofImpacTEMControl"))
    proc_method = st.text_area("Processing Method", height=60,
        value=(r"D:\Methods\MS_Methods\DIA\Farah\TimsControl methods"
               r"\DIA_PASEF_Var_windows_test4_pydiAID_300to1200_80PASEF_scans_-05shift"
               r".proteoscape.m?DataAnalysis"))

# Data paths — editable, pre-filled from date + initials
_year      = date[:4]
_month_str = datetime.datetime.strptime(date, "%Y%m%d").strftime("%m %B")
with st.expander("Data paths (verify and edit if needed)", expanded=True):
    cp1, cp2 = st.columns(2)
    sample_path = cp1.text_input("Sample / K562 / Supermix path",
                                 value=rf"D:\Data\{_year}\{_month_str}\Sample\{initials}")
    blank_path  = cp2.text_input("Blank path",
                                 value=rf"D:\Data\{_year}\{_month_str}\Blank")

# File upload
uploaded = st.file_uploader("Upload sample list CSV", type=["csv"])

if uploaded and uploaded.name != st.session_state.last_file:
    st.session_state.results   = None
    st.session_state.zip_bytes = None
    st.session_state.last_file = uploaded.name

if uploaded:
    stem = re.sub(r'^\d{8}', date, uploaded.name.replace(".csv", ""))
    st.caption(f"Output stem: `{stem}`")

    if st.button("Generate queue", type="primary"):
        params = dict(
            date=date, initials=initials, lc_short=lc_short, ms_short=ms_short,
            sample_load=sample_load, k562_load=k562_load, supermix_load=supermix_load,
            use_k562=use_k562, use_supermix=use_supermix,
            sep_method=sep_method, inj_method=inj_method,
            ms_method=ms_method, proc_method=proc_method,
            sample_path=sample_path, blank_path=blank_path, stem=stem,
        )
        with st.spinner("Generating..."):
            res = build_queue_core(uploaded.getvalue(), params)
            st.session_state.results   = res
            st.session_state.zip_bytes = build_zip(res)

# ============================================================
# SIDEBAR — settings summary (read-only)
# ============================================================
with st.sidebar:
    st.header("Current settings")
    st.markdown(f"**Date:** {date}")
    st.markdown(f"**Initials:** {initials}")
    st.markdown(f"**LC:** {lc_short or '—'}")
    st.markdown(f"**MS:** {ms_short or '—'}")
    st.divider()
    st.markdown(f"**Sample load:** {sample_load}")
    if use_k562:
        st.markdown(f"**K562 load:** {k562_load}")
    else:
        st.markdown("**K562:** off")
    if use_supermix:
        st.markdown(f"**Supermix load:** {supermix_load}")
    else:
        st.markdown("**Supermix:** off")
    st.divider()
    st.markdown(f"**Sample path:**\n\n`{sample_path}`")
    st.markdown(f"**Blank path:**\n\n`{blank_path}`")

# ============================================================
# RESULTS
# ============================================================
if st.session_state.results:
    res  = st.session_state.results
    stem = res["stem"]
    c    = res["counts"]

    st.success(
        f"Done — {res['n_queue']} queue rows | "
        f"K562: {c['K562']} (+{res['k562_spares']} spare) | "
        f"Supermix: {c['Supermix']} (+{res['supermix_spares']} spare) | "
        f"Blank: {c['Blank']} (+{res['blank_spares']} spare)"
    )

    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Slot 1 — Samples")
        st.image(res["slot1_png"], use_container_width=True)
    with col2:
        st.subheader("Slot 2 — Controls")
        st.image(res["slot2_png"], use_container_width=True)

    st.subheader("Downloads")
    st.download_button("⬇ Download all (zip)", st.session_state.zip_bytes,
                       file_name=f"{stem}_ms_queue.zip", mime="application/zip",
                       type="primary")
    st.markdown("---")

    cols = st.columns(6)
    cols[0].download_button("Queue CSV",   res["queue_csv"],
                            file_name=f"{stem}_queue.csv",  mime="text/csv")
    if res["queue_xlsx"]:
        cols[1].download_button("Queue XLSX", res["queue_xlsx"],
                                file_name=f"{stem}_queue.xlsx",
                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    cols[2].download_button("Slot1 CSV",  res["slot1_csv"],
                            file_name=f"{stem}_slot1.csv",  mime="text/csv")
    cols[3].download_button("Slot2 CSV",  res["slot2_csv"],
                            file_name=f"{stem}_slot2.csv",  mime="text/csv")
    cols[4].download_button("Slot1 PNG",  res["slot1_png"],
                            file_name=f"{stem}_slot1.png",  mime="image/png")
    cols[5].download_button("Slot2 PNG",  res["slot2_png"],
                            file_name=f"{stem}_slot2.png",  mime="image/png")
