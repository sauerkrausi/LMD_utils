"""
L2 MS Queue
=============================
Input:  CSV with at minimum columns: sample_name, well
        (well in 96-well format, e.g. A1 .. H12)
Output: Bruker timsTOF queue XLSX + plate map PNGs (zip)

Run:
    cd L2_msQ
    pip install -r requirements.txt
    streamlit run app.py
"""

import csv
import datetime
import io
import math
import random
import re
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import openpyxl
import pandas as pd
import streamlit as st

st.set_page_config(page_title="L2 MS Queue", layout="wide")

# ============================================================
# CONSTANTS
# ============================================================
ROWS      = list("ABCDEFGH")
COLS      = list(range(1, 13))
GROUP_SIZE = 6

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

QUEUE_COLS = [
    "Vial", "Sample ID", "Method Set", "Separation Method",
    "Injection Method", "MS Method", "Processing Method",
    "Sample Type", "Volume [µl]", "Data Path", "Run Automated Processing",
]

CTRL_COLORS       = {"K562": "#4C9BE8", "Supermix": "#F4A261", "Blank": "#B7E4C7"}
CTRL_COLORS_SPARE = {"K562": "#C5DDF7", "Supermix": "#FAE0C8", "Blank": "#E6F7EC"}


# ============================================================
# HELPERS
# ============================================================
def suggest_group(name: str) -> str:
    key = re.sub(r'[\s_-]*\d+\s*$', '', name.strip()).strip()
    return key or name.strip()


def well_sort_key(w):
    try:
        return (ord(w[0].upper()) - ord('A'), int(w[1:]))
    except Exception:
        return (99, 99)


def index_to_well(index):
    i = index - 1
    return ROWS[i // 12], (i % 12) + 1


def well_to_vial(well_id: str, slot: str) -> str:
    row = ord(well_id[0].upper()) - ord('A')
    col = int(well_id[1:])
    return f"{slot}:{row * 12 + col}"


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


def make_row(vial, sid, data_path, sep, inj, ms, proc):
    return {
        "Vial": vial, "Sample ID": sid, "Method Set": "",
        "Separation Method": sep, "Injection Method": inj,
        "MS Method": ms, "Processing Method": proc,
        "Sample Type": "Sample", "Volume [µl]": 1,
        "Data Path": data_path, "Run Automated Processing": "False",
    }


# ============================================================
# PLATE VISUALIZATION
# ============================================================
def plot_plate_png(grid, color_map, title, label_map=None) -> bytes:
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
                fsize = 4.5 if len(display) > 10 else 5.5
                ax.text(x, y, display, ha="center", va="center",
                        fontsize=fsize, zorder=3, color="black", clip_on=True)
    for r_idx, r in enumerate(ROWS):
        ax.text(-0.55, 7 - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")
    seen = {grid[r][c]: color_map.get(grid[r][c], "white")
            for r in ROWS for c in COLS if grid[r][c]}
    if seen:
        patches = [mpatches.Patch(color=col, label=(label_map or {}).get(lbl, lbl))
                   for lbl, col in sorted(seen.items())]
        fig.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5, -0.02),
                   ncol=min(len(patches), 6), fontsize=7, framealpha=0.9)
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def plot_sample_plate(samples_ordered, title) -> bytes:
    """samples_ordered: list of {well, name, group}"""
    well_to = {s["well"]: s for s in samples_ordered}
    groups  = sorted({s["group"] for s in samples_ordered})
    palette = cm.tab20
    gcol    = {g: palette(i / max(len(groups), 1)) for i, g in enumerate(groups)}

    grid      = {r: {c: "" for c in COLS} for r in ROWS}
    color_map = {}
    for s in samples_ordered:
        w = s["well"]
        if len(w) >= 2 and w[0] in ROWS:
            grid[w[0]][int(w[1:])] = s["name"]
            color_map[s["name"]] = gcol.get(s["group"], "white")

    fig, ax = plt.subplots(figsize=(14, 8))
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(-0.5, 8.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    for r_idx, r in enumerate(ROWS):
        for c_idx, c in enumerate(COLS):
            x, y  = c_idx, 7 - r_idx
            name  = grid[r][c]
            color = color_map.get(name, "whitesmoke")
            edge  = "#444444" if name else "#aaaaaa"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.8, zorder=2))
            if name:
                ax.text(x, y, name, ha="center", va="center",
                        fontsize=4.5, zorder=3, color="black", clip_on=True)
    for r_idx, r in enumerate(ROWS):
        ax.text(-0.55, 7 - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(COLS):
        ax.text(c_idx, 8.0, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")
    patches = [mpatches.Patch(color=gcol[g], label=g) for g in groups]
    if patches:
        fig.legend(handles=patches, loc="lower center", bbox_to_anchor=(0.5, -0.02),
                   ncol=min(len(patches), 6), fontsize=7, framealpha=0.9)
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# ============================================================
# QUEUE BUILDER
# ============================================================
def build_queue(samples, group_assignments, p) -> dict:
    """
    samples: list of {name, well, group (original)}
    group_assignments: {name: group_label}
    p: params dict
    Returns: {queue_xlsx, sample_png, ctrl_png, n_queue, counts, spares}
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
    sample_slot   = p.get("sample_slot", "Slot1")
    ctrl_slot_start = int(p.get("ctrl_slot_start", 2))
    randomize_order = p.get("randomize_order", False)
    block_size    = int(p.get("block_size", GROUP_SIZE))
    run_seed      = int(p.get("run_seed", 42))

    # Build group map
    groups_seen, group_map = [], {}
    for s in samples:
        grp = group_assignments.get(s["name"], suggest_group(s["name"]))
        if grp not in group_map:
            group_map[grp] = []
            groups_seen.append(grp)
        group_map[grp].append(s)

    if randomize_order:
        rng = random.Random(run_seed)
        for grp in groups_seen:
            rng.shuffle(group_map[grp])

    # Count controls
    n_groups      = len(groups_seen)
    blank_used    = sum(1 + len(split_groups(group_map[g], block_size)) for g in groups_seen)
    k562_used     = n_groups if use_k562     else 0
    supermix_used = n_groups if use_supermix else 0
    k562_sp       = max(3, math.ceil(k562_used     * 0.10)) if use_k562     else 0
    supermix_sp   = max(3, math.ceil(supermix_used * 0.10)) if use_supermix else 0
    blank_sp      = max(3, math.ceil(blank_used    * 0.10))

    k562_rows     = math.ceil((k562_used + k562_sp)     / 12) if use_k562     else 0
    supermix_rows = math.ceil((supermix_used + supermix_sp) / 12) if use_supermix else 0
    blank_offset  = (k562_rows + supermix_rows) * 12

    counts      = {"K562": 0, "Supermix": 0, "Blank": 0}
    ctrl_entries = []

    def _ctrl_vial(ctype):
        offsets = {"K562": 0, "Supermix": k562_rows * 12, "Blank": blank_offset}
        counts[ctype] += 1
        abs_pos  = offsets[ctype] + counts[ctype]
        slot     = ctrl_slot_start + (abs_pos - 1) // 96
        slot_pos = ((abs_pos - 1) % 96) + 1
        return f"Slot{slot}:{slot_pos}", slot, slot_pos, ctype

    queue_rows = []

    def add_ctrl(ctype, sid, in_queue=True):
        vial, slot, pos, ct = _ctrl_vial(ctype)
        ctrl_entries.append((slot, pos, ct, sid, in_queue))
        if in_queue:
            dp = blank_path if ctype == "Blank" else sample_path
            queue_rows.append(make_row(vial, sid, dp, sep_method, inj_method, ms_method, proc_method))
        return vial

    for grp in groups_seen:
        n_k = counts["K562"]
        n_s = counts["Supermix"]
        n_b = counts["Blank"]
        if use_k562:
            add_ctrl("K562",
                f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{n_k+1}")
        if use_supermix:
            add_ctrl("Supermix",
                f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{n_s+1}")
        add_ctrl("Blank",
            f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{n_b+1}")
        for batch in split_groups(group_map[grp], block_size):
            for s in batch:
                sid  = f"{date}_{initials}_{lc_short}_{ms_short}_{sample_load}_{s['name']}"
                vial = well_to_vial(s["well"], sample_slot)
                queue_rows.append(make_row(vial, sid, sample_path,
                                           sep_method, inj_method, ms_method, proc_method))
            nb = counts["Blank"]
            add_ctrl("Blank",
                f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{nb+1}")

    # Spares
    for i in range(1, k562_sp + 1):
        n = counts["K562"] + i
        _, slot, pos, ct = _ctrl_vial("K562")
        ctrl_entries.append((slot, pos, "K562",
            f"{date}_{initials}_{lc_short}_{ms_short}_{k562_load}_K562_{n}_spare", False))
    for i in range(1, supermix_sp + 1):
        n = counts["Supermix"] + i
        _, slot, pos, ct = _ctrl_vial("Supermix")
        ctrl_entries.append((slot, pos, "Supermix",
            f"{date}_{initials}_{lc_short}_{ms_short}_{supermix_load}_Supermix_{n}_spare", False))
    for i in range(1, blank_sp + 1):
        n = counts["Blank"] + i
        _, slot, pos, ct = _ctrl_vial("Blank")
        ctrl_entries.append((slot, pos, "Blank",
            f"{date}_{initials}_{lc_short}_{ms_short}_Blank_{n}_spare", False))

    # XLSX
    wb = openpyxl.Workbook()
    ws = wb.active
    ws.append(QUEUE_COLS)
    for row in queue_rows:
        ws.append([row.get(col, "") for col in QUEUE_COLS])
    xbuf = io.BytesIO()
    wb.save(xbuf)

    # Sample plate PNG
    sample_png = plot_sample_plate(
        [{"well": s["well"], "name": s["name"],
          "group": group_assignments.get(s["name"], suggest_group(s["name"]))}
         for s in samples],
        f"{sample_slot} - Samples"
    )

    # Control plate PNGs (one per slot)
    ctrl_grids = {}
    ctrl_cmap  = {}
    ctrl_lmap  = {}
    for slot, pos, ctype, sid, in_queue in ctrl_entries:
        if slot not in ctrl_grids:
            ctrl_grids[slot] = {r: {c: "" for c in COLS} for r in ROWS}
        r, c = index_to_well(pos)
        ctrl_grids[slot][r][c] = sid
        ctrl_cmap[sid] = (CTRL_COLORS if in_queue else CTRL_COLORS_SPARE).get(ctype, "white")
        ctrl_lmap[sid] = f"{ctype}\n{sid.split('_')[-1]}"

    ctrl_pngs = {}
    for slot_num, grid in ctrl_grids.items():
        ctrl_pngs[slot_num] = plot_plate_png(grid, ctrl_cmap, f"Slot{slot_num} - Controls",
                                             label_map=ctrl_lmap)

    return {
        "queue_xlsx":  xbuf.getvalue(),
        "queue_rows":  queue_rows,
        "sample_png":  sample_png,
        "ctrl_pngs":   ctrl_pngs,
        "n_queue":     len(queue_rows),
        "counts":      counts,
        "spares":      {"K562": k562_sp, "Supermix": supermix_sp, "Blank": blank_sp},
    }


def build_zip(res, stem) -> bytes:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr(f"{stem}_queue.xlsx", res["queue_xlsx"])
        z.writestr(f"{stem}_samples.png", res["sample_png"])
        for slot_num, png in res["ctrl_pngs"].items():
            z.writestr(f"{stem}_slot{slot_num}.png", png)
    return buf.getvalue()


# ============================================================
# APP
# ============================================================
st.title("L2 MS Queue")
st.caption("Standalone MS queue generator. Input: CSV with `sample_name` and `well` columns.")

for key in ("msq_results", "msq_zip", "msq_csv_hash", "msq_group_assignments"):
    if key not in st.session_state:
        st.session_state[key] = None

date = datetime.date.today().strftime("%Y%m%d")

# Parameters
st.subheader("Instrument Parameters")
c1, c2, c3 = st.columns([1, 2, 2])
initials = c1.text_input("Initials", value="FK", key="initials")

lc_sel   = c2.selectbox("LC method", LC_OPTIONS, key="lc_sel")
lc_short = c2.text_input("Custom LC short name", value="", key="lc_custom") \
           if lc_sel == "Custom" else lc_sel

ms_sel   = c3.selectbox("MS method", MS_OPTIONS, key="ms_sel")
ms_short = c3.text_input("Custom MS short name", value="", key="ms_custom") \
           if ms_sel == "Custom" else ms_sel

c4, c5, c6, c7, c8 = st.columns(5)
sample_load   = c4.text_input("Sample load",  value="20ul", key="sload")
use_k562      = c5.checkbox("K562",           value=True,   key="k562")
k562_load     = c6.text_input("K562 load",    value="1ng",  key="k562l") if use_k562     else ""
use_supermix  = c7.checkbox("Supermix",       value=True,   key="smix")
supermix_load = c8.text_input("Supermix load",value="20ng", key="smixl") if use_supermix else ""

with st.expander("Instrument method paths", expanded=False):
    default_sep            = LC_METHODS.get(lc_short, "")
    default_ms, default_proc = MS_METHODS.get(ms_short, ("", ""))
    sep_method  = st.text_area("Separation Method", value=default_sep,    height=60, key="sep")
    inj_method  = st.text_input("Injection Method",  value="Standard",               key="inj")
    ms_method   = st.text_area("MS Method",          value=default_ms,    height=60, key="ms")
    proc_method = st.text_area("Processing Method",  value=default_proc,  height=60, key="proc")

_year      = date[:4]
_month_str = datetime.datetime.strptime(date, "%Y%m%d").strftime("%m %B")
with st.expander("Data paths", expanded=True):
    cp1, cp2    = st.columns(2)
    sample_path = cp1.text_input("Sample / control path",
                                 value=rf"D:\Data\{_year}\{_month_str}\Sample\{initials}",
                                 key="spath")
    blank_path  = cp2.text_input("Blank path",
                                 value=rf"D:\Data\{_year}\{_month_str}\Blank",
                                 key="bpath")

SLOT_OPTIONS = [f"Slot{i}" for i in range(1, 13)]
with st.expander("Slot assignment", expanded=True):
    sa1, sa2 = st.columns(2)
    sample_slot     = sa1.selectbox("Sample slot",        SLOT_OPTIONS, index=0, key="s_slot")
    ctrl_slot_start = int(sa2.selectbox("Controls start slot", SLOT_OPTIONS, index=1,
                                        key="c_slot").replace("Slot", ""))

st.divider()

# CSV upload
st.subheader("Sample List")
st.caption("Upload a CSV with at minimum columns `sample_name` and `well` (e.g. A1, H12)."
           " Optionally include a `group` column.")

uploaded = st.file_uploader("Upload CSV", type=["csv"], key="csv_upload")

if not uploaded:
    st.info("Upload a sample CSV to continue.")
    st.stop()

csv_bytes = uploaded.getvalue()
stem      = uploaded.name.replace(".csv", "")

csv_hash = hash(csv_bytes)
if st.session_state.msq_csv_hash != csv_hash:
    st.session_state.msq_group_assignments = None
    st.session_state.msq_results           = None
    st.session_state.msq_csv_hash          = csv_hash

df = pd.read_csv(io.BytesIO(csv_bytes))
df.columns = [c.strip().lower() for c in df.columns]

if "sample_name" not in df.columns or "well" not in df.columns:
    st.error("CSV must have `sample_name` and `well` columns.")
    st.stop()

samples = []
for _, row in df.iterrows():
    name = str(row["sample_name"]).strip()
    well = str(row["well"]).strip()
    grp  = str(row.get("group", suggest_group(name))).strip() if "group" in df.columns \
           else suggest_group(name)
    if name and well:
        samples.append({"name": name, "well": well, "group": grp})

samples.sort(key=lambda s: well_sort_key(s["well"]))

st.success(f"{len(samples)} samples loaded from `{uploaded.name}`.")

# Group assignment editor
st.subheader("Group Assignments")
st.caption("Groups auto-suggested from sample names. Edit **Group** to reassign.")

confirmed = st.session_state.msq_group_assignments or {}
init_data = [
    {"Sample": s["name"], "Well": s["well"],
     "Group": confirmed.get(s["name"]) or s["group"]}
    for s in samples
]

group_counts = {}
for row in init_data:
    g = row["Group"]
    group_counts[g] = group_counts.get(g, 0) + 1
gcols = st.columns(max(len(group_counts), 1))
for i, (g, n) in enumerate(sorted(group_counts.items())):
    gcols[i % len(gcols)].metric(g, f"{n} ROIs")

# Divide helper
_d1, _d2, _d3 = st.columns([2, 1, 1])
split_grp = _d1.selectbox("Divide group", sorted(group_counts.keys()), key="split_grp")
split_n   = _d2.number_input("Into N parts", min_value=2, max_value=20, value=2, step=1,
                              key="split_n")
if _d3.button("Divide", key="split_btn", use_container_width=True):
    rois = [r["Sample"] for r in init_data if r["Group"] == split_grp]
    chunk = math.ceil(len(rois) / int(split_n))
    updated = {r["Sample"]: r["Group"] for r in init_data}
    for i, roi in enumerate(rois):
        updated[roi] = f"{split_grp}{chr(ord('a') + i // chunk)}"
    st.session_state.msq_group_assignments = updated
    st.session_state.msq_results           = None

edited = st.data_editor(
    pd.DataFrame(init_data),
    column_config={
        "Sample": st.column_config.TextColumn("Sample", disabled=True),
        "Well":   st.column_config.TextColumn("Well",   disabled=True),
        "Group":  st.column_config.TextColumn("Group"),
    },
    hide_index=True, use_container_width=True, key="grp_editor"
)
new_assignments = dict(zip(edited["Sample"], edited["Group"]))

if st.button("Confirm grouping", type="primary", key="confirm"):
    st.session_state.msq_group_assignments = new_assignments
    st.session_state.msq_results           = None
    st.rerun()

if not st.session_state.msq_group_assignments:
    st.info("Confirm grouping above to continue.")
    st.stop()

group_assignments = st.session_state.msq_group_assignments

st.divider()
st.subheader("Run Order")
_ro1, _ro2, _ro3 = st.columns(3)
randomize_order = _ro1.checkbox("Randomize sample order within groups", value=False,
                                key="randomize")
block_size = _ro2.number_input("Blank interval (samples per block)", min_value=1,
                               max_value=50, value=GROUP_SIZE, step=1, key="block_size")
run_seed = _ro3.number_input("Random seed", min_value=0, max_value=9999, value=42,
                             step=1, key="run_seed", disabled=not randomize_order)

st.divider()
if st.button("Generate queue", type="primary", key="gen"):
    params = dict(
        date=date, initials=initials, lc_short=lc_short, ms_short=ms_short,
        sample_load=sample_load, k562_load=k562_load, supermix_load=supermix_load,
        use_k562=use_k562, use_supermix=use_supermix,
        sep_method=sep_method, inj_method=inj_method,
        ms_method=ms_method, proc_method=proc_method,
        sample_path=sample_path, blank_path=blank_path,
        sample_slot=sample_slot, ctrl_slot_start=ctrl_slot_start,
        randomize_order=randomize_order, block_size=block_size, run_seed=run_seed,
    )
    with st.spinner("Generating..."):
        res = build_queue(samples, group_assignments, params)
        st.session_state.msq_results = res
        st.session_state.msq_zip     = build_zip(res, stem)

res = st.session_state.msq_results
if not res:
    st.stop()

c = res["counts"]
sp = res["spares"]
st.success(
    f"{res['n_queue']} queue rows | "
    f"K562: {c['K562']} (+{sp['K562']} spare) | "
    f"Supermix: {c['Supermix']} (+{sp['Supermix']} spare) | "
    f"Blank: {c['Blank']} (+{sp['Blank']} spare)"
)

# Plate maps
img_cols = st.columns(1 + len(res["ctrl_pngs"]))
with img_cols[0]:
    st.subheader("Sample plate")
    st.image(res["sample_png"], use_container_width=True)
for i, (slot_num, png) in enumerate(res["ctrl_pngs"].items()):
    with img_cols[i + 1]:
        st.subheader(f"Slot{slot_num} Controls")
        st.image(png, use_container_width=True)

st.subheader("Downloads")
dl1, dl2 = st.columns(2)
dl1.download_button("Download queue (zip)", st.session_state.msq_zip,
                    file_name=f"{stem}_ms_queue.zip", mime="application/zip", type="primary")
dl2.download_button("Download queue XLSX", res["queue_xlsx"],
                    file_name=f"{stem}_queue.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
