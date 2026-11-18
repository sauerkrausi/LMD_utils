"""
ms_queue_utils.py
=================
Tab 4: MS Sample Queue Generator.

Ported from create_ms_queue_streamlit.py.
Accepts sample list CSV piped from Tab 3 or via direct upload.

Instruments: Bruker timsTOF Ultra (HyStar XLSX) and Thermo Astral Zoom
(Xcalibur CSV) with Evosep One or Vanquish Neo.
Outputs: queue file + slot plate CSVs + plate map PNGs + run order CSV (zip).
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
PLATE_FORMATS = {
    "96":  (list("ABCDEFGH"),         list(range(1, 13))),
    "384": (list("ABCDEFGHIJKLMNOP"), list(range(1, 25))),
}
WELL_RE = re.compile(r"^\s*([A-Za-z])\s*0*(\d+)\s*$")
GROUP_SIZE = 6
SUPERMIX_INTERVAL = 12

QUEUE_COLS = [
    "Vial", "Sample ID", "Method Set", "Separation Method",
    "Injection Method", "MS Method", "Processing Method",
    "Sample Type", "Volume [µl]", "Data Path", "Run Automated Processing",
]

# ============================================================
# INSTRUMENTS
# ============================================================
INST_ULTRA = "timsTOF Ultra"
INST_ZOOM  = "Astral Zoom"
LC_EVOSEP  = "Evosep One"
LC_NEO     = "Vanquish Neo"

# Xcalibur sequence columns, names must match exactly
XCAL_COLS = [
    "Sample Type", "File Name", "Sample ID", "Path", "Instrument Method",
    "Process Method", "Calibration File", "Position", "Inj Vol", "Level",
    "Sample Wt", "Sample Vol", "ISTD Amt", "Dil Factor", "L1 Study",
    "L2 Client", "L3 Laboratory", "L4 Company", "L5 Phone", "Comment", "Sample Name",
]
XCAL_HEADER = "Bracket Type=4"

# Dummy presets from methods.csv, replace once the real methods are set
ZOOM_MS_SHORT     = "Zoom"
ZOOM_LC_NAME      = {LC_EVOSEP: "Evosep", LC_NEO: "VanquishNeo"}
ZOOM_GRADIENT     = {LC_EVOSEP: "40SPD",  LC_NEO: "60SPD"}
ZOOM_MS_METHODS = {
    "defaultDIA": {
        "noFAIMS": r"D:\MS method\22min_Thermo_60SPDgrad_DIA_MS1OT_MS2AST_20251125",
        "FAIMS":   r"D:\MS method\ROI opt\22min_adj5_grad_FAIMS_-50CV_DIA_MS1OT_MS2AST_20261102",
    },
}

# Accepted spellings of the sample name column, first match wins
SAMPLE_COLS = ["ROI", "Sample", "sample_name", "Sample Name", "Sample_Name", "sample"]

EVOSEP_SLOTS  = [f"Slot{i}" for i in range(1, 7)]
NEO_TRAYS     = ["R", "G", "B", "Y"]
NEO_FIXED_DEFAULT = {"Blank": "R:F9", "K562": "R:A1", "Supermix": "R:A2"}
NEO_POS_RE    = re.compile(r"^[RGBY]:[A-H](1[0-2]|[1-9])$")

CTRL_COLORS       = {"K562": "#4C9BE8", "Supermix": "#F4A261", "Blank": "#B7E4C7"}
CTRL_COLORS_SPARE = {"K562": "#C5DDF7", "Supermix": "#FAE0C8", "Blank": "#E6F7EC"}
GREY_REPLACEMENTS = {14: "#777777", 15: "#555555"}


# ============================================================
# HELPERS
# ============================================================
def well_to_slot1(well_id, n_cols=12):
    row = ord(well_id[0].upper()) - ord('A')
    col = int(well_id[1:])
    return f"Slot1:{row * n_cols + col}"


def check_well(raw, rows, cols):
    """Normalized well or None plus the reason."""
    m = WELL_RE.match("" if raw is None else str(raw))
    if not m:
        return None, "not a well id"
    r, c = m.group(1).upper(), int(m.group(2))
    if r not in rows:
        return None, f"row {r} outside plate"
    if c not in cols:
        return None, f"column {c} outside plate"
    return f"{r}{c}", None


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


def normalize_well(well_id):
    """A01 to A1."""
    w = well_id.strip().upper()
    return f"{w[0]}{int(w[1:])}"


def faims_token(on, cv):
    """FAIMS plus optional CV(s), - written as neg, e.g. -40,-60 to FAIMS_neg40_neg60."""
    if not on:
        return ""
    parts = [c.strip() for c in re.split(r"[,;\s/]+", cv or "") if c.strip()]
    cvs   = "_".join(re.sub(r"[^A-Za-z0-9]", "", c.replace("-", "neg").replace("+", "pos")) for c in parts)
    return "FAIMS" + (f"_{cvs}" if cvs else "")


def name_prefix(p):
    """Ultra: date_initials_LC_MS. Zoom: date_initials_MS_LC_FAIMS_Acq_Gradient. Tags last, empties skipped."""
    if p.get("inst", "ultra") == "ultra":
        parts = [p["date"], p["initials"], p["lc_short"], p["ms_short"]]
    else:
        z = p.get("zoom_name", {})
        parts = [p["date"], p["initials"], z.get("ms", ""), z.get("lc", ""),
                 faims_token(z.get("faims"), z.get("cv")), z.get("acq", ""), z.get("gradient", "")]
        parts = [x for x in parts if str(x).strip()]
    if p.get("tags"):
        parts.append(p["tags"])
    return "_".join(str(x).strip() for x in parts)


# Column aliases, matched case insensitively, so the shared templates load here and in msQ
COL_ALIASES = {
    "Well_ID": ["well_id", "well id", "wellid", "well", "location"],
    "Group":   ["group"],
    "Load":    ["load"],
    "Inj Vol": ["inj vol", "inj_vol", "injvol", "inj vol (µl)", "inj_vol_ul"],
    "Plate":   ["plate"],
}


def _key(h):
    return re.sub(r"\s+", " ", h.strip().lower())


def read_sample_rows(csv_bytes):
    """Parse the sample list, map the sample name column to ROI and known aliases to their
    canonical names. Returns rows, fieldnames, source column."""
    text   = csv_bytes.decode("utf-8-sig")
    reader = csv.DictReader(io.StringIO(text))
    fields = [h.strip() for h in (reader.fieldnames or [])]
    reader.fieldnames = fields
    rows   = list(reader)
    by_key = {_key(f): f for f in fields}
    for canon, aliases in COL_ALIASES.items():
        if canon in fields:
            continue
        src_f = next((by_key[a] for a in aliases if a in by_key), None)
        if src_f:
            for r in rows:
                r[canon] = r.get(src_f, "")
            fields = fields + [canon]
    by_key = {_key(f): f for f in fields}
    src    = next((by_key[_key(c)] for c in SAMPLE_COLS if _key(c) in by_key), None)

    # CellenONE export: Location plus Event Id, name becomes Location_EventId
    ev_col = next((by_key[a] for a in ("event id", "event_id", "eventid") if a in by_key), None)
    if not src and ev_col and "Well_ID" in fields:
        for r in rows:
            ev  = str(r.get(ev_col, "")).strip()
            ev  = ev[:-2] if ev.endswith(".0") else ev
            w   = r.get("Well_ID", "").strip()
            w   = normalize_well(w) if WELL_RE.match(w) else w
            r["ROI"] = f"{w}_{ev}"
        return rows, fields + ["ROI"], f"{ev_col} (CellenONE)"
    if src and src != "ROI":
        for r in rows:
            r["ROI"] = r.get(src, "")
        fields = fields + ["ROI"]
    return rows, fields, src


EXAMPLE_CSV = b"sample_name,well\nSample_1,A1\nSample_2,A2\nSample_3,B1\n"


def find_templates():
    """Example CSVs: any CSV in a templates folder (any case) at the repo root or in an app
    folder. Falls back to a built in example so a download is always offered."""
    from pathlib import Path
    root = Path(__file__).resolve().parents[2]
    dirs = [d for d in [root, *root.iterdir()] if d.is_dir()]
    found = {}
    for base in dirs:
        for d in base.iterdir():
            if d.is_dir() and d.name.lower() in ("templates", "template"):
                for f in sorted(d.glob("*.csv")) + sorted(d.glob("*.CSV")):
                    found.setdefault(f.name, f.read_bytes())
    return found or {"sample_list_example.csv": EXAMPLE_CSV}


def apply_sample_table(rows, table):
    """Per-sample Load and Inj Vol from the editor. A sample listed again becomes a repeat
    injection from the same well, attached to its original so it runs right after it."""
    if not table:
        return rows
    by_name = {}
    for r in rows:
        by_name.setdefault(r["ROI"].strip(), r)
        r["_repeats"] = []
    seen = set()
    for t in table:
        name = str(t.get("Sample") or "").strip()
        orig = by_name.get(name)
        if orig is None:
            continue
        load = str(t.get("Load") or "").strip()
        vol  = t.get("Inj Vol")
        vol  = None if vol is None or (isinstance(vol, float) and math.isnan(vol)) else float(vol)
        if name not in seen:
            seen.add(name)
            orig["_load"], orig["_inj_vol"] = load, vol
        else:
            rr = dict(orig, _load=load, _inj_vol=vol, _repeats=[])
            rr["_rep"] = len(orig["_repeats"]) + 2
            orig["_repeats"].append(rr)
    return rows


def _num(v):
    return int(v) if float(v).is_integer() else v


def sample_vial(well_id, sample_slot, inst, n_cols=12):
    """Ultra and Evosep: SlotN:pos, row wise. Vanquish: tray:well."""
    if inst == "zoom_neo":
        return f"{sample_slot}:{normalize_well(well_id)}"
    return well_to_slot1(normalize_well(well_id), n_cols).replace("Slot1", sample_slot)


def _xcal_cell(v):
    s = "" if v is None else str(v)
    return f'"{s}"' if ("," in s or '"' in s) else s


def build_xcalibur_csv(queue_rows, inj_vol):
    """Xcalibur sequence: header line, 21 columns, trailing comma, CRLF."""
    lines = [XCAL_HEADER + ",", ",".join(XCAL_COLS) + ","]
    for r in queue_rows:
        rec = {
            "Sample Type": "Unknown", "File Name": r["Sample ID"], "Sample ID": 1,
            "Path": r["Data Path"], "Instrument Method": r["MS Method"],
            "Position": r["Vial"],
            "Inj Vol": _num(r["_inj_vol"]) if r.get("_inj_vol") is not None else inj_vol.get(r.get("_type", "Sample"), 1),
            "Sample Wt": 0, "Sample Vol": 0, "ISTD Amt": 0, "Dil Factor": 1,
            "Comment": r.get("_comment", ""),
        }
        lines.append(",".join(_xcal_cell(rec.get(c, "")) for c in XCAL_COLS) + ",")
    return ("\r\n".join(lines) + "\r\n").encode("utf-8")


def make_row(vial, sample_id, data_path, sep_method, inj_method, ms_method, proc_method,
             ctype="Sample"):
    return {
        "_type": ctype,
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
               label_color_map=None, rows=None, cols=None) -> bytes:
    rows = rows or ROWS
    cols = cols or COLS
    n_r, n_c = len(rows), len(cols)
    small = n_c > 12
    fig, ax = plt.subplots(figsize=(14, 8) if not small else (18, 12))
    ax.set_xlim(-0.5, n_c + 0.5)
    ax.set_ylim(-0.5, n_r + 0.5)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)

    for r_idx, r in enumerate(rows):
        for c_idx, c in enumerate(cols):
            x, y  = c_idx, (n_r - 1) - r_idx
            label = grid[r][c]
            color = color_map.get(label, "white") if label else "whitesmoke"
            edge  = "#aaaaaa" if not label else "#444444"
            ax.add_patch(plt.Circle((x, y), 0.42, color=color, ec=edge, lw=0.8, zorder=2))
            display    = (label_map or {}).get(label, label)
            font_color = (label_color_map or {}).get(label, "black")
            if display:
                fsize = (4.5 if len(display) > 10 else 5.5) * (0.6 if small else 1)
                ax.text(x, y, display, ha="center", va="center",
                        fontsize=fsize, zorder=3, color=font_color, clip_on=True)

    for r_idx, r in enumerate(rows):
        ax.text(-0.55, (n_r - 1) - r_idx, r, ha="right", va="center", fontsize=9, fontweight="bold")
    for c_idx, c in enumerate(cols):
        ax.text(c_idx, n_r, str(c), ha="center", va="bottom", fontsize=9, fontweight="bold")

    seen = {}
    for r in rows:
        for c in cols:
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

    # Too many entries (e.g. one per cell) would cover the plate
    if patches and len(patches) <= 30:
        if small:
            fig.legend(handles=patches, loc="upper center", bbox_to_anchor=(0.5, 0.0),
                       ncol=min(len(patches), 8), fontsize=7, framealpha=0.9, title="Legend", title_fontsize=8)
        else:
            fig.legend(handles=patches, loc="lower center",
                       bbox_to_anchor=(0.5, -0.02), ncol=min(len(patches), 6),
                       fontsize=7, framealpha=0.9, title="Legend", title_fontsize=8)

    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


# ============================================================
# GROUPING / BLOCKING
# ============================================================
def suggest_group(name: str) -> str:
    """Strip trailing numbers/spaces to get a run-group key."""
    key = re.sub(r'[\s_-]*\d+\s*$', '', name.strip()).strip()
    return key or name.strip()


def suggest_block(name: str) -> str:
    """First token before - or _ is the block key (e.g. 'H20' from 'H20-015340-A3-Core32-A')."""
    parts = re.split(r'[-_]', name.strip(), maxsplit=1)
    return parts[0].strip() if parts and parts[0].strip() else name.strip()


# ============================================================
# RUN SEQUENCE
# ============================================================
def collapse_repeats(events):
    """Drop a control injection that immediately repeats the previous one."""
    out = []
    for ev in events:
        if ev[0] != "Sample" and out and out[-1][0] == ev[0]:
            continue
        out.append(ev)
    return out


def apply_overrides(events, overrides):
    """Manual add or remove of one control, anchored after a named ROI.

    overrides: [{"action": "add"|"remove", "type": "K562"|"Supermix"|"Blank",
                 "after": roi_name}]
    """
    if not overrides:
        return events
    adds, drops = {}, {}
    for ov in overrides:
        roi = (ov.get("after") or "").strip()
        ctl = ov.get("type")
        if not roi or ctl not in ("K562", "Supermix", "Blank"):
            continue
        target = drops if ov.get("action") == "remove" else adds
        target.setdefault(roi, []).append(ctl)

    out, pending = [], []
    for kind, row in events:
        if kind == "Sample":
            out.append((kind, row))
            roi = row["ROI"].strip()
            out.extend((c, None) for c in adds.get(roi, []))
            pending = list(drops.get(roi, []))
            continue
        if pending and kind in pending:
            pending.remove(kind)
            continue
        out.append((kind, row))
    return out


def control_gaps(events):
    """Max samples between Blanks and between Supermix injections."""
    max_b = max_s = since_b = since_s = 0
    for kind, _ in events:
        if kind == "Sample":
            since_b += 1
            since_s += 1
            max_b, max_s = max(max_b, since_b), max(max_s, since_s)
        elif kind == "Blank":
            since_b = 0
        elif kind == "Supermix":
            since_s = 0
    return max_b, max_s


def build_run_events(groups_seen, group_map, p, block_assignments=None, closing_set=False):
    """Ordered run events as (kind, row), kind in Sample, K562, Supermix, Blank.

    Block start always fires the full K562 + Supermix + Blank set and resets both
    counters, so no interval carries across a block boundary. Inside a block the samples
    of all its groups are batched together, evenly, at most block_size per batch, so a
    group can span a Blank. A Blank fires between batches, a Supermix when the next batch
    would push past smix_interval and is always followed by a Blank. Only samples count.
    """
    use_k562     = p["use_k562"]
    use_supermix = p["use_supermix"]
    block_size   = int(p.get("block_size", GROUP_SIZE))
    smix_every   = int(p.get("smix_interval", SUPERMIX_INTERVAL))
    k562_w_smix  = bool(p.get("k562_with_smix")) and use_k562
    std_on       = use_supermix or k562_w_smix

    events = []
    since_smix = since_blank = 0

    def push_controls(n_next):
        """Controls due before the next batch of n_next samples."""
        nonlocal since_smix, since_blank
        if not n_next:
            return
        if std_on and since_smix + n_next > smix_every:
            if k562_w_smix:
                events.append(("K562", None))
            if use_supermix:
                events.append(("Supermix", None))
            events.append(("Blank", None))
            since_smix = since_blank = 0
        elif since_blank + n_next > block_size:
            events.append(("Blank", None))
            since_blank = 0

    # Samples of a block run together, in group order
    block_order, block_samples = [], {}
    for grp in groups_seen:
        blk = (block_assignments or {}).get(grp, suggest_block(grp))
        if blk not in block_order:
            block_order.append(blk)
            block_samples[blk] = []
        block_samples[blk].extend(group_map[grp])

    for bk, blk in enumerate(block_order):
        if use_k562:
            events.append(("K562", None))
        if use_supermix:
            events.append(("Supermix", None))
        events.append(("Blank", None))
        since_smix = since_blank = 0

        batches = split_groups(block_samples[blk], max_size=block_size)
        for bi, batch in enumerate(batches):
            events.extend(("Sample", row) for row in batch)
            since_smix  += len(batch)
            since_blank += len(batch)
            if closing_set and bk == len(block_order) - 1 and bi == len(batches) - 1:
                continue
            # Nothing after the last batch of a block, the next block start fires its set
            push_controls(len(batches[bi + 1]) if bi + 1 < len(batches) else 0)

    if closing_set:
        if use_k562:
            events.append(("K562", None))
        if use_supermix:
            events.append(("Supermix", None))
        events.append(("Blank", None))

    return collapse_repeats(apply_overrides(events, p.get("overrides")))


# ============================================================
# CONTROL COUNTING (from run sequence)
# ============================================================
def count_controls(events, use_k562, use_supermix):
    """Control totals plus spares, counted from the built run sequence."""
    k562_used     = sum(1 for k, _ in events if k == "K562")
    supermix_used = sum(1 for k, _ in events if k == "Supermix")
    blank_used    = sum(1 for k, _ in events if k == "Blank")

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
# CONTROL ALLOCATOR: row-banded, dynamic multi-slot
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


class FixedAllocator:
    """Vanquish Neo: every control of one type is drawn from the same vial."""
    def __init__(self, positions: dict):
        self.positions = positions
        self.entries   = []

    def add(self, ctype: str, sid: str, in_queue: bool = True) -> str:
        return self.positions[ctype]


# ============================================================
# QUEUE BUILDER
# ============================================================
def _build_plate_queue(rows_for_plate, group_assignments, p,
                       sample_slot, ctrl_start_slot,
                       counts, alloc, block_assignments=None, closing_set=False):
    """Build queue rows for a single sample plate, mutating counts and alloc in place."""
    use_k562      = p["use_k562"]
    use_supermix  = p["use_supermix"]
    date          = p["date"]
    initials      = p["initials"]
    lc_short      = p["lc_short"]
    ms_short      = p["ms_short"]
    sample_load   = p["sample_load"]
    k562_load     = p.get("k562_load", "")
    supermix_load = p.get("supermix_load", "")
    sep_method    = p["sep_method"]
    inj_method    = p["inj_method"]
    ms_method     = p["ms_method"]
    proc_method   = p["proc_method"]
    sample_path   = p["sample_path"]
    blank_path    = p["blank_path"]

    std_path        = p.get("std_path") or sample_path
    inst            = p.get("inst", "ultra")
    prefix          = name_prefix(p)
    rand_mode       = p.get("randomize_mode") or ("Within group" if p.get("randomize_order") else "Off")
    rng_seed        = int(p.get("run_seed", 42))

    samples = [r for r in rows_for_plate
               if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    groups_seen, group_map, roi_group = [], {}, {}
    for row in samples:
        roi = row["ROI"].strip()
        grp = group_assignments.get(roi, suggest_group(roi))
        roi_group[roi] = grp
        if grp not in group_map:
            group_map[grp] = []
            groups_seen.append(grp)
        group_map[grp].append(row)

    orig_blocks = dict(block_assignments or {})

    # Randomization, each sample keeps its own well
    import random as _random
    rng = _random.Random(rng_seed)
    if rand_mode == "Within group":
        for grp in groups_seen:
            rng.shuffle(group_map[grp])
    elif rand_mode == "Whole plate" and samples:
        flat = list(samples)
        rng.shuffle(flat)
        groups_seen, group_map = ["All samples"], {"All samples": flat}
        block_assignments = {"All samples": "All samples"}

    events = build_run_events(groups_seen, group_map, p, block_assignments, closing_set)
    # Repeats run directly after their original, controls placed after the original follow the repeats
    expanded = []
    for kind, row in events:
        expanded.append((kind, row))
        if kind == "Sample":
            expanded.extend(("Sample", rr) for rr in row.get("_repeats", []))
    events = expanded
    ctrl_counts = count_controls(events, use_k562, use_supermix)
    fixed_ctrl  = inst == "zoom_neo"
    if fixed_ctrl:
        alloc_pl = FixedAllocator(p["fixed_positions"])
    else:
        alloc_pl = ControlAllocator(ctrl_counts, use_k562, use_supermix, start_slot=ctrl_start_slot)
    queue_pl = []

    def add_ctrl(ctype, load, path):
        counts[ctype] += 1
        sid  = f"{prefix}_{load}_{ctype}_{counts[ctype]}" if load else f"{prefix}_{ctype}_{counts[ctype]}"
        vial = alloc_pl.add(ctype, sid)
        row  = make_row(vial, sid, path, sep_method, inj_method, ms_method, proc_method, ctype=ctype)
        queue_pl.append(row)

    for kind, row in events:
        if kind == "Sample":
            roi      = row["ROI"].strip()
            well_pos = sample_vial(row["Well_ID"], sample_slot, inst,
                                   len(PLATE_FORMATS.get(str(p.get("plate_format", "96")), (ROWS, COLS))[1]))
            load     = row.get("_load") or sample_load
            sid      = "_".join(x for x in (prefix, load, roi) if x) + (f"_r{row['_rep']}" if row.get("_rep") else "")
            qrow     = make_row(well_pos, sid, sample_path, sep_method, inj_method, ms_method, proc_method)
            if row.get("_inj_vol") is not None:
                qrow["_inj_vol"]     = row["_inj_vol"]
                qrow["Volume [µl]"]  = _num(row["_inj_vol"])
            qrow.update({"_well": normalize_well(row["Well_ID"]), "_group": roi_group.get(roi, ""),
                         "_block": orig_blocks.get(roi_group.get(roi, ""), "")})
            queue_pl.append(qrow)
        elif kind == "K562":
            add_ctrl("K562", k562_load, std_path)
        elif kind == "Supermix":
            add_ctrl("Supermix", supermix_load, std_path)
        else:
            add_ctrl("Blank", "", blank_path)

    # Spares (not in queue, just placed in ctrl plate), none for fixed vials
    if fixed_ctrl:
        k562_sp = smix_sp = blnk_sp = 0
    else:
        k562_sp = ctrl_counts["k562_spares"]
        smix_sp = ctrl_counts["supermix_spares"]
        blnk_sp = ctrl_counts["blank_spares"]
    for ctype, load, n_sp in (("K562", k562_load, k562_sp), ("Supermix", supermix_load, smix_sp), ("Blank", "", blnk_sp)):
        for i in range(1, n_sp + 1):
            n   = counts[ctype] + i
            sid = f"{prefix}_{load}_{ctype}_{n}_spare" if load else f"{prefix}_{ctype}_{n}_spare"
            alloc_pl.add(ctype, sid, in_queue=False)

    alloc.entries.extend(alloc_pl.entries)

    max_blank_gap, max_smix_gap = control_gaps(events)

    return {
        "queue_rows":    queue_pl,
        "max_blank_gap": max_blank_gap,
        "max_smix_gap":  max_smix_gap,
        "groups_seen":   groups_seen,
        "samples":       samples,
        "ctrl_counts":   ctrl_counts,
        "alloc":         alloc_pl,
        "sample_slot":   sample_slot,
        "k562_spares":   k562_sp,
        "supermix_spares": smix_sp,
        "blank_spares":  blnk_sp,
    }


def build_queue_core(csv_bytes: bytes, group_assignments: dict, p: dict,
                     block_assignments: dict = None) -> dict:
    """
    group_assignments: {roi_name: group_label}
    plate_slot_map in p: {plate_label: {"sample_slot": "Slot1", "ctrl_start_slot": 2}}
    If CSV has no Plate column, treated as single plate.
    """
    date          = p["date"]
    initials      = p["initials"]
    lc_short      = p["lc_short"]
    ms_short      = p["ms_short"]
    use_k562      = p["use_k562"]
    use_supermix  = p["use_supermix"]
    stem          = p["stem"]
    plate_slot_map = p.get("plate_slot_map", {})
    # fallback globals (used when no plate_slot_map entry exists)
    default_sample_slot   = p.get("sample_slot", "Slot1")
    default_ctrl_start    = p.get("ctrl_start_slot", 2)

    all_rows, fields, _ = read_sample_rows(csv_bytes)
    all_rows  = apply_sample_table(all_rows, p.get("sample_table"))
    has_plate = "Plate" in fields

    # Split rows by plate
    plate_rows = {}
    for row in all_rows:
        pl = row.get("Plate", "Plate1").strip() if has_plate else "Plate1"
        plate_rows.setdefault(pl, []).append(row)
    plate_order = sorted(plate_rows.keys())

    # If no plate_slot_map provided, auto-assign: Plate1->Slot1+Slot2, Plate2->Slot3+Slot4
    if not plate_slot_map:
        for i, pl in enumerate(plate_order):
            plate_slot_map[pl] = {
                "sample_slot":    f"Slot{i * 2 + 1}",
                "ctrl_start_slot": i * 2 + 2,
            }

    # Process plates
    all_queue       = []
    all_groups_seen = []
    counts          = {"K562": 0, "Supermix": 0, "Blank": 0}
    combined_alloc  = ControlAllocator.__new__(ControlAllocator)
    combined_alloc.entries = []

    plate_results  = {}
    for pl in plate_order:
        slot_info = plate_slot_map.get(pl, {
            "sample_slot":    default_sample_slot,
            "ctrl_start_slot": default_ctrl_start,
        })
        pr = _build_plate_queue(
            plate_rows[pl], group_assignments, p,
            slot_info["sample_slot"], slot_info["ctrl_start_slot"],
            counts, combined_alloc, block_assignments=block_assignments,
            closing_set=(pl == plate_order[-1]),
        )
        plate_results[pl] = pr
        all_queue.extend(pr["queue_rows"])
        all_groups_seen.extend(pr["groups_seen"])

    # Run tracking, written to the Xcalibur Comment column and run_order.csv
    inst      = p.get("inst", "ultra")
    rand_mode = p.get("randomize_mode") or ("Within group" if p.get("randomize_order") else "Off")
    seed_tag  = f"_seed{int(p.get('run_seed', 42))}" if rand_mode != "Off" else ""
    run_order = []
    for n, r in enumerate(all_queue, start=1):
        well = r.get("_well", "")
        r["_comment"] = f"run{n}_well{well}{seed_tag}" if well else f"run{n}"
        run_order.append({"Run": n, "Sample ID": r["Sample ID"], "Type": r.get("_type", ""),
                          "Position": r["Vial"], "Well": well, "Group": r.get("_group", ""),
                          "Block": r.get("_block", ""), "Randomization": rand_mode,
                          "Seed": int(p.get("run_seed", 42)) if rand_mode != "Off" else ""})
    ro_buf = io.StringIO()
    pd.DataFrame(run_order).to_csv(ro_buf, index=False)
    run_order_csv = ro_buf.getvalue().encode("utf-8")

    # Queue file, HyStar XLSX for the Ultra, Xcalibur CSV for the Zoom
    queue_xlsx = queue_csv = None
    if inst == "ultra":
        try:
            import openpyxl
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.append(QUEUE_COLS)
            for row in all_queue:
                ws.append([row.get(col, "") for col in QUEUE_COLS])
            xlsx_buf = io.BytesIO()
            wb.save(xlsx_buf)
            queue_xlsx = xlsx_buf.getvalue()
        except ImportError:
            queue_xlsx = None
    else:
        queue_csv = build_xcalibur_csv(all_queue, p.get("inj_vol", {}))

    # Per-plate sample slot grids + PNGs
    sample_slot_outputs = {}   # {sample_slot_label: {grid, png, csv}}
    for pl in plate_order:
        pr        = plate_results[pl]
        sslot     = pr["sample_slot"]
        rows_pl   = plate_rows[pl]
        p_rows, p_cols = PLATE_FORMATS.get(str(p.get("plate_format", "96")), (ROWS, COLS))
        dropout_w = {normalize_well(r["Well_ID"]) for r in rows_pl
                     if r["Well_ID"].strip() and r.get("Dropout {Y/N}", "").strip().upper() == "Y"}
        grid = {r: {c: "" for c in p_cols} for r in p_rows}
        well_to_group = {}
        for row in rows_pl:
            w_id = row["Well_ID"].strip()
            if w_id:
                w_id = normalize_well(w_id)
                roi  = row["ROI"].strip()
                grid[w_id[0]][int(w_id[1:])] = roi
                well_to_group[roi] = group_assignments.get(roi, suggest_group(roi))

        groups_pl = pr["groups_seen"]
        if groups_pl == ["All samples"]:
            groups_pl = list(dict.fromkeys(well_to_group.values()))
        n_g       = max(len(groups_pl), 1)
        group_colors = {}
        for i, g in enumerate(groups_pl):
            idx = round(i / n_g * 20)
            group_colors[g] = GREY_REPLACEMENTS.get(idx, cm.tab20(i / n_g))

        color_map  = {}
        label_color = {}
        for r in p_rows:
            for c in p_cols:
                roi = grid[r][c]
                if roi:
                    grp = well_to_group.get(roi, "")
                    if f"{r}{c}" in dropout_w:
                        color_map[roi]   = "#dddddd"
                        label_color[roi] = "red"
                    else:
                        color_map[roi]   = group_colors.get(grp, "white")
                        label_color[roi] = "black"

        roi_to_group = {roi: well_to_group.get(roi, roi) for roi in color_map}
        png = plot_plate(grid, color_map, f"{sslot} - {pl} Samples ({stem})",
                         legend_group_map=roi_to_group, label_color_map=label_color,
                         rows=p_rows, cols=p_cols)
        buf = io.StringIO()
        wtr = csv.writer(buf)
        wtr.writerow([""] + p_cols)
        for r in p_rows:
            wtr.writerow([r] + [grid[r][c] for c in p_cols])
        key = f"{sslot}_{pl}"
        sample_slot_outputs[key] = {
            "grid": grid, "png": png,
            "csv":  buf.getvalue().encode("utf-8"),
            "sample_slot": sslot, "plate": pl,
        }

    # Control slot grids + PNGs (combined from all plates)
    ctrl_color_map = {}
    ctrl_label_map = {}
    for _, _, ctype, sid, in_queue in combined_alloc.entries:
        ctrl_color_map[sid] = (CTRL_COLORS if in_queue else CTRL_COLORS_SPARE).get(ctype, "white")
        num = sid.split("_")[-1]
        ctrl_label_map[sid] = f"{ctype}\n{num}"
    sid_to_type = {sid: ctype for _, _, ctype, sid, _ in combined_alloc.entries}

    # Rebuild plates dict from combined_alloc
    ctrl_grids = {}
    for slot, pos, ctype, sid, _ in combined_alloc.entries:
        if slot not in ctrl_grids:
            ctrl_grids[slot] = {r: {c: "" for c in COLS} for r in ROWS}
        r, c = index_to_well(pos)
        ctrl_grids[slot][r][c] = sid

    ctrl_slot_outputs = {}
    for slot_num, grid in ctrl_grids.items():
        png = plot_plate(grid, ctrl_color_map, f"Slot{slot_num} - Controls ({stem})",
                         label_map=ctrl_label_map, legend_group_map=sid_to_type)
        buf = io.StringIO()
        wtr = csv.writer(buf)
        wtr.writerow([""] + COLS)
        for r in ROWS:
            wtr.writerow([r] + [grid[r][c] for c in COLS])
        ctrl_slot_outputs[slot_num] = {"png": png, "csv": buf.getvalue().encode("utf-8")}

    total_k562_sp = sum(pr["k562_spares"]     for pr in plate_results.values())
    total_smix_sp = sum(pr["supermix_spares"] for pr in plate_results.values())
    total_blnk_sp = sum(pr["blank_spares"]    for pr in plate_results.values())

    max_ctrl_slot = max((slot for slot, *_ in combined_alloc.entries), default=0)

    return {
        "queue_xlsx":         queue_xlsx,
        "queue_csv":          queue_csv,
        "run_order_csv":      run_order_csv,
        "inst":               inst,
        "max_ctrl_slot":      max_ctrl_slot,
        "sample_slot_outputs": sample_slot_outputs,
        "ctrl_slots":         ctrl_slot_outputs,
        "queue_rows":         all_queue,
        "n_queue":            len(all_queue),
        "max_blank_gap":      max((r["max_blank_gap"] for r in plate_results.values()), default=0),
        "max_smix_gap":       max((r["max_smix_gap"]  for r in plate_results.values()), default=0),
        "counts":             counts,
        "k562_spares":        total_k562_sp,
        "supermix_spares":    total_smix_sp,
        "blank_spares":       total_blnk_sp,
        "stem":               stem,
        "groups":             all_groups_seen,
        "plate_order":        plate_order,
        "plate_slot_map":     plate_slot_map,
    }


def build_zip(res: dict) -> bytes:
    buf  = io.BytesIO()
    stem = res["stem"]
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        if res["queue_xlsx"]:
            z.writestr(f"{stem}_queue.xlsx", res["queue_xlsx"])
        if res.get("queue_csv"):
            z.writestr(f"{stem}_queue_xcalibur.csv", res["queue_csv"])
        if res.get("run_order_csv"):
            z.writestr(f"{stem}_run_order.csv", res["run_order_csv"])
        for key, data in res["sample_slot_outputs"].items():
            z.writestr(f"{stem}_{key}_samples.csv", data["csv"])
            z.writestr(f"{stem}_{key}_samples.png", data["png"])
        for slot_num, data in res["ctrl_slots"].items():
            z.writestr(f"{stem}_slot{slot_num}.csv", data["csv"])
            z.writestr(f"{stem}_slot{slot_num}.png", data["png"])
    return buf.getvalue()


# ============================================================
# STREAMLIT TAB
# ============================================================
# Run mode presets, every value stays editable afterwards
RUN_PRESETS = {
    "Bulk": {"msq_rand_mode": "Off", "msq_block_size": GROUP_SIZE, "msq_smix_interval": SUPERMIX_INTERVAL,
             "msq_k562_with_smix": False, "msq_plate_fmt": "96"},
    "SCP":  {"msq_rand_mode": "Whole plate", "msq_block_size": 12, "msq_smix_interval": 24,
             "msq_k562_with_smix": True, "msq_plate_fmt": "384"},
}


def _apply_run_preset():
    for k, v in RUN_PRESETS[st.session_state.msq_run_mode].items():
        st.session_state[k] = v


def render_ms_queue_tab(standalone=False):
    """standalone: msQ app, no Tab 3 pipe and no all steps zip."""
    if not standalone:
        st.header("MS Sample Queue")
    for k, v in RUN_PRESETS["Bulk"].items():
        st.session_state.setdefault(k, v)

    for key in ("msq_results", "msq_zip", "msq_last", "msq_group_assignments",
                "msq_block_labels", "msq_csv_hash", "msq_overrides", "msq_last_key",
                "msq_sample_table"):
        if key not in st.session_state:
            st.session_state[key] = None

    date = datetime.date.today().strftime("%Y%m%d")

    def _sample_list_section(title, first):
        """Upload block. msQ: on top with example CSVs. L2-LMD: after the settings, Tab 3 list used."""
        if not first:
            st.divider()
        st.subheader(title)
        pipe_csv  = None if standalone else st.session_state.get("t3_sample_list")
        pipe_stem = None if standalone else st.session_state.get("t3_stem")
        if standalone:
            st.markdown(
                "- **Required**: a sample name column (`ROI`, `Sample`, `sample_name` or `Sample Name`) and a well "
                "column (`Well_ID` or `well`), wells A1 to H12, or A1 to P24 on 384 plates, `A01` style is fine.\n"
                "- **Optional**: `Group`, `Load`, `Inj Vol`, `Plate`, `Dropout {Y/N}`. Column names are not case sensitive.\n"
                "- **CellenONE**: an export with `Location` and `Event Id` works too, the name becomes `Location_EventId`.\n"
            )
            _tpl = find_templates()
            _tc  = st.columns(min(len(_tpl), 4))
            for i, (fn, data) in enumerate(_tpl.items()):
                _tc[i % len(_tc)].download_button(f"Example CSV: {fn}", data, file_name=fn, mime="text/csv",
                                                  key=f"msq_tpl_{fn}")
        uploaded = st.file_uploader("Upload sample list CSV", type=["csv"], key="msq_upload")

        if uploaded:
            csv_b = uploaded.getvalue()
            stm   = re.sub(r'^\d{8}', date, uploaded.name.replace(".csv", ""))
            if st.session_state.msq_last != uploaded.name:
                st.session_state.msq_results          = None
                st.session_state.msq_last             = uploaded.name
                st.session_state.msq_group_assignments = None
            return csv_b, stm
        if pipe_csv is not None:
            st.info("Using sample list piped from Tab 3.")
            return pipe_csv, pipe_stem or "sample_list"
        return None, None

    # Section numbers depend on the app: msQ has the sample list first
    _sec = iter(range(1, 30))
    if standalone:
        csv_bytes, stem = _sample_list_section(f"{next(_sec)}. Sample list", first=True)

    # Instrument
    if standalone:
        st.divider()
    st.subheader(f"{next(_sec)}. Instrument")
    i1, i2 = st.columns(2)
    instrument = i1.radio("Instrument", [INST_ULTRA, INST_ZOOM], horizontal=True, key="msq_instrument")
    lc_system  = i2.radio("LC system", [LC_EVOSEP, LC_NEO], horizontal=True, key="msq_lc_system") \
                 if instrument == INST_ZOOM else LC_EVOSEP
    inst = "ultra" if instrument == INST_ULTRA else ("zoom_neo" if lc_system == LC_NEO else "zoom_evo")
    is_zoom = inst != "ultra"
    i3, i4 = st.columns(2)
    i3.radio("Run mode", list(RUN_PRESETS), horizontal=True, key="msq_run_mode", on_change=_apply_run_preset,
             help="Sets defaults below. SCP: whole plate randomized, blank every 12, K562 and Supermix "
                  "every 24, 384 plates on the Ultra. Everything stays editable.")
    if is_zoom:
        plate_format = "96"
        i4.caption("Plate format: 96 well (Evosep and Vanquish).")
    else:
        plate_format = i4.radio("Sample plate format", list(PLATE_FORMATS), horizontal=True, key="msq_plate_fmt",
                                help="Controls plate is always 96 well.")
    p_rows, p_cols = PLATE_FORMATS[plate_format]
    st.caption({
        "ultra":    "Bruker HyStar queue (XLSX). Positions SlotN:1 to 96 (384: 1 to 384), controls one well per injection.",
        "zoom_evo": "Evosep Slot1 to Slot6, positions SlotN:1 to 96 (A1=1, B1=13), "
                    "controls one tip per injection.",
        "zoom_neo": "Xcalibur sequence (CSV). Vanquish trays R, G, B, Y, positions A1 to H12, "
                    "controls drawn from fixed vials.",
    }[inst])

    # Naming and loads
    st.divider()
    st.subheader(f"{next(_sec)}. Naming and loads")
    if is_zoom:
        st.markdown(
            "- **Pattern**: `date_initials_MS_LC_FAIMS_Acquisition_Gradient_Tags_load_sample`. Empty fields are skipped.\n"
            "- **MS**: instrument, e.g. `Zoom`. **Acquisition**: MS method, e.g. `defaultDIA`, `hybridDIA`.\n"
            "- **FAIMS**: tick to add `FAIMS`. CV `-50` gives `FAIMS_neg50`, a sweep `-40,-60` gives "
            "`FAIMS_neg40_neg60`, empty CV gives `FAIMS`.\n"
            "- **LC**: `VanquishNeo` or `Evosep`. **Gradient**: e.g. `60SPD`, `40SPD`.\n"
            "- **Tags**: free text, e.g. `peptide`.\n"
            "- **Load**: free text, e.g. `500pg`, `2ul`. Empty is left out. Per sample values go in the sample table.\n"
            "- **Inj Vol**: µl per injection type. Per sample values go in the sample table.\n"
            "- **Sample**: name column of the CSV (`ROI`, `Sample`, `sample_name` or `Sample Name`). Repeats get `_r2`, `_r3`."
        )
    else:
        st.markdown(
            "- **Pattern**: `date_initials_LC_MS_Tags_load_sample`.\n"
            "- **LC** and **MS**: method short names, paths are filled in section 3.\n"
            "- **Tags**: free text, e.g. `peptide`.\n"
            "- **Load**: free text, e.g. `20ul`. Per sample values go in the sample table.\n"
            "- **Sample**: name column of the CSV (`ROI`, `Sample`, `sample_name` or `Sample Name`). Repeats get `_r2`, `_r3`."
        )

    t1, t2 = st.columns(2)
    initials = t1.text_input("Initials", value="FK", key="msq_initials", help="Your initials, second field of the name.")
    tags     = t2.text_input("Tags (optional)", value="", key="msq_tags",
                             help="Free text added after the method fields, e.g. peptide. Spaces become _.")
    tags     = re.sub(r"\s+", "_", tags.strip())

    zoom_name = {}
    col_ms, col_lc = st.columns(2)
    with col_ms:
        st.markdown("**MS**")
        if is_zoom:
            zoom_name["ms"] = st.text_input("MS", value=ZOOM_MS_SHORT, key="msq_z_ms",
                                            help="Instrument short name, e.g. Zoom.")
            acq_opts = list(ZOOM_MS_METHODS.keys()) + ["Custom"]
            acq_sel  = st.selectbox("Acquisition", acq_opts, key="msq_z_acq_sel",
                                    help="MS method short name, also picks the default instrument method path.")
            zoom_name["acq"] = st.text_input("Custom acquisition name", value="", key="msq_z_acq_custom") \
                               if acq_sel == "Custom" else acq_sel
            f1, f2 = st.columns([1, 2])
            zoom_name["faims"] = f1.checkbox("FAIMS", value=False, key="msq_z_faims",
                                             help="Adds FAIMS to the name and uses the FAIMS method path.")
            zoom_name["cv"] = f2.text_input("FAIMS CV (optional)", value="", key="msq_z_cv",
                                            help="e.g. -50, or a sweep -40,-60. Written as neg50.") \
                              if zoom_name["faims"] else ""
        else:
            ms_sel   = st.selectbox("MS method", MS_OPTIONS, key="msq_ms_sel_ultra",
                                    help="MS method short name, also picks the default method paths.")
            ms_short = st.text_input("Custom MS short name", value="", key="msq_ms_custom_ultra") \
                       if ms_sel == "Custom" else ms_sel
    with col_lc:
        st.markdown("**LC**")
        if is_zoom:
            zoom_name["lc"] = st.text_input("LC", value=ZOOM_LC_NAME[lc_system], key=f"msq_z_lc_{inst}",
                                            help="LC system short name, set by the LC system in section 1.")
            zoom_name["gradient"] = st.text_input("Gradient", value=ZOOM_GRADIENT[lc_system], key=f"msq_z_grad_{inst}",
                                                  help="Gradient or throughput, e.g. 60SPD or 40SPD.")
        else:
            lc_sel   = st.selectbox("LC method", LC_OPTIONS, key="msq_lc_sel_ultra",
                                    help="LC method short name, also picks the default separation method path.")
            lc_short = st.text_input("Custom LC short name", value="", key="msq_lc_custom_ultra") \
                       if lc_sel == "Custom" else lc_sel
    if is_zoom:
        lc_short, ms_short = zoom_name["lc"], zoom_name["acq"]

    # Standards: one row per injection type
    st.markdown("**Standards and loads**")
    _w = [1, 2, 2]
    inj_vol = {}
    h = st.columns(_w)
    h[1].caption("Load (in the name)")
    if is_zoom:
        h[2].caption("Inj Vol (µl)")

    l1 = st.columns(_w)
    l1[0].markdown("Sample")
    sample_load = l1[1].text_input("Sample load", value="" if is_zoom else "20ul",
                                   key=f"msq_sload_{'zoom' if is_zoom else 'ultra'}", label_visibility="collapsed",
                                   placeholder="optional, e.g. 500pg").strip()

    l2 = st.columns(_w)
    use_k562  = l2[0].checkbox("K562", value=True, key="msq_k562")
    k562_load = l2[1].text_input("K562 load", value="1ng", key="msq_k562l", label_visibility="collapsed",
                                 disabled=not use_k562)
    k562_load = k562_load if use_k562 else ""

    l3 = st.columns(_w)
    use_supermix  = l3[0].checkbox("Supermix", value=True, key="msq_smix")
    supermix_load = l3[1].text_input("Supermix load", value="20ng", key="msq_smixl", label_visibility="collapsed",
                                     disabled=not use_supermix)
    supermix_load = supermix_load if use_supermix else ""

    l4 = st.columns(_w)
    l4[0].markdown("Blank")

    if is_zoom:
        inj_vol["Sample"]   = l1[2].number_input("Inj Vol sample",   min_value=0.0, value=1.0, step=0.1, key="msq_iv_s", label_visibility="collapsed")
        inj_vol["K562"]     = l2[2].number_input("Inj Vol K562",     min_value=0.0, value=1.0, step=0.1, key="msq_iv_k", label_visibility="collapsed", disabled=not use_k562)
        inj_vol["Supermix"] = l3[2].number_input("Inj Vol Supermix", min_value=0.0, value=1.0, step=0.1, key="msq_iv_m", label_visibility="collapsed", disabled=not use_supermix)
        inj_vol["Blank"]    = l4[2].number_input("Inj Vol blank",    min_value=0.0, value=1.0, step=0.1, key="msq_iv_b", label_visibility="collapsed")

    # Live example name
    _pre = name_prefix(dict(date=date, initials=initials, lc_short=lc_short, ms_short=ms_short,
                            tags=tags, inst=inst, zoom_name=zoom_name))
    st.markdown("**Example names**")
    st.code("\n".join(["_".join(x for x in (_pre, sample_load, "H20_core1") if x)]
                      + ([f"{_pre}_{k562_load}_K562_1"] if use_k562 else [])
                      + [f"{_pre}_Blank_1"]), language=None)

    # Methods and data paths
    st.divider()
    st.subheader(f"{next(_sec)}. Methods and data paths")
    with st.expander("Instrument method paths", expanded=False):
        if is_zoom:
            sep_method = inj_method = proc_method = ""
            _fkey      = "FAIMS" if zoom_name["faims"] else "noFAIMS"
            ms_method  = st.text_area("Instrument Method", value=ZOOM_MS_METHODS.get(ms_short, {}).get(_fkey, ""),
                                      height=60, key=f"msq_ms_{inst}_{ms_short}_{_fkey}")
            st.caption("Dummy presets from methods.csv until the real methods are set.")
        else:
            default_sep              = LC_METHODS.get(lc_short, "")
            default_ms, default_proc = MS_METHODS.get(ms_short, ("", ""))
            sep_method  = st.text_area("Separation Method", value=default_sep,   height=60, key="msq_sep")
            inj_method  = st.text_input("Injection Method",  value="Standard",               key="msq_inj")
            ms_method   = st.text_area("MS Method",          value=default_ms,   height=60, key="msq_ms")
            proc_method = st.text_area("Processing Method",  value=default_proc, height=60, key="msq_proc")

    _year      = date[:4]
    _month_str = datetime.datetime.strptime(date, "%Y%m%d").strftime("%m %B")
    _root      = rf"D:\Data\{_year}\{_month_str}"
    with st.expander("Data paths", expanded=True):
        if is_zoom:
            cp1, cp2, cp3 = st.columns(3)
            sample_path = cp1.text_input("Sample path",   value=rf"{_root}\Sample\{initials}", key="msq_spath_zoom")
            std_path    = cp2.text_input("Standard path", value=rf"{_root}\Standard",          key="msq_stdpath_zoom")
            blank_path  = cp3.text_input("Blank path",    value=rf"{_root}\Blank",             key="msq_bpath_zoom")
        else:
            cp1, cp2    = st.columns(2)
            sample_path = cp1.text_input("Sample / K562 / Supermix path",
                                         value=rf"{_root}\Sample\{initials}", key="msq_spath")
            blank_path  = cp2.text_input("Blank path", value=rf"{_root}\Blank", key="msq_bpath")
            std_path    = sample_path

    # Sample list in L2-LMD, after the settings
    if not standalone:
        csv_bytes, stem = _sample_list_section(f"{next(_sec)}. Sample list", first=False)
        if csv_bytes is None:
            return
        st.markdown("**Sample table**")
    else:
        st.divider()
        st.subheader(f"{next(_sec)}. Sample table")
        if csv_bytes is None:
            st.info("Upload a sample list in step 1 to continue.")
            return

    # Reset groups if CSV content changed
    csv_hash = hash(csv_bytes)
    if st.session_state.msq_csv_hash != csv_hash:
        st.session_state.msq_group_assignments = None
        st.session_state.msq_results           = None
        st.session_state.msq_block_labels      = None
        st.session_state.msq_sample_table      = None
        st.session_state.msq_csv_hash          = csv_hash

    st.caption(f"Output stem: `{stem}`")

    # Parse samples, sample name column can be ROI, Sample, sample_name or Sample Name
    all_rows, fields, name_col = read_sample_rows(csv_bytes)
    if not name_col:
        st.error("No sample name column found. Expected one of: " + ", ".join(SAMPLE_COLS) + ".")
        return
    if "Well_ID" not in fields:
        st.error("No `Well_ID` column found in the sample list.")
        return
    if name_col != "ROI":
        st.caption(f"Sample names read from column `{name_col}`.")
    samples  = [r for r in all_rows if r.get("Dropout {Y/N}", "").strip().upper() != "Y"]

    if not samples:
        st.error("No active samples in CSV (all marked as dropouts).")
        return

    bad_wells = []
    for i, r in enumerate(all_rows):
        w, why = check_well(r.get("Well_ID", ""), p_rows, p_cols)
        if w is None:
            bad_wells.append({"CSV row": i + 2, "Sample": r.get("ROI", ""), "Well": r.get("Well_ID", ""), "Problem": why})
    if bad_wells:
        st.error(f"{len(bad_wells)} row(s) have a well that is not valid for a {plate_format} well plate.")
        st.dataframe(pd.DataFrame(bad_wells), hide_index=True, use_container_width=True)
        return
    _wells = [normalize_well(r["Well_ID"]) + "|" + r.get("Plate", "") for r in samples]
    _dupes = sorted({w.split("|")[0] for w in _wells if _wells.count(w) > 1})
    if _dupes:
        st.warning("Wells used more than once: " + ", ".join(_dupes))

    # Detect plates from CSV
    has_plate_col = "Plate" in fields
    plate_order   = sorted({r.get("Plate", "Plate1").strip()
                             for r in all_rows}) if has_plate_col else ["Plate1"]

    # Per-sample table, applied with the Update queue button
    st.caption("Empty Load or Inj Vol uses the global value. To repeat a sample, add a row and pick "
               "it again: it runs from the same well directly after the original, named `_r2`, `_r3`. "
               "Edits apply when you press Update queue.")
    names    = [r["ROI"].strip() for r in samples]
    well_of  = {r["ROI"].strip(): r.get("Well_ID", "").strip() for r in samples}
    plate_of = {r["ROI"].strip(): r.get("Plate", "Plate1").strip() for r in samples}
    if st.session_state.get("msq_sample_table") is None:
        def _csv_vol(r):
            v = str(r.get("Inj Vol", "") or "").strip()
            try:
                return float(v)
            except ValueError:
                return None
        st.session_state.msq_sample_table = [
            {"Sample": n, "Load": str(r.get("Load", "") or "").strip(), "Inj Vol": _csv_vol(r)}
            for n, r in zip(names, samples)
        ]
    table_df = pd.DataFrame(
        [{"Sample": t["Sample"], "Well": well_of.get(t["Sample"], ""), "Plate": plate_of.get(t["Sample"], ""),
          "Load": t["Load"], "Inj Vol": t["Inj Vol"]} for t in st.session_state.msq_sample_table],
        columns=["Sample", "Well", "Plate", "Load", "Inj Vol"],
    )
    table_df["Inj Vol"] = table_df["Inj Vol"].astype("float64")
    with st.form("msq_sample_form", border=False):
        edited_tbl = st.data_editor(
            table_df,
            column_config={
                "Sample":  st.column_config.SelectboxColumn("Sample", options=names, required=True),
                "Well":    st.column_config.TextColumn("Well", disabled=True),
                "Plate":   st.column_config.TextColumn("Plate", disabled=True),
                "Load":    st.column_config.TextColumn("Load", help="e.g. 500pg, 2ul. Empty uses the global load."),
                "Inj Vol": st.column_config.NumberColumn("Inj Vol (µl)", min_value=0.0, step=0.1, format="%.2f",
                                                         help="Empty uses the global injection volume."),
            },
            column_order=["Sample", "Well"] + (["Plate"] if has_plate_col else []) + ["Load", "Inj Vol"],
            num_rows="dynamic", hide_index=True, use_container_width=True, key="msq_sample_editor",
        )
        submitted = st.form_submit_button("Update queue", type="primary")
    if submitted:
        new_tbl = []
        for _, r in edited_tbl.iterrows():
            if pd.isna(r["Sample"]) or not str(r["Sample"]).strip():
                continue
            new_tbl.append({"Sample": str(r["Sample"]).strip(),
                            "Load": "" if pd.isna(r["Load"]) else str(r["Load"]).strip(),
                            "Inj Vol": None if pd.isna(r["Inj Vol"]) else float(r["Inj Vol"])})
        st.session_state.msq_sample_table = new_tbl
        st.session_state.msq_results      = None
        st.rerun()
    sample_table = st.session_state.msq_sample_table
    _listed  = {t["Sample"] for t in sample_table}
    _n_rep   = len(sample_table) - len(_listed)
    _missing = [n for n in names if n not in _listed]
    st.caption(f"{len(names)} samples, {_n_rep} repeat injection(s).")
    if _missing:
        st.warning(f"{len(_missing)} sample(s) removed from the table still run with global values. "
                   f"Mark dropouts in the CSV to exclude them: " + ", ".join(_missing[:10])
                   + (" ..." if len(_missing) > 10 else ""))

    # Randomization
    st.divider()
    st.subheader(f"{next(_sec)}. Randomization")
    st.caption("Samples are shuffled first, blanks and standards are then inserted by the rules below. "
               "Each sample keeps its own well. Run number, well and seed go to run_order.csv"
               + (" and the Xcalibur Comment column." if is_zoom else "."))
    _ra1, _ra2 = st.columns([2, 1])
    randomize_mode = _ra1.radio("Randomize run order", ["Off", "Within group", "Whole plate"],
                                horizontal=True, key="msq_rand_mode",
                                help="Whole plate ignores groups and blocks: one block per plate.")
    run_seed = _ra2.number_input("Random seed", min_value=0, max_value=9999, value=42, step=1,
                                 key="msq_run_seed", disabled=randomize_mode == "Off")

    # Plates and control positions
    st.divider()
    st.subheader(f"{next(_sec)}. Plates and control positions")
    fixed_positions = {}
    saved_slot_map  = st.session_state.get(f"msq_plate_slot_map_{inst}") or {}

    if inst == "zoom_neo":
        st.caption("One tray per sample plate. Controls are drawn from fixed vials, edit if needed.")
        free_trays = [t for t in NEO_TRAYS if t != "R"] + ["R"]
        slot_df = pd.DataFrame([
            {"Plate": pl, "Tray": saved_slot_map.get(pl, {}).get("sample_slot", free_trays[i % 4])}
            for i, pl in enumerate(plate_order)
        ])
        edited_slots = st.data_editor(
            slot_df,
            column_config={
                "Plate": st.column_config.TextColumn("Plate", disabled=True),
                "Tray":  st.column_config.SelectboxColumn("Tray", options=NEO_TRAYS),
            },
            hide_index=True, use_container_width=True, key=f"msq_slot_editor_{inst}",
        )
        plate_slot_map = {row["Plate"]: {"sample_slot": row["Tray"], "ctrl_start_slot": 0}
                          for _, row in edited_slots.iterrows()}

        f1, f2, f3 = st.columns(3)
        fixed_positions["Blank"]    = f1.text_input("Blank vial",    value=NEO_FIXED_DEFAULT["Blank"],    key="msq_fix_b").strip().upper()
        fixed_positions["K562"]     = f2.text_input("K562 vial",     value=NEO_FIXED_DEFAULT["K562"],     key="msq_fix_k",
                                                    disabled=not use_k562).strip().upper()
        fixed_positions["Supermix"] = f3.text_input("Supermix vial", value=NEO_FIXED_DEFAULT["Supermix"], key="msq_fix_m",
                                                    disabled=not use_supermix).strip().upper()
        bad = [f"{k} `{v}`" for k, v in fixed_positions.items() if not NEO_POS_RE.match(v)]
        if bad:
            st.error("Invalid vial position: " + ", ".join(bad) + ". Use tray:well, e.g. R:F9.")
            return
        if len(plate_order) > len(NEO_TRAYS):
            st.error(f"{len(plate_order)} plates but only {len(NEO_TRAYS)} trays.")
            return

        # Sample wells that collide with a fixed control vial
        tray_of = {pl: v["sample_slot"] for pl, v in plate_slot_map.items()}
        fixed_set = {v: k for k, v in fixed_positions.items()
                     if k == "Blank" or (k == "K562" and use_k562) or (k == "Supermix" and use_supermix)}
        clashes = []
        for r in samples:
            if not r.get("Well_ID", "").strip():
                continue
            pl  = r.get("Plate", "Plate1").strip() if has_plate_col else "Plate1"
            pos = f"{tray_of.get(pl, '')}:{normalize_well(r['Well_ID'])}"
            if pos in fixed_set:
                clashes.append(f"{r['ROI'].strip()} at {pos} ({fixed_set[pos]})")
        if clashes:
            st.error("Sample wells collide with control vials: " + "; ".join(clashes))
            return
        dup_trays = {t for t in tray_of.values() if list(tray_of.values()).count(t) > 1}
        if dup_trays:
            st.error("Two plates on the same tray: " + ", ".join(sorted(dup_trays)))
            return
    else:
        slot_opts = EVOSEP_SLOTS if inst == "zoom_evo" else [f"Slot{i}" for i in range(1, 13)]
        st.caption("Auto-assigned consecutively (Plate1: Slot1+Slot2, Plate2: Slot3+Slot4, ...). Edit to override.")
        slot_data = []
        for i, pl in enumerate(plate_order):
            saved = saved_slot_map.get(pl, {})
            slot_data.append({
                "Plate":               pl,
                "Sample slot":         saved.get("sample_slot",         f"Slot{i * 2 + 1}"),
                "Controls start slot": saved.get("ctrl_start_slot_str", f"Slot{i * 2 + 2}"),
            })
        edited_slots = st.data_editor(
            pd.DataFrame(slot_data),
            column_config={
                "Plate":               st.column_config.TextColumn("Plate", disabled=True),
                "Sample slot":         st.column_config.SelectboxColumn("Sample slot",         options=slot_opts),
                "Controls start slot": st.column_config.SelectboxColumn("Controls start slot", options=slot_opts),
            },
            hide_index=True, use_container_width=True, key=f"msq_slot_editor_{inst}",
        )
        plate_slot_map = {}
        for _, row in edited_slots.iterrows():
            plate_slot_map[row["Plate"]] = {
                "sample_slot":         row["Sample slot"],
                "ctrl_start_slot":     int(row["Controls start slot"].replace("Slot", "")),
                "ctrl_start_slot_str": row["Controls start slot"],
            }
    st.session_state[f"msq_plate_slot_map_{inst}"] = plate_slot_map


    # Read Group column from CSV if present (piped from Tab 2/3)
    has_group_col = "Group" in fields
    csv_group_map = {}
    if has_group_col:
        for r in all_rows:
            roi = r.get("ROI", "").strip()
            grp = r.get("Group", "").strip()
            if roi and grp:
                csv_group_map[roi] = grp

    confirmed_assignments = st.session_state.msq_group_assignments or {}

    # Auto-number groups sequentially (Group 1, 2, ...) in order of first appearance
    _seen_sg, _roi_to_group_label = [], {}
    for r in samples:
        roi = r["ROI"].strip()
        sg  = suggest_group(roi)
        if sg not in _seen_sg:
            _seen_sg.append(sg)
        _roi_to_group_label[roi] = f"Group {_seen_sg.index(sg) + 1}"

    init_data = [
        {
            "ROI":   r["ROI"].strip(),
            "Well":  r.get("Well_ID", "").strip(),
            "Group": confirmed_assignments.get(r["ROI"].strip()) or _roi_to_group_label[r["ROI"].strip()],
            "Name":  suggest_group(r["ROI"].strip()),
        }
        for r in samples
    ]
    group_df = pd.DataFrame(init_data)

    # Block editor, one row per block key
    st.divider()
    st.subheader(f"{next(_sec)}. Block assignment")
    st.caption(
        "K562 + Supermix + Blank at each block start, counters reset there. "
        "Inside a block the samples of all its groups are batched together, a Blank runs "
        "at least every 6 samples and a Supermix at least every 12. "
        "Run ends on K562, Supermix, Blank. "
        "Blocks come from the first token of the sample name. Edit a label to rename a "
        "block, or give two rows the same label to merge them."
    )

    # Group order and representative sample name
    seen_order, group_label_to_name = [], {}
    for row in init_data:
        if row["Group"] not in seen_order:
            seen_order.append(row["Group"])
            group_label_to_name[row["Group"]] = row["Name"]

    # Block key per group, numbered by first appearance
    group_prefix, seen_prefixes = {}, []
    for g in seen_order:
        pre = suggest_block(group_label_to_name[g])
        group_prefix[g] = pre
        if pre not in seen_prefixes:
            seen_prefixes.append(pre)
    auto_block = {pre: f"Block {i + 1}" for i, pre in enumerate(seen_prefixes)}

    n_groups, n_rois = {}, {}
    for g in seen_order:
        n_groups[group_prefix[g]] = n_groups.get(group_prefix[g], 0) + 1
    for row in init_data:
        pre = group_prefix[row["Group"]]
        n_rois[pre] = n_rois.get(pre, 0) + 1

    saved_labels = st.session_state.msq_block_labels or {}
    block_data   = [
        {
            "Block":  saved_labels.get(pre, auto_block[pre]),
            "Sample": pre,
            "Groups": n_groups[pre],
            "ROIs":   n_rois[pre],
        }
        for pre in seen_prefixes
    ]

    edited_blocks = st.data_editor(
        pd.DataFrame(block_data),
        column_config={
            "Block":  st.column_config.TextColumn("Block", help="Rename, or reuse a label to merge blocks"),
            "Sample": st.column_config.TextColumn("Sample prefix", disabled=True),
            "Groups": st.column_config.NumberColumn("Groups", disabled=True),
            "ROIs":   st.column_config.NumberColumn("ROIs",   disabled=True),
        },
        column_order=["Block", "Sample", "Groups", "ROIs"],
        hide_index=True, use_container_width=True, key="msq_block_editor",
    )

    new_labels = dict(zip(seen_prefixes, edited_blocks["Block"]))
    if new_labels != saved_labels:
        st.session_state.msq_block_labels = new_labels
        st.session_state.msq_results      = None
    block_assignments = {g: new_labels.get(group_prefix[g], auto_block[group_prefix[g]])
                         for g in seen_order}

    # Block summary
    block_summary = {}
    for g in seen_order:
        block_summary.setdefault(block_assignments[g], []).append(g)
    st.info("Blocks: " + "  |  ".join(
        f"**{b}** ({len(gs)} groups)" for b, gs in block_summary.items()
    ))

    # Control intervals
    st.divider()
    st.subheader(f"{next(_sec)}. Control intervals")
    _ro1, _ro2, _ro3, _ro4 = st.columns(4)
    block_size = _ro1.number_input("Blank interval (samples)", min_value=1, max_value=96, step=1,
                                   key="msq_block_size",
                                   help="Max samples between Blank injections.")
    k562_with_smix = _ro3.checkbox("K562 with every Supermix", key="msq_k562_with_smix", disabled=not use_k562,
                                   help="Inside a block K562 also runs at the standard interval, as in SCP. "
                                        "Off: K562 only at block starts.")
    _std_on = use_supermix or (k562_with_smix and use_k562)
    smix_interval = _ro2.number_input("Standard interval (samples)", min_value=1, max_value=192, step=1,
                                      key="msq_smix_interval", disabled=not _std_on,
                                      help="Max samples between Supermix injections, plus K562 if ticked.")
    run_len_min = _ro4.number_input("Run length (min)", min_value=1, max_value=240,
                                    value=22 if is_zoom else 40, step=1, key=f"msq_run_len_{'zoom' if is_zoom else 'ultra'}",
                                    help="Per injection, used for the runtime estimate.")

    if _std_on and int(smix_interval) < int(block_size):
        st.warning("Standard interval is below the Blank interval. Standards can only be "
                   "placed at a batch boundary, so the interval cannot be met.")

    group_assignments = {row["ROI"]: row["Group"] for row in init_data}
    overrides         = st.session_state.msq_overrides or []

    params = dict(
        date=date, initials=initials, lc_short=lc_short, ms_short=ms_short, tags=tags,
        sample_load=sample_load, k562_load=k562_load, supermix_load=supermix_load,
        use_k562=use_k562, use_supermix=use_supermix,
        sep_method=sep_method, inj_method=inj_method,
        ms_method=ms_method, proc_method=proc_method,
        sample_path=sample_path, blank_path=blank_path, std_path=std_path,
        stem=stem, inst=inst, inj_vol=inj_vol, fixed_positions=fixed_positions,
        plate_format=plate_format, k562_with_smix=bool(k562_with_smix),
        zoom_name=zoom_name, sample_table=sample_table,
        plate_slot_map=plate_slot_map,
        randomize_mode=randomize_mode,
        block_size=block_size,
        smix_interval=smix_interval,
        run_seed=run_seed,
        overrides=overrides,
    )

    # Rebuild only when something that changes the queue changed
    build_key = (
        csv_hash, inst, tags, std_path, tuple(sorted(zoom_name.items())), plate_format, bool(k562_with_smix),
        tuple((t["Sample"], t["Load"], t["Inj Vol"]) for t in sample_table),
        tuple(sorted(inj_vol.items())), tuple(sorted(fixed_positions.items())),
        tuple(sorted(group_assignments.items())),
        tuple(sorted(block_assignments.items())),
        tuple(sorted((pl, v["sample_slot"], v["ctrl_start_slot"])
                     for pl, v in plate_slot_map.items())),
        randomize_mode, int(block_size), int(smix_interval), int(run_seed),
        bool(use_k562), bool(use_supermix),
        sample_load, k562_load, supermix_load, initials, lc_short, ms_short,
        sep_method, inj_method, ms_method, proc_method, sample_path, blank_path, stem,
        tuple(sorted((o["action"], o["type"], o["after"]) for o in overrides)),
    )

    # Run sequence and adjustments
    st.divider()
    st.subheader(f"{next(_sec)}. Run sequence and adjustments")
    if st.session_state.msq_last_key != build_key or not st.session_state.msq_results:
        with st.spinner("Building queue..."):
            res = build_queue_core(csv_bytes, group_assignments, params,
                                   block_assignments=block_assignments)
            st.session_state.msq_results  = res
            st.session_state.msq_zip      = build_zip(res)
            st.session_state.msq_last_key = build_key

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
        + ("Fixed vials: " + ", ".join(f"{k} {v}" for k, v in fixed_positions.items())
           if inst == "zoom_neo" else f"Control slots used: {n_cslots}")
    )
    if inst == "zoom_evo" and res.get("max_ctrl_slot", 0) > len(EVOSEP_SLOTS):
        st.error(f"Controls spill into Slot{res['max_ctrl_slot']}, the Evosep has "
                 f"{len(EVOSEP_SLOTS)} slots.")
    _total_min = res["n_queue"] * int(run_len_min)
    st.caption(f"Estimated runtime: {res['n_queue']} injections x {int(run_len_min)} min = "
               f"{_total_min / 60:.1f} h ({_total_min / 1440:.1f} days)")

    if res.get("max_blank_gap", 0) > int(block_size):
        st.warning(f"Blank gap reaches {res['max_blank_gap']} samples, above the "
                   f"{int(block_size)} sample setting.")
    if use_supermix and res.get("max_smix_gap", 0) > int(smix_interval):
        st.warning(f"Supermix gap reaches {res['max_smix_gap']} samples, above the "
                   f"{int(smix_interval)} sample setting.")
    _n_removed = sum(1 for o in overrides if o.get("action") == "remove")
    if _n_removed:
        st.warning(f"{_n_removed} manual removal(s) active. Removals bypass the automatic "
                   f"rules, including the control set at a block start.")

    with st.expander("Run sequence", expanded=True):
        seq_rows = [
            {"#": n + 1, ("Position" if is_zoom else "Vial"): r.get("Vial", ""),
             ("File Name" if is_zoom else "Sample ID"): r.get("Sample ID", ""),
             "Inj Vol": r["_inj_vol"] if r.get("_inj_vol") is not None
                        else (inj_vol.get(r.get("_type", "Sample"), 1) if is_zoom else r.get("Volume [µl]", "")),
             "Comment": r.get("_comment", "")}
            for n, r in enumerate(res.get("queue_rows", []))
        ]
        st.dataframe(pd.DataFrame(seq_rows), hide_index=True,
                     use_container_width=True, height=360)

    # Adjust
    st.markdown("**Adjust groups and controls**")
    st.caption(
        "Groups are auto-suggested by stripping trailing numbers from sample names. "
        "Edit the grouping or force a control injection, the queue rebuilds on change."
    )

    gsum = {}
    for row in init_data:
        gsum[row["Group"]] = gsum.get(row["Group"], 0) + 1
    st.info("Groups: " + "  |  ".join(f"**{g}** ({n})" for g, n in gsum.items()))

    _sp1, _sp2, _sp3 = st.columns([2, 1, 1])
    split_grp = _sp1.selectbox("Divide group", list(gsum.keys()), key="msq_split_grp")
    split_n   = _sp2.number_input("Into N parts", min_value=2, max_value=20, value=2, step=1,
                                  key="msq_split_n")
    if _sp3.button("Divide", key="msq_split_btn", use_container_width=True):
        rois_in_grp = [row["ROI"] for row in init_data if row["Group"] == split_grp]
        chunk       = max(1, math.ceil(len(rois_in_grp) / int(split_n)))
        updated     = {row["ROI"]: row["Group"] for row in init_data}
        for n, roi in enumerate(rois_in_grp):
            updated[roi] = f"{split_grp}{chr(ord('a') + n // chunk)}"
        st.session_state.msq_group_assignments = updated
        st.rerun()

    edited_df = st.data_editor(
        group_df,
        column_config={
            "Group": st.column_config.TextColumn("Group", help="Edit to reassign to a different group"),
            "Name":  st.column_config.TextColumn("Name",  disabled=True),
            "ROI":   st.column_config.TextColumn("ROI",   disabled=True),
            "Well":  st.column_config.TextColumn("Well",  disabled=True),
        },
        column_order=["Group", "Name", "ROI", "Well"],
        hide_index=True, use_container_width=True, key="msq_group_editor",
    )
    new_assignments = dict(zip(edited_df["ROI"], edited_df["Group"]))
    if new_assignments != group_assignments:
        st.session_state.msq_group_assignments = new_assignments
        st.rerun()

    st.markdown("**Control overrides**")
    st.caption("Add or remove one injection after a chosen sample, on top of the automatic rules.")
    ov_cols = ["Action", "Control", "After ROI"]
    ov_df   = pd.DataFrame(
        [{"Action": o["action"], "Control": o["type"], "After ROI": o["after"]}
         for o in overrides],
        columns=ov_cols,
    )
    edited_ov = st.data_editor(
        ov_df,
        column_config={
            "Action":    st.column_config.SelectboxColumn("Action",  options=["add", "remove"]),
            "Control":   st.column_config.SelectboxColumn("Control", options=["K562", "Supermix", "Blank"]),
            "After ROI": st.column_config.SelectboxColumn("After ROI",
                                                          options=[r["ROI"] for r in init_data]),
        },
        num_rows="dynamic", hide_index=True, use_container_width=True,
        key="msq_override_editor",
    )
    new_ov = [
        {"action": r["Action"], "type": r["Control"], "after": r["After ROI"]}
        for _, r in edited_ov.iterrows()
        if all(pd.notna(r[c]) and str(r[c]).strip() for c in ov_cols)
    ]
    if new_ov != overrides:
        st.session_state.msq_overrides = new_ov
        st.rerun()

    # Plate maps
    st.divider()
    st.subheader(f"{next(_sec)}. Plate maps")
    all_plate_items = []
    for key, data in res["sample_slot_outputs"].items():
        all_plate_items.append((f"{data['sample_slot']} - {data['plate']} Samples", data["png"]))
    for sn, data in res["ctrl_slots"].items():
        all_plate_items.append((f"Slot{sn} - Controls", data["png"]))
    for i in range(0, len(all_plate_items), 2):
        chunk = all_plate_items[i:i + 2]
        cols  = st.columns(len(chunk))
        for col, (title, png) in zip(cols, chunk):
            with col:
                st.markdown(f"**{title}**")
                st.image(png, use_container_width=True)

    # Downloads
    st.divider()
    st.subheader(f"{next(_sec)}. Downloads")
    stem_out = res["stem"]

    def build_all_zip(res, msq_zip_bytes):
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
            if st.session_state.get("t1_geojson") and st.session_state.get("t1_stem"):
                z.writestr(f"{st.session_state.t1_stem}_reclassified.geojson",
                           st.session_state.t1_geojson)
            if st.session_state.get("t2_zip") and st.session_state.get("t2_stem"):
                # Unpack Tab 2 zip into all_steps zip
                with zipfile.ZipFile(io.BytesIO(st.session_state.t2_zip)) as t2z:
                    for name in t2z.namelist():
                        z.writestr(name, t2z.read(name))
            if st.session_state.get("proc_zip"):
                with zipfile.ZipFile(io.BytesIO(st.session_state.proc_zip)) as p3z:
                    for name in p3z.namelist():
                        z.writestr(name, p3z.read(name))
            if st.session_state.get("t3_sample_list"):
                z.writestr(f"{stem_out}_sample_list.csv", st.session_state.t3_sample_list)
            with zipfile.ZipFile(io.BytesIO(msq_zip_bytes)) as msq_z:
                for name in msq_z.namelist():
                    z.writestr(name, msq_z.read(name))
        return buf.getvalue()

    dl_all, dl_ms = st.columns(2)
    if not standalone:
        dl_all.download_button(
            "Download all steps (zip)", build_all_zip(res, st.session_state.msq_zip),
            file_name=f"{stem_out}_all_steps.zip", mime="application/zip", type="primary"
        )
    dl_ms.download_button(
        "Download MS queue (zip)" if standalone else "Download MS queue only (zip)", st.session_state.msq_zip,
        file_name=f"{stem_out}_ms_queue.zip", mime="application/zip", type="primary" if standalone else "secondary"
    )

    st.divider()

    # Individual file downloads
    dl_items = []
    if res.get("queue_csv"):
        dl_items.append(("Queue CSV (Xcalibur)", res["queue_csv"], f"{stem_out}_queue_xcalibur.csv", "text/csv"))
    if res.get("run_order_csv"):
        dl_items.append(("Run order CSV", res["run_order_csv"], f"{stem_out}_run_order.csv", "text/csv"))
    if res["queue_xlsx"]:
        dl_items.append(("Queue XLSX", res["queue_xlsx"],
                         f"{stem_out}_queue.xlsx",
                         "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"))
    for key, data in res["sample_slot_outputs"].items():
        dl_items.append((f"{key} CSV", data["csv"], f"{stem_out}_{key}_samples.csv", "text/csv"))
        dl_items.append((f"{key} PNG", data["png"], f"{stem_out}_{key}_samples.png", "image/png"))
    for sn, data in res["ctrl_slots"].items():
        dl_items.append((f"Slot{sn} CSV", data["csv"], f"{stem_out}_slot{sn}.csv", "text/csv"))
        dl_items.append((f"Slot{sn} PNG", data["png"], f"{stem_out}_slot{sn}.png", "image/png"))

    cols = st.columns(min(len(dl_items), 5))
    for i, (label, data, fname, mime) in enumerate(dl_items):
        cols[i % 5].download_button(label, data, file_name=fname, mime=mime)
