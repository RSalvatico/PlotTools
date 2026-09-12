"""Shared CMS/mplhep styling for the rabbit fit plots (stack, scan, impacts).

One place for:
  * the CMS label convention already used elsewhere in this repo
    (plotUnstacked.py: "Work in Progress", 110 fb^-1, sqrt(s) = 13.6 TeV);
  * the physics grouping of the 17 raw templates into 9 legend entries
    (ttLF, tt2b, ttbb, ttbj, tt2c, ttcc, ttcj, ttWcb, others) with the fixed
    cmsstyle p10 colors, so every plot (any category, any fit) uses the same
    color for the same process group;
  * display labels for the 8 categories on the unrolled stack plot.
"""

import mplhep as hep
import numpy as np

hep.style.use(hep.style.CMS)

LUMI = "110"
COM = "13.6"


def cms_label(ax, data, loc=2, **kwargs):
    """CMS exp label ('Work in Progress'). data=False prepends 'Simulation'."""
    hep.cms.label("Work in Progress", loc=loc, ax=ax, data=data, lumi=LUMI,
                  com=COM, **kwargs)


# https://cmsstyle.readthedocs.io/en/latest/reference/#cmsstyle.cmsstyle.p10
COLOURS = {
    "ttWcb": "#ffa90e",
    "ttLF": "#e76300",
    "ttcj": "#b9ac70",
    "tt2c": "#94a4a2",
    "ttcc": "#717581",
    "ttbj": "#92dadd",
    "tt2b": "#3f90da",
    "ttbb": "#832db6",
    "others": "#a96b59",
}

# raw process -> plotted group, in stacking order (bottom to top).
# (key into COLOURS, legend label, member raw processes)
PROCESS_GROUPS = [
    ("ttLF", r"$\mathrm{t\bar{t}}$ LF", ["ttLF"]),
    ("tt2b", r"$\mathrm{t\bar{t}}$ 2b", ["tt2b", "tt2b-dps"]),
    ("ttbb", r"$\mathrm{t\bar{t}}$ bb", ["ttbb", "ttbb-dps"]),
    ("ttbj", r"$\mathrm{t\bar{t}}$ bj", ["ttbj", "ttbj-dps"]),
    #("tt2b", r"$\mathrm{t\bar{t}}$ 2b", ["tt2b"]),
    #("ttbb", r"$\mathrm{t\bar{t}}$ bb", ["ttbb"]),
    #("ttbj", r"$\mathrm{t\bar{t}}$ bj", ["ttbj"]),
    ("tt2c", r"$\mathrm{t\bar{t}}$ 2c", ["tt2c"]),
    ("ttcc", r"$\mathrm{t\bar{t}}$ cc", ["ttcc"]),
    ("ttcj", r"$\mathrm{t\bar{t}}$ cj", ["ttcj"]),
    ("ttWcb", r"$\mathrm{t\bar{t}}$ (Vcb)", ["tt-vcb"]),
    # everything else: single top, W+jets, ttZ/ttW, diboson, ttH
    ("others", "Other bkg.", ["singletop", "wjets", "ttZ", "ttW", "diboson",
                              "ttHbb", "ttHcc"]),
]

GROUP_COLORS = {label: COLOURS[key] for key, label, _ in PROCESS_GROUPS}


# category (from configs.model.CATEGORIES) -> display label for the unrolled
# stack plot.
REGION_LABELS = {
    "Vcb_catWcb_SR": "SR",
    "Vcb_catBB_CR": r"$\mathrm{t\bar{t}}$+bb CR",
    "Vcb_catBJ_CR": r"$\mathrm{t\bar{t}}$+bj CR",
    "Vcb_cat2B_CR": r"$\mathrm{t\bar{t}}$+2b CR",
    "Vcb_catCC_CR": r"$\mathrm{t\bar{t}}$+cc CR",
    "Vcb_catCJ_CR": r"$\mathrm{t\bar{t}}$+cj CR",
    "Vcb_cat2C_CR": r"$\mathrm{t\bar{t}}$+2c CR",
    "Vcb_catLF_CR": r"$\mathrm{t\bar{t}}$+LF CR",
}


def region_label(channel):
    return REGION_LABELS[channel]


def group_stack(stack):
    """{raw process: array} -> (labels, colors, arrays) in draw order."""
    known = {p for _, _, procs in PROCESS_GROUPS for p in procs}
    assert set(stack) <= known, f"process(es) not in PROCESS_GROUPS: {set(stack) - known}"
    labels, colors, arrays = [], [], []
    for _, label, procs in PROCESS_GROUPS:
        present = [stack[p] for p in procs if p in stack]
        if not present:
            continue
        labels.append(label)
        colors.append(GROUP_COLORS[label])
        arrays.append(np.sum(present, axis=0))
    return labels, colors, arrays
