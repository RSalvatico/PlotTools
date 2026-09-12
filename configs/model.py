"""Fit model definition for the Vcb analysis.

Single source of truth for everything that is *not* in the merged shapes file:
category -> distribution map, process lists, lnN normalisation uncertainties,
free rate parameters, nuisance groups, and the per-systematic metadata the
tensor pipeline needs (statistical correlation class, symmetrisation policy).

Consumed by analysis/prepareTensor.py. prepareDatacards.py still carries its own
copy; wiring it to this module is tracked in COMBINE_FIXES.md.
"""

CHANNEL = "SL"

# Bin content used by fixNegativeBins.py / the shape production to stand in for
# an empty or negative bin. Combine cannot take a literal 0 for a process with
# shape systematics; rabbit can, so the tensor writer undoes this substitution.
SENTINEL = 1e-6


# --------------------------------------------------------------------------
# categories
# --------------------------------------------------------------------------

CATEGORIES = {
    "Vcb_catWcb_SR": "score_tt_Wcb",
    "Vcb_catBB_CR": "fscore_ttbb",
    "Vcb_catBJ_CR": "fscore_ttbj",
    "Vcb_cat2B_CR": "fscore_tt2b",
    "Vcb_catCC_CR": "fscore_ttcc",
    "Vcb_catCJ_CR": "fscore_ttcj",
    "Vcb_cat2C_CR": "fscore_tt2c",
    "Vcb_catLF_CR": "fscore_ttLF",
}

SIGNAL_REGIONS = [c for c in CATEGORIES if c.endswith("_SR")]
CONTROL_REGIONS = [c for c in CATEGORIES if c.endswith("_CR")]


# --------------------------------------------------------------------------
# processes
# --------------------------------------------------------------------------

SIGNAL = ["tt-vcb"]

BACKGROUNDS = [
    "singletop", "ttbb", "ttbj", "tt2b", "ttbb-dps", "ttbj-dps", "tt2b-dps",
    #"singletop", "ttbb", "ttbj", "tt2b",
    "ttcc", "ttcj", "tt2c", "ttLF", "wjets", "ttZ", "ttW", "diboson",
    "ttHbb", "ttHcc",
]

ALL_PROCESSES = SIGNAL + BACKGROUNDS

TT_COMPONENTS_MAIN_BKG_NODPS = [
    "ttbb", "ttbj", "tt2b", "ttcc", "ttcj", "tt2c", "ttLF",
]
TT_COMPONENTS_BBDPS = ["ttbb-dps", "ttbj-dps", "tt2b-dps"]
TTH_COMPONENTS = ["ttHbb", "ttHcc"]


def disable_dps():
    """Drop the tt+bb DPS split (ttbb-dps, ttbj-dps, tt2b-dps) from the process
    list, in place.

    Call this once, before read_shapes(), when the shapes file's ttbar
    simulation is 5FS-only: 5FS MC has no separate double-parton-shower
    component, so the merged shapes file carries no histograms for the -dps
    processes at all, and validate() would otherwise reject them as declared
    processes with no nominal. norm_ttbb-dps then targets an empty process
    list and lnN_systematics() drops it automatically.
    """
    global BACKGROUNDS, ALL_PROCESSES, TT_COMPONENTS_BBDPS
    BACKGROUNDS = [p for p in BACKGROUNDS if p not in TT_COMPONENTS_BBDPS]
    ALL_PROCESSES = SIGNAL + BACKGROUNDS
    TT_COMPONENTS_BBDPS = []


# --------------------------------------------------------------------------
# normalisation
# --------------------------------------------------------------------------

# Processes whose normalisation floats freely. These become POIs of the
# FreeNorm param model, mirroring the `xsec_<proc> rateParam` lines and the
# `--redefineSignalPOIs xsec_tt*` used in the Combine CR fits.
FREE_NORM_PROCESSES = list(TT_COMPONENTS_MAIN_BKG_NODPS)


def lnN_systematics(year):
    """name -> (kappa, [processes]). Mirrors prepareDatacards.py."""
    if str(year) != "2024":
        raise ValueError(f"lnN systematics are only defined for 2024, got {year}")
    all_systs = {
        # Values copied from Datacards_*/orig/Vcb_SL_2024.txt. They were wrong
        # here for norm_singletop, norm_wjets, norm_diboson and norm_ttZ, which
        # removed real freedom from the fit (singletop is 7.4% of the CR
        # prediction). A kappa may be a scalar or a (kappa_up, kappa_down) pair.
        f"CMS_lumi_13p6TeV_{year}": (1.016, list(ALL_PROCESSES)),
        "norm_singletop": (1.25, ["singletop"]),
        "norm_ttW": (1.068, ["ttW"]),
        "norm_ttZ": ((1.096, 1.085), ["ttZ"]),
        "norm_ttbb-dps": (1.50, list(TT_COMPONENTS_BBDPS)),
        "norm_ttH": (1.20, list(TTH_COMPONENTS)),
        "norm_diboson": (1.30, ["diboson"]),
        "norm_wjets": (1.30, ["wjets"]),
    }
    # A kappa targeting no process (e.g. norm_ttbb-dps after disable_dps())
    # is a nuisance with nothing to do; drop it rather than declare a no-op.
    return {name: v for name, v in all_systs.items() if v[1]}


# Processes whose normalisation the fit can move appreciably. Used as the safety
# catch in the relevance cut: a large coherent variation on one of these is never
# dropped for being small relative to the *current* total, because the fit is
# free to scale the process up.
def processes_with_loose_normalisation(year):
    loose = set(FREE_NORM_PROCESSES)
    for name, (kappa, procs) in lnN_systematics(year).items():
        k = max(kappa) if isinstance(kappa, tuple) else kappa
        if k >= 1.15:
            loose.update(procs)
    return sorted(loose)


# --------------------------------------------------------------------------
# systematic metadata
# --------------------------------------------------------------------------

# How a shape variation relates statistically to the nominal:
#
#   "independent" : separate MC sample (CR1/CR2/tune_CP5)
#   "migration"   : same events, jets/MET recomputed, events move between bins
#   "correlated"  : same events, same kinematics, only the event weight changes
#
# For "independent", var(varied - nominal) = var_nom + var_syst is exact. For the
# other two the shared events cancel in the difference and the exact variance is
# NOT recoverable from the stored sumw2 of the two histograms -- it needs
# sum(w^2 (r-1)^2), which is only available at fill time (see COMBINE_FIXES.md).
#
# So the default here is deliberately conservative: use var_nom + var_syst for
# every class. That over-estimates the uncertainty on the shift, hence
# under-estimates its significance, so a variation is only ever dropped when it
# is compatible with zero even under an inflated error. The Barlow form for
# nested samples, |var_syst - var_nom|, is available but is not obviously better
# here: for a small weight variation it goes to zero and the significance
# diverges.

_ALT_SAMPLE_SYSTS = {"CR", "CR1", "CR2", "tune_CP5"}

_OBJECT_SYSTS_BASE = {
    "jer", "met",
    "jes_Absolute", "jes_BBEC1", "jes_FlavorQCD", "jes_RelativeBal",
}
_OBJECT_SYSTS_YEARLY = {
    "jes_Absolute", "jes_BBEC1", "jes_RelativeSample",
}


def correlation_class(syst, year):
    """Return 'independent', 'migration' or 'correlated' for a systematic."""
    if syst in _ALT_SAMPLE_SYSTS:
        return "independent"
    if syst in _OBJECT_SYSTS_BASE:
        return "migration"
    if syst.endswith(f"_{year}") and syst[: -len(f"_{year}")] in _OBJECT_SYSTS_YEARLY:
        return "migration"
    return "correlated"


def _empirical_ratio_noise(nom, leg, var_nom=None):
    """Robust per-bin variance of (leg - nom), from the high-frequency scatter.

    For a weight-based variation neither var_nom + var_leg nor |var_leg -
    var_nom| estimates var(leg - nom): both overestimate it by orders of
    magnitude (factors 2/eps^2 and 2/eps for a variation of relative size eps).
    Neither is recoverable from the stored sumw2 alone -- the exact quantity is
    sum(w^2 (r-1)^2), available only at fill time.

    What is recoverable is the scatter of the ratio r = leg/nom around a smooth
    curve. Second differences kill any linear trend, so the MAD of
    r[i+1] - 2 r[i] + r[i-1] measures noise and not structure. The sqrt(6)
    normalises the second difference of independent points to a per-point sigma.

    Two properties matter and were both wrong before:
      * the second difference must use ADJACENT bins. Taken on the compressed
        array it leaps over empty bins, and that step is structure, not noise.
      * the ratio noise scales as 1/sqrt(n_eff), so a single flat relative error
        is wrong. Each second difference is divided by its own 1/sqrt(n_eff)
        factor, which leaves one dimensionless constant c, and the per-bin
        variance is then  var_i = c^2 * nom_i^2 / n_eff_i.
        With n_eff_i = nom_i^2 / var_nom_i this is  var_i = c^2 * var_nom_i.
    """
    import numpy as _np

    m = nom > 0
    if m.sum() < 4:
        return None

    r = _np.zeros_like(nom, dtype=float)
    r[m] = leg[m] / nom[m]
    trip = m[:-2] & m[1:-1] & m[2:]
    if trip.sum() < 2:
        return None
    d2 = (r[2:] - 2.0 * r[1:-1] + r[:-2])[trip]

    if var_nom is None:
        # no statistics available: fall back to one flat relative error
        mad = _np.median(_np.abs(d2 - _np.median(d2)))
        c = 1.4826 * mad / _np.sqrt(6.0)
        if not _np.isfinite(c) or c <= 0:
            return None
        out = _np.zeros_like(nom, dtype=float)
        out[m] = (c * nom[m]) ** 2
        return out

    # relative statistical error of the nominal, per bin
    rel = _np.zeros_like(nom, dtype=float)
    ok = m & (var_nom > 0)
    rel[ok] = _np.sqrt(var_nom[ok]) / nom[ok]

    # the expected scale of one second difference, from the three bins in it
    scale = _np.sqrt(rel[:-2] ** 2 + 4.0 * rel[1:-1] ** 2 + rel[2:] ** 2)[trip]
    use = scale > 0
    if use.sum() < 2:
        return None
    z = d2[use] / scale[use]

    c = 1.4826 * _np.median(_np.abs(z - _np.median(z)))
    if not _np.isfinite(c) or c <= 0:
        return None

    out = _np.zeros_like(nom, dtype=float)
    out[ok] = (c ** 2) * var_nom[ok]
    return out


def shift_variance(nom, var_nom, leg, var_leg, cls):
    """Estimated variance of (varied - nominal) per bin.

    Exact for cls == 'independent'; for the same-event classes the empirical
    ratio-noise estimate is the only one not biased by orders of magnitude
    (see the note above).
    """
    if cls == "independent":
        return var_nom + var_leg
    emp = _empirical_ratio_noise(nom, leg, var_nom)
    if emp is not None:
        return emp
    return abs(var_leg - var_nom)  # fall back when there are too few bins

def raw_shift_variance(var_nom, var_leg):
    """Sumw2(nominal) + Sumw2(leg) -- the classic Combine chain's variance for
    (varied - nominal) (analysis/smoothing.py get_smoothed_scale_factor,
    total_var = nom_hist.variances() + var_hist.variances()).

    NOT statistically correct for a same-event weight variation ('migration' or
    'correlated' class): nominal and varied share the same events, so the sums
    do not cancel that correlation. It OVER-ESTIMATES var(leg - nominal) by a
    factor of order 2/eps^2 for a weight shift of relative size eps (eps ~ 1%
    here -> up to ~10^4x, see COMBINE_FIXES.md item 7). shift_variance() above
    is the fix: it uses the empirical ratio-noise estimate instead.

    Kept only so the least-squares rescale in the smoother can be run with the
    classic weighting, to A/B test against the old chain (--rescale-variance
    raw). Never use this as the default: it reintroduces a documented bias.
    """
    return var_nom + var_leg


# --------------------------------------------------------------------------
# nuisance groups
# --------------------------------------------------------------------------

# --------------------------------------------------------------------------
# leg-preprocessing policy
# --------------------------------------------------------------------------
# Applied to the histograms in analysis/prepareTensor.py BEFORE the significance
# tests, so the repaired / one-sided variation goes through the same drop, norm
# and smoothing decisions as any other systematic. Matching is by substring
# (mirror-up) or exact name (two-point), first match wins.

# flavTag: the stored down leg is unreliable, so rebuild it by mirroring the up
# leg about the nominal (down = nom^2 / up). This is what the old Combine chain
# did (analysis/smoothing.py forceSymm, up leg preferred).
MIRROR_UP = ["CMS_flavTag_"]


def is_mirror_up(syst):
    """True if the down leg must be rebuilt by mirroring the up leg."""
    return any(entry in syst for entry in MIRROR_UP)


def is_two_point(syst):
    """Alternative-sample (two-point) systematic: one real sample vs nominal.

    Treated one-sided -- the up leg is the alternative, the down leg is set to
    the nominal. It is still statistics-limited, so it stays eligible for
    smoothing and for conversion to a normalisation, but it is never smoothed
    into an antisymmetric shape it does not have.
    """
    return syst in _ALT_SAMPLE_SYSTS


def nuisance_groups(all_nuisances):
    """Reproduce the `group` lines of the Combine datacard."""
    allnp = sorted(all_nuisances)
    return {
        "systematics": allnp,
        "allbutflavor": [n for n in allnp if not n.startswith("CMS_flavTag_")],
    }
