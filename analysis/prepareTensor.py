#!/usr/bin/env python3
"""Build a rabbit input tensor from the merged Combine shapes file.

Same inputs as the Combine path: the merged <card>_shapes.root plus the model
definition in configs/model.py. Replaces ValidateDatacards -> simplifyDatacards
-> text2workspace with a single pass:

    0)   drop exact no-ops
    0.2) validate; apply the sentinel policy
    1)   Test A -- is the shift significant against the MC statistical error?
    2)   Test B -- is the shape significant against flat? if not, keep as norm
    3)   smooth the entries still tagged as shape
    3.5) re-validate after smoothing
    4)   drop what failed Test A or the relevance cut
    5)   symmetrise only same-sign and one-sided variations
    6)   write the tensor

A JSON report of every decision is written next to the tensor.
"""

import argparse
import csv
import json
import os
import re
from collections import defaultdict

import hist
import numpy as np
import uproot
from scipy import stats

from analysis.templateSmoothing import smooth_variation, classic_forcesymm_smooth
from configs import model as M

UP, DOWN = "Up", "Down"


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def _split_key(name, processes):
    """'ttbb_CMS_pileup_2024Up' -> ('ttbb', 'CMS_pileup_2024', 'Up')."""
    for leg in (UP, DOWN):
        if name.endswith(leg):
            stem = name[: -len(leg)]
            for proc in processes:
                if stem.startswith(proc + "_"):
                    return proc, stem[len(proc) + 1:], leg
            return None, None, None
    return (name, None, None) if name in processes else (None, None, None)


def read_shapes(path):
    f = uproot.open(path)
    procs = sorted(M.ALL_PROCESSES, key=len, reverse=True)

    edges, data, nominal = {}, {}, {}
    variations = defaultdict(dict)
    unknown = []

    for key in f.keys(cycle=False):
        if "/" not in key:
            continue
        cat, name = key.split("/", 1)
        h = f[key]
        vals, vars_ = h.values(), h.variances()
        if cat not in edges:
            edges[cat] = h.axis().edges()

        if name == "data_obs":
            data[cat] = [vals, vars_]
            continue

        proc, syst, leg = _split_key(name, procs)
        if proc is None:
            unknown.append(key)
        elif syst is None:
            nominal[(cat, proc)] = [vals, vars_]
        else:
            variations[(cat, proc, syst)][leg] = [vals, vars_]

    return edges, data, nominal, dict(variations), unknown


# --------------------------------------------------------------------------
# 0.2 validate
# --------------------------------------------------------------------------

def validate(edges, data, nominal, variations, unknown, args, report):
    errors, warnings = [], []

    file_cats = set(edges)
    model_cats = set(M.CATEGORIES)
    if file_cats != model_cats:
        if model_cats - file_cats:
            errors.append(f"categories missing from shapes file: {sorted(model_cats - file_cats)}")
        if file_cats - model_cats:
            errors.append(f"categories in shapes file not declared in model: {sorted(file_cats - model_cats)}")

    if unknown:
        errors.append(f"{len(unknown)} histograms could not be matched to a declared process, e.g. {unknown[:3]}")

    for cat in sorted(file_cats & model_cats):
        nbins = len(edges[cat]) - 1
        if cat not in data:
            errors.append(f"{cat}: no data_obs")
        for proc in M.ALL_PROCESSES:
            if (cat, proc) not in nominal:
                errors.append(f"{cat}: no nominal for declared process {proc}")
                continue
            if len(nominal[(cat, proc)][0]) != nbins:
                errors.append(f"{cat}/{proc}: {len(nominal[(cat, proc)][0])} bins, category has {nbins}")

    for (cat, proc), (v, e) in nominal.items():
        if not np.all(np.isfinite(v)) or not np.all(np.isfinite(e)):
            errors.append(f"{cat}/{proc}: non-finite bins in nominal")
        if (v < 0).any():
            errors.append(f"{cat}/{proc}: negative nominal bins")

    for (cat, proc, syst), legs in variations.items():
        for leg, (v, e) in legs.items():
            if not np.all(np.isfinite(v)) or not np.all(np.isfinite(e)):
                errors.append(f"{cat}/{proc}_{syst}{leg}: non-finite bins")

    # undo the 1e-6 placeholder Combine needs: set those bins back to 0
    n_sentinel_nom = n_sentinel_var = 0
    sentinel_bins = {}
    for (cat, proc), arrs in nominal.items():
        m = np.isclose(arrs[0], M.SENTINEL, rtol=1e-3, atol=0)
        if m.any():
            arrs[0] = np.where(m, 0.0, arrs[0])
            arrs[1] = np.where(m, 0.0, arrs[1])
            sentinel_bins[(cat, proc)] = m
            n_sentinel_nom += int(m.sum())
    for (cat, proc, syst), legs in variations.items():
        dead = sentinel_bins.get((cat, proc))
        for leg, arrs in legs.items():
            m = np.isclose(arrs[0], M.SENTINEL, rtol=1e-3, atol=0)
            if dead is not None:
                m = m | dead
            if m.any():
                arrs[0] = np.where(m, 0.0, arrs[0])
                arrs[1] = np.where(m, 0.0, arrs[1])
                n_sentinel_var += int(m.sum())
    warnings.append(f"restored {n_sentinel_nom} nominal and {n_sentinel_var} variation sentinel bins to 0")

    # total prediction must be positive wherever there is data
    for cat in sorted(file_cats & model_cats):
        if cat not in data:
            continue
        tot = np.zeros(len(edges[cat]) - 1)
        for proc in M.ALL_PROCESSES:
            if (cat, proc) in nominal:
                tot = tot + nominal[(cat, proc)][0]
        bad = (tot <= 0) & (data[cat][0] > 0)
        if bad.any():
            errors.append(f"{cat}: zero total prediction in bin(s) {list(np.where(bad)[0] + 1)} where data > 0")

    # warnings only
    low_neff = 0
    for (cat, proc), (v, e) in nominal.items():
        m = (v > 0) & (e > 0)
        if m.any():
            low_neff += int((v[m] ** 2 / e[m] < args.low_neff_warn).sum())
    if low_neff:
        warnings.append(f"{low_neff} nominal bins with n_eff < {args.low_neff_warn}")

    clamped = 0
    for legs in variations.values():
        for arrs in legs.values():
            m = arrs[0] < 0
            if m.any():
                arrs[0] = np.where(m, 0.0, arrs[0])
                clamped += int(m.sum())
    if clamped:
        warnings.append(f"{clamped} negative variation bins clamped to 0")

    report["validation"] = {"errors": errors, "warnings": warnings,
                            "sentinel_nominal_bins": n_sentinel_nom,
                            "sentinel_variation_bins": n_sentinel_var}
    for w in warnings:
        print(f"  [warn]  {w}")
    if errors:
        for e in errors:
            print(f"  [ERROR] {e}")
        raise SystemExit(f"validation failed with {len(errors)} error(s)")


# --------------------------------------------------------------------------
# 0.3 process yield scaling (togglable; off at factor 1.0)
# --------------------------------------------------------------------------

def scale_processes(nominal, variations, args, report):
    """Multiply the yield of the tt+bb DPS processes by ``args.dps_scale``.

    Emulates a change to the cross section of ``ttbb-dps``, ``ttbj-dps`` and
    ``tt2b-dps`` in the simulation: every nominal template and every systematic
    variation leg of those processes is scaled by the same constant. Bin
    variances are scaled by ``factor**2`` (a larger event weight, same MC
    statistics), so the relative MC-stat error and every relative systematic on
    these processes are preserved -- only their absolute normalisation moves.
    The related ``norm_ttbb-dps`` lnN is a relative kappa and needs no change.

    Runs after validate() (so sentinel bins are already restored to 0 and are
    not perturbed) and before every downstream test, so the scaled templates go
    through no-op detection, classification, smoothing and the relevance cut
    exactly like any other input. Factor 1.0 is a no-op and leaves the tensor
    bit-for-bit identical to a run without the flag.
    """
    f = args.dps_scale
    procs = set(M.TT_COMPONENTS_BBDPS)
    if f == 1.0:
        report["process_scaling"] = {"factor": 1.0, "processes": [],
                                     "note": "disabled (--dps-scale 1.0)"}
        print("  process yield scaling disabled (--dps-scale 1.0)")
        return
    if f <= 0.0:
        raise SystemExit(f"--dps-scale must be positive, got {f}")

    n_nom = n_leg = 0
    for (cat, proc), arrs in nominal.items():
        if proc not in procs:
            continue
        arrs[0] = arrs[0] * f
        arrs[1] = arrs[1] * f * f
        n_nom += 1
    for (cat, proc, syst), legs in variations.items():
        if proc not in procs:
            continue
        for arrs in legs.values():
            arrs[0] = arrs[0] * f
            arrs[1] = arrs[1] * f * f
            n_leg += 1

    report["process_scaling"] = {"factor": f, "processes": sorted(procs),
                                 "nominal_templates_scaled": n_nom,
                                 "variation_legs_scaled": n_leg}
    print(f"  scaled {sorted(procs)} by {f}: "
          f"{n_nom} nominal templates, {n_leg} variation legs")


# --------------------------------------------------------------------------
# 0 no-ops
# --------------------------------------------------------------------------

def drop_noops(nominal, variations, report, dec):
    dropped = []
    for key in list(variations):
        cat, proc, syst = key
        legs = variations[key]
        if UP not in legs or DOWN not in legs:
            continue
        nom = nominal[(cat, proc)][0]
        if np.allclose(legs[UP][0], nom) and np.allclose(legs[DOWN][0], nom):
            dropped.append(key)
            dec[key].update(decision="dropped_noop",
                            reason="both legs identical to nominal in every bin")
            del variations[key]
    report["noops_dropped"] = len(dropped)
    print(f"  dropped {len(dropped)} exact no-op entries")
    return dropped


# --------------------------------------------------------------------------
# 0.5 leg preprocessing (before the significance tests)
# --------------------------------------------------------------------------

def preprocess_legs(nominal, variations, args, report, dec):
    """Repair / normalise legs BEFORE classify, smooth and drop, so the result
    goes through the same decisions as any other systematic.

    * mirror-up (flavTag): the stored down leg is unreliable, so rebuild it by
      mirroring the up leg about the nominal, down = nom^2 / up where up > 0.
    * two-point (alt-sample): one alternative sample, so keep the up leg as the
      single variation and set the down leg to the nominal. It stays eligible
      for smoothing / norm conversion but is never forced into an antisymmetric
      shape it does not have.

    Skipped with --no-symmetrisation, so the verbatim baseline is unchanged.
    Default is --no-flavtag-mirror's old behaviour (stored flavTag legs kept
    as-is); pass --flavtag-mirror to rebuild the down leg instead
    (down = nom^2/up), the old default, kept for comparison/debugging.
    """
    if args.no_symmetrisation:
        report["leg_preprocess"] = "skipped (--no-symmetrisation)"
        print("  leg preprocess skipped (--no-symmetrisation)")
        return
    n_mirror = n_two = n_kept = 0
    for (cat, proc, syst), legs in variations.items():
        nom, nomvar = nominal[(cat, proc)]
        if M.is_mirror_up(syst):
            if not args.flavtag_mirror:
                # keep whatever the shapes file stored, the same thing the
                # Combine noFlavTagSymm card does
                dec[(cat, proc, syst)]["leg_preprocess"] = "flavTag stored legs (default)"
                n_kept += 1
                continue
            if args.flavtag_classic:
                # deferred to smooth_flavtag_classic(), which runs on the
                # untouched (possibly broken) stored legs
                dec[(cat, proc, syst)]["leg_preprocess"] = "flavTag classic (deferred to smoothing)"
                continue
            if UP not in legs:
                continue
            up = legs[UP][0]
            ok = (nom > 0) & (up > 0)
            dn = nom.copy()
            dn[ok] = nom[ok] * nom[ok] / up[ok]
            legs[DOWN] = [dn, nomvar.copy()]
            dec[(cat, proc, syst)]["leg_preprocess"] = "mirror-up (down = nom^2/up)"
            n_mirror += 1
        elif M.is_two_point(syst):
            if UP not in legs and DOWN in legs:
                legs[UP] = legs.pop(DOWN)
            legs[DOWN] = [nom.copy(), nomvar.copy()]
            dec[(cat, proc, syst)]["leg_preprocess"] = "two-point one-sided (down = nominal)"
            n_two += 1
    report["leg_preprocess"] = {"mirror_up": n_mirror, "two_point_one_sided": n_two,
                                "flavtag_stored_legs": n_kept}
    print(f"  leg preprocess: {n_mirror} flavTag mirror-up, {n_two} two-point one-sided"
          + (f", {n_kept} flavTag stored legs kept (default)" if n_kept else ""))


def smooth_flavtag_classic(edges, nominal, variations, args, report, dec):
    """Test-only alternative to the mirror-up in preprocess_legs(): apply the
    classic analysis/smoothing.py forceSymm treatment to every flavTag entry
    with both legs, UNCONDITIONALLY (matching the old chain's
    bypassSmoothing=True -- this never waits for classify()'s shape/norm
    decision). Only runs when --flavtag-classic is set.
    """
    centers = {c: 0.5 * (e[1:] + e[:-1]) for c, e in edges.items()}
    n_done = n_skipped = 0
    for key, legs in variations.items():
        cat, proc, syst = key
        if not M.is_mirror_up(syst):
            continue
        if UP not in legs or DOWN not in legs:
            n_skipped += 1
            continue
        nom, var_nom = nominal[(cat, proc)]
        res = classic_forcesymm_smooth(nom, legs[UP][0], legs[DOWN][0], var_nom,
                                       legs[UP][1], legs[DOWN][1], centers[cat],
                                       lowess_frac=args.lowess_frac)
        if res is None:
            n_skipped += 1
            dec[key]["smoothed"] = "no (classic: too few usable bins)"
            continue
        legs[UP][0] = res.up
        legs[DOWN][0] = res.down
        dec[key]["smoothed"] = "yes (classic forceSymm)"
        n_done += 1
    report["flavtag_classic"] = {"smoothed": n_done, "skipped": n_skipped}
    print(f"  flavTag classic forceSymm: {n_done} entries smoothed, {n_skipped} skipped")


def _smooth_tags(tags, args):
    """tags for the generic smooth() call: when --flavtag-classic already
    smoothed the flavTag entries (unconditionally, before this runs), mask
    them out here so smooth() does not smooth them a second time."""
    if not args.flavtag_classic:
        return tags
    return {k: (v if not M.is_mirror_up(k[2]) else "done (flavtag classic)")
            for k, v in tags.items()}


# --------------------------------------------------------------------------
# 1 + 2 significance tests
# --------------------------------------------------------------------------

def _leg_pvalues(nom, var_nom, leg, var_leg, cls):
    """(p_shift, p_shape) for one leg."""
    s2 = M.shift_variance(nom, var_nom, leg, var_leg, cls)
    # A bin with no nominal content carries no information about a shift; drop it
    # so an empty bin cannot dominate the chi2.
    usable = (s2 > 0) & np.isfinite(s2) & (nom > 0)
    if usable.sum() == 0:
        return 1.0, 1.0          # no information -> not significant
    d = (leg - nom)[usable]
    chi2_shift = float(np.sum(d * d / s2[usable]))
    p_shift = float(stats.chi2.sf(chi2_shift, usable.sum()))

    # k must come from the same bins the chi2 runs over, otherwise the
    # constraint sum(d_shape) = 0 does not hold on that set and ndf = n-1 is
    # wrong.
    ntot, ltot = nom[usable].sum(), leg[usable].sum()
    if ltot <= 0 or ntot <= 0 or usable.sum() < 2:
        return p_shift, 1.0
    k = ntot / ltot
    d_shape = (leg * k - nom)[usable]
    chi2_shape = float(np.sum(d_shape * d_shape / (s2[usable] * k * k)))
    p_shape = float(stats.chi2.sf(chi2_shape, usable.sum() - 1))
    return p_shift, p_shape


def classify(nominal, variations, args, report, dec, totals=None):
    """Tag every entry as 'shape' or 'norm'.

    A shape is replaced by a pure normalisation only when BOTH conditions hold:
    the flat fit is acceptable, and the shape that this throws away is small
    against the same yardstick that decides whether a systematic is worth
    keeping at all (--thr-total).

    The p-value alone is not enough. The flat test has little power, so a shape
    worth 0.5 % of the total prediction can pass it with p = 0.5. Keeping a
    systematic that moves the total by 0.1 %, and at the same time discarding a
    shape worth 0.5 %, is not consistent.
    """
    tags, pvals = {}, {}
    for key, legs in variations.items():
        cat, proc, syst = key
        nom, var_nom = nominal[(cat, proc)]
        cls = M.correlation_class(syst, args.year)

        ps, psh = [], []
        for leg in (UP, DOWN):
            if leg not in legs:
                continue
            a, b = _leg_pvalues(nom, var_nom, legs[leg][0], legs[leg][1], cls)
            ps.append(a)
            psh.append(b)

        # Taking min() over the two legs is two tries at one threshold, so the
        # true false-positive rate is about 2*alpha. Fisher's rule combines them
        # into one p-value with the right calibration.
        def _combine(vals):
            vals = [min(max(v, 1e-300), 1.0) for v in vals]
            if not vals:
                return 1.0
            if len(vals) == 1:
                return vals[0]
            x = -2.0 * sum(np.log(v) for v in vals)
            return float(stats.chi2.sf(x, 2 * len(vals)))

        p_shift = _combine(ps)
        p_shape = _combine(psh)
        pvals[key] = (p_shift, p_shape)

        # how much shape the conversion would throw away, against the total
        disc = 0.0
        tot = None if totals is None else totals.get(cat)
        if tot is not None:
            ok = tot > 0
            for leg in legs.values():
                tn, tl = nom.sum(), leg[0].sum()
                if tn <= 0 or tl <= 0:
                    continue
                k = tl / tn
                d = np.abs(leg[0] - nom * k)
                if ok.any():
                    disc = max(disc, float(np.max(d[ok] / tot[ok])))
        dec[key]["discarded_shape"] = disc

        if p_shape > args.p_shape and disc < args.thr_total:
            tags[key] = "norm"
            why = (f"flat-line fit acceptable (p_shape={p_shape:.3g} > {args.p_shape}) "
                   f"and the discarded shape is {disc:.2e} of the total, below "
                   f"thr-total={args.thr_total}")
        elif p_shape > args.p_shape:
            tags[key] = "shape"
            why = (f"flat-line fit acceptable (p_shape={p_shape:.3g}) but the shape "
                   f"it would discard is {disc:.2e} of the total, at or above "
                   f"thr-total={args.thr_total}; the flat test has little power")
        else:
            tags[key] = "shape"
            why = f"shape needed: p_shape={p_shape:.3g} <= {args.p_shape}"
        dec[key].update(correlation_class=cls, p_shift=p_shift, p_shape=p_shape,
                        tag=tags[key], tag_reason=why,
                        one_sided=int(len(legs) < 2))

    counts = {t: sum(1 for v in tags.values() if v == t) for t in ("shape", "norm")}
    report["classification"] = counts
    print(f"  shape {counts['shape']}   norm-only {counts['norm']}")
    return tags, pvals


# --------------------------------------------------------------------------
# 3 smoothing
# --------------------------------------------------------------------------

def smooth(edges, nominal, variations, tags, args, report, dec):
    centers = {c: 0.5 * (e[1:] + e[:-1]) for c, e in edges.items()}
    n_done = n_skipped = 0
    notes = []
    for key, tag in tags.items():
        if tag != "shape":
            continue
        cat, proc, syst = key
        legs = variations[key]
        if UP not in legs or DOWN not in legs:
            n_skipped += 1
            dec[key]["smoothed"] = "no (one-sided)"
            continue
        nom, var_nom = nominal[(cat, proc)]
        cls = M.correlation_class(syst, args.year)
        if args.rescale_variance == "raw":
            s2_up = M.raw_shift_variance(var_nom, legs[UP][1])
            s2_dn = M.raw_shift_variance(var_nom, legs[DOWN][1])
        else:
            s2_up = M.shift_variance(nom, var_nom, legs[UP][0], legs[UP][1], cls)
            s2_dn = M.shift_variance(nom, var_nom, legs[DOWN][0], legs[DOWN][1], cls)

        res = smooth_variation(nom, legs[UP][0], legs[DOWN][0], s2_up, s2_dn,
                               centers[cat], smoothing_factor=args.smoothing_factor,
                               antisymmetrise=not M.is_two_point(syst),
                               method=args.smoothing_method, lowess_frac=args.lowess_frac)
        if res is None:
            n_skipped += 1
            dec[key]["smoothed"] = "no (smoothing not applicable)"
            continue
        legs[UP][0] = res.up
        legs[DOWN][0] = res.down
        n_done += 1
        dec[key]["smoothed"] = "yes"
        if res.note:
            notes.append(f"{cat}/{proc}_{syst}: {res.note}")
            dec[key]["smoothing_note"] = res.note

    report["smoothing"] = {"smoothed": n_done, "skipped": n_skipped, "notes": notes[:50]}
    print(f"  smoothed {n_done} entries, skipped {n_skipped}")
    if notes:
        print(f"  {len(notes)} smoothing notes (first in report)")


def revalidate(nominal, variations, report):
    problems = []
    for (cat, proc, syst), legs in variations.items():
        for leg, (v, _e) in legs.items():
            if not np.all(np.isfinite(v)):
                problems.append(f"{cat}/{proc}_{syst}{leg}: non-finite after smoothing")
            elif (v < 0).any():
                problems.append(f"{cat}/{proc}_{syst}{leg}: negative after smoothing")
    report["post_smoothing_problems"] = problems
    if problems:
        for p in problems[:10]:
            print(f"  [ERROR] {p}")
        raise SystemExit(f"post-smoothing validation failed ({len(problems)} problems)")
    print("  post-smoothing validation clean")


# --------------------------------------------------------------------------
# 4 relevance
# --------------------------------------------------------------------------

def apply_drops(edges, nominal, variations, tags, args, report, dec):
    """Keep or drop a whole nuisance parameter, on its coherent effect.

    A nuisance moves every process at the same time, and the shifts add. Judging
    one (category, process) entry against the category total therefore drops a
    systematic that moves all 18 processes by 0.05 % each, although its real
    effect on the prediction is close to 1 %. The comparison now uses the SIGNED
    sum over processes, which is the physical shift of the prediction, and the
    decision is taken once per systematic over all categories. "Drop" then means
    "this parameter does nothing anywhere", which is the only defensible drop.
    """
    loose = set(M.processes_with_loose_normalisation(args.year))
    totals = {}
    for cat in edges:
        tot = np.zeros(len(edges[cat]) - 1)
        for proc in M.ALL_PROCESSES:
            if (cat, proc) in nominal:
                tot = tot + nominal[(cat, proc)][0]
        totals[cat] = tot

    # coherent shift of the prediction, per (category, systematic) and leg
    coherent = defaultdict(lambda: defaultdict(float))   # syst -> cat -> rel
    rel_proc_of = {}
    summed = defaultdict(dict)                           # (cat, syst) -> leg -> array
    for (cat, proc, syst), legs in variations.items():
        nom = nominal[(cat, proc)][0]
        m = nom > 0
        rp = 0.0
        for leg, arr in legs.items():
            d = arr[0] - nom
            acc = summed[(cat, syst)].get(leg)
            summed[(cat, syst)][leg] = d if acc is None else acc + d
            if m.any():
                rp = max(rp, float(np.max(np.abs(d[m]) / nom[m])))
        rel_proc_of[(cat, proc, syst)] = rp

    for (cat, syst), legs in summed.items():
        tot = totals[cat]
        ok = tot > 0
        best = 0.0
        for d in legs.values():
            if ok.any():
                best = max(best, float(np.max(np.abs(d[ok]) / tot[ok])))
        coherent[syst][cat] = best

    rel_syst = {sy: max(c.values()) for sy, c in coherent.items()}

    # largest per-process effect of each systematic, over all (cat, proc)
    rel_syst_proc = defaultdict(float)
    for (cat, proc, syst), rp in rel_proc_of.items():
        rel_syst_proc[syst] = max(rel_syst_proc[syst], rp)

    # a large variation on a freely floating process protects its systematic
    protected = set()
    for key, rp in rel_proc_of.items():
        if key[1] in loose and rp > args.safety_proc:
            protected.add(key[2])

    keep = {sy for sy in rel_syst
            if rel_syst[sy] >= args.thr_total
            or rel_syst_proc[sy] >= args.thr_proc
            or sy in protected}

    dropped_irrel, kept_by_catch, kept_by_proc = [], [], []
    for key in list(variations):
        cat, proc, syst = key
        rel_tot = coherent[syst].get(cat, 0.0)
        dec[key].update(rel_total=rel_tot, rel_process=rel_proc_of[key],
                        rel_systematic=rel_syst[syst])
        if syst in keep:
            if rel_syst[syst] >= args.thr_total:
                dec[key].update(
                    decision=tags[key],
                    reason=f"kept: coherent shift of this systematic is "
                           f"{rel_syst[syst]:.2e} >= thr-total={args.thr_total} "
                           f"({dec[key].get('tag_reason', '')})")
            elif rel_syst_proc[syst] >= args.thr_proc:
                kept_by_proc.append(key)
                dec[key].update(
                    decision=tags[key],
                    reason=f"kept: this systematic reaches {rel_syst_proc[syst]:.2e} "
                           f"on a single process >= thr-proc={args.thr_proc}, though "
                           f"its coherent total is only {rel_syst[syst]:.2e} "
                           f"({dec[key].get('tag_reason', '')})")
            else:
                kept_by_catch.append(key)
                dec[key].update(decision=tags[key], kept_by_safety_catch=1,
                                reason=f"coherent shift {rel_syst[syst]:.2e} < "
                                       f"{args.thr_total}, but a freely floating "
                                       f"process moves by more than {args.safety_proc}")
        else:
            dropped_irrel.append(key)
            dec[key].update(
                decision="dropped_irrelevant",
                reason=f"the whole systematic moves the total prediction by at "
                       f"most {rel_syst[syst]:.2e} (thr-total={args.thr_total}) and any "
                       f"single process by at most {rel_syst_proc[syst]:.2e} "
                       f"(thr-proc={args.thr_proc})")
            del variations[key]

    n_syst_kept = len({k[2] for k in variations})
    report["drops"] = {
        "irrelevant": len(dropped_irrel),
        "kept_by_process": len(kept_by_proc),
        "kept_by_safety_catch": len(kept_by_catch),
        "remaining": len(variations),
        "systematics_kept": n_syst_kept,
        "systematics_dropped": len(rel_syst) - n_syst_kept,
    }
    print(f"  dropped {len(dropped_irrel)} irrelevant, kept {len(kept_by_proc)} by per-process "
          f"-> {len(variations)} entries, {n_syst_kept} of {len(rel_syst)} systematics kept")


# --------------------------------------------------------------------------
# 5 symmetrisation policy
# --------------------------------------------------------------------------

def symmetrisation(nominal, variations, report, dec, no_symmetrisation=False,
                   flavtag_mirror=False):
    """Rabbit `symmetrize` policy.

    All leg shaping is done on the histograms in preprocess_legs (step 0.5):
    flavTag is mirror-up, two-point is one-sided, everything else is verbatim.
    So rabbit receives the final legs and interpolates up/down directly and the
    policy is None for every entry. This only records, per systematic, which
    preprocessing was applied.
    """
    policy, syst_rows, seen = {}, [], {}
    for key in variations:
        policy[key] = None
        dec[key]["symmetrize"] = "None"
    for (cat, proc, syst) in variations:
        if syst in seen:
            continue
        if no_symmetrisation:
            why = "disabled with --no-symmetrisation"
        elif M.is_mirror_up(syst):
            why = ("flavTag: down rebuilt by mirroring up (step 0.5)"
                   if flavtag_mirror
                   else "flavTag: stored legs kept (default)")
        elif M.is_two_point(syst):
            why = "two-point: one-sided, down set to nominal (step 0.5)"
        else:
            why = "asymmetric; rabbit interpolates up and down directly"
        seen[syst] = why
        syst_rows.append({"systematic": syst, "symmetrize": "None", "reason": why})
    report["symmetrisation"] = {"policy": "None (histograms preshaped in step 0.5)",
                                "systematics": len(seen)}
    print(f"  rabbit symmetrize=None for {len(policy)} entries "
          f"(shaping done in [0.5]); {len(seen)} systematics")
    return policy, syst_rows


# --------------------------------------------------------------------------
# 6 write
# --------------------------------------------------------------------------

def _mkhist(edges, values, variances, name):
    h = hist.Hist(hist.axis.Variable(edges, name=name, overflow=False, underflow=False),
                  storage=hist.storage.Weight())
    view = h.view()
    view["value"] = values
    view["variance"] = variances
    return h


def write_tensor(edges, data, nominal, variations, tags, policy, args, report, dec):
    from rabbit import tensorwriter

    writer = tensorwriter.TensorWriter(
        sparse=args.sparse,
        clip_syst_variations=args.clip_syst_variations,
        zero_syst_low_neff=args.zero_syst_low_neff,
    )

    masked = [c for c in M.CATEGORIES if args.mask and re.search(args.mask, c)]
    if masked:
        print(f"  masked channels: {masked}")

    for cat, dist in M.CATEGORIES.items():
        axis = hist.axis.Variable(edges[cat], name=dist, overflow=False, underflow=False)
        writer.add_channel([axis], cat, masked=cat in masked)
        if cat not in masked:
            writer.add_data(_mkhist(edges[cat], data[cat][0], data[cat][1], dist), cat)
        for proc in M.ALL_PROCESSES:
            v, e = nominal[(cat, proc)]
            writer.add_process(_mkhist(edges[cat], v, e, dist), proc, cat,
                               signal=proc in M.SIGNAL)

    all_nuisances = sorted({k[2] for k in variations} | set(M.lnN_systematics(args.year)))
    groups = M.nuisance_groups(all_nuisances)
    groups_of = defaultdict(list)
    for gname, members in groups.items():
        for n in members:
            groups_of[n].append(gname)

    n_shape = n_norm = 0
    for (cat, proc, syst), legs in variations.items():
        dist = M.CATEGORIES[cat]
        gr = groups_of.get(syst, [syst])
        if tags[(cat, proc, syst)] == "norm":
            nomv = nominal[(cat, proc)][0]
            nomvar = nominal[(cat, proc)][1]
            tot = nomv.sum()
            if tot <= 0:
                dec[(cat, proc, syst)]["written"] = "NOT written (zero nominal yield)"
                continue
            ku = legs[UP][0].sum() / tot if UP in legs else 1.0
            kd = legs[DOWN][0].sum() / tot if DOWN in legs else 1.0
            # Written as a flat shape variation rather than via
            # add_norm_systematic: mathematically identical to an asymmetric
            # lnN, but the asymmetric branch of add_norm_systematic passes the
            # process *list* into _compute_asym_syst instead of the current
            # process, so any symmetrize other than "average" raises there.
            hu = _mkhist(edges[cat], nomv * ku, nomvar, dist)
            hd = _mkhist(edges[cat], nomv * kd, nomvar, dist)
            writer.add_systematic([hu, hd], syst, proc, cat,
                                  symmetrize=policy[(cat, proc, syst)], groups=gr)
            dec[(cat, proc, syst)].update(written="norm (flat lnN-like)",
                                          kappa_up=ku, kappa_down=kd)
            n_norm += 1
        else:
            if UP not in legs or DOWN not in legs:
                # A one-sided entry used to be dropped here while the table said
                # it was kept. Mirror the single leg about the nominal instead,
                # which is what the norm branch already does.
                one = legs[UP] if UP in legs else legs[DOWN]
                nomv, nomvar = nominal[(cat, proc)]
                mirrored = np.where(nomv > 0, nomv * nomv / np.where(one[0] > 0, one[0], 1.0), nomv)
                hu = _mkhist(edges[cat], one[0], one[1], dist)
                hd = _mkhist(edges[cat], mirrored, nomvar, dist)
                writer.add_systematic([hu, hd], syst, proc, cat,
                                      symmetrize=None, groups=gr)
                dec[(cat, proc, syst)]["written"] = "shape (one-sided, mirrored)"
                n_shape += 1
                continue
            hu = _mkhist(edges[cat], legs[UP][0], legs[UP][1], dist)
            hd = _mkhist(edges[cat], legs[DOWN][0], legs[DOWN][1], dist)
            writer.add_systematic([hu, hd], syst, proc, cat,
                                  symmetrize=policy[(cat, proc, syst)], groups=gr)
            dec[(cat, proc, syst)]["written"] = "shape"
            n_shape += 1

    for name, (kappa, procs) in M.lnN_systematics(args.year).items():
        # one kappa per process: add_norm_systematic zips process against
        # uncertainty, so a scalar with a process list silently applies to the
        # first process only. A tuple is an asymmetric (up, down) lnN.
        unc = [tuple(kappa) if isinstance(kappa, tuple) else kappa] * len(procs)
        for cat in M.CATEGORIES:
            writer.add_norm_systematic(name, procs, cat, unc,
                                       groups=groups_of.get(name, [name]))

    report["written"] = {"shape_systematics": n_shape, "norm_systematics": n_norm,
                         "lnN": len(M.lnN_systematics(args.year)), "masked_channels": masked}
    print(f"  wrote {n_shape} shape + {n_norm} norm systematics + {len(M.lnN_systematics(args.year))} lnN")

    os.makedirs(args.output, exist_ok=True)
    writer.write(outfolder=args.output, outfilename=args.outname)
    return os.path.join(args.output, args.outname + ".hdf5")


# --------------------------------------------------------------------------
# decision tables
# --------------------------------------------------------------------------

COLUMNS = ["category", "process", "systematic", "correlation_class", "decision",
           "tag", "written", "p_shift", "p_shape", "rel_total", "rel_process",
           "rel_systematic", "discarded_shape", "smoothed", "symmetrize", "leg_preprocess", "one_sided",
           "kept_by_safety_catch",
           "kappa_up", "kappa_down", "smoothing_note", "reason"]


def write_tables(dec, syst_rows, args):
    base = os.path.join(args.output, args.outname)

    rows = []
    for (cat, proc, syst), d in sorted(dec.items()):
        r = {"category": cat, "process": proc, "systematic": syst}
        r.update(d)
        rows.append({c: r.get(c, "") for c in COLUMNS})
    with open(base + "_decisions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS)
        w.writeheader()
        w.writerows(rows)

    with open(base + "_symmetrisation.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(syst_rows[0]) if syst_rows else ["systematic"])
        w.writeheader()
        w.writerows(syst_rows)

    summary = defaultdict(int)
    for r in rows:
        summary[r["decision"] or "unknown"] += 1
    print("  decision census: " + "  ".join(f"{k}={v}" for k, v in sorted(summary.items())))
    print(f"  tables : {base}_decisions.csv , {base}_symmetrisation.csv")
    return rows, dict(summary)


# --------------------------------------------------------------------------

def make_parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("shapes", help="merged Combine shapes ROOT file")
    p.add_argument("-o", "--output", default="./", help="output directory")
    p.add_argument("--outname", default="Vcb_tensor", help="output file name without extension")
    p.add_argument("--year", default="2024")

    p.add_argument("--dps-scale", type=float, default=1.0,
                   help="scale the yield of the tt+bb DPS processes "
                        f"({', '.join(M.TT_COMPONENTS_BBDPS)}) -- nominal AND every "
                        "systematic variation leg -- by this factor, emulating a "
                        "change to their cross section in the simulation. 1.0 "
                        "(default) is a no-op. MC-stat variance scales as "
                        "factor**2, so the relative MC-stat error and every "
                        "relative systematic on these processes are unchanged; "
                        "only their absolute normalisation moves.")
    p.add_argument("--no-dps", action="store_true",
                   help="the shapes file's ttbar simulation is 5FS-only, so it has "
                        f"no separate DPS components ({', '.join(M.TT_COMPONENTS_BBDPS)}). "
                        "Drops them from the declared process list (and norm_ttbb-dps) "
                        "before validation, instead of failing with 'no nominal for "
                        "declared process ttbb-dps'.")

    p.add_argument("--p-shape", type=float, default=0.05,
                   help="treat as normalisation only if the shape is flat with p > this")
    p.add_argument("--thr-total", type=float, default=1e-3,
                   help="relevance: max |variation| relative to the total prediction in a bin")
    p.add_argument("--thr-proc", type=float, default=1e-3,
                   help="relevance: keep a systematic if its |variation| relative to a "
                        "single process reaches this in any bin, even when the coherent "
                        "effect on the total is below --thr-total (0 disables)")
    p.add_argument("--safety-proc", type=float, default=0.30,
                   help="never drop a variation this large relative to a freely-floating process")
    p.add_argument("--smoothing-factor", type=float, default=1.0)
    p.add_argument("--low-neff-warn", type=float, default=1.0)

    p.add_argument("--mask", default=None, help="regex; matching channels are added as masked")
    p.add_argument("--sparse", action="store_true")
    p.add_argument("--clip-syst-variations", type=float, default=10.0,
                   help="clip per-bin variation factors to [1/x, x]; 0 disables")
    p.add_argument("--zero-syst-low-neff", type=float, default=1.0,
                   help="zero systematic levers in bins with n_eff below this; 0 disables")
    p.add_argument("--no-smoothing", action="store_true")
    p.add_argument("--smoothing-method", choices=["spline", "lowess"], default="lowess",
                   help="shape smoother: 'lowess' (locally-weighted regression, as in "
                        "analysis/smoothing.py) or 'spline' (local cubic, chi2/ndf~1)")
    p.add_argument("--lowess-frac", type=float, default=0.9,
                   help="fraction of points in each local lowess fit (only used "
                        "with --smoothing-method lowess)")
    p.add_argument("--rescale-variance", choices=["model", "raw"], default="model",
                   help="variance used to weight the least-squares rescale onto each "
                        "leg: 'model' (shift_variance, the corrected estimate) or "
                        "'raw' (Sumw2(nominal)+Sumw2(leg), the classic chain's weight "
                        "-- NOT statistically correct for a same-event weight "
                        "systematic, see raw_shift_variance() in configs/model.py; "
                        "test-only, never the default)")
    p.add_argument("--no-symmetrisation", action="store_true",
                   help="keep every variation asymmetric, ignoring configs/model.py")
    p.add_argument("--flavtag-mirror", action="store_true",
                   help="rebuild the flavTag down leg as down = nom^2/up instead of "
                        "keeping the stored value (the old default; off by default now, "
                        "matching the Combine noFlavTagSymm treatment). Only affects "
                        "configs/model.py MIRROR_UP entries; two-point handling is "
                        "untouched (unlike --no-symmetrisation).")
    p.add_argument("--legacy-order", action="store_true",
                   help="smooth ALL two-sided entries unconditionally, before "
                        "classify()/apply_drops() -- matches the old chain's order "
                        "(smoothing.py runs before ValidateCards/simplifyDatacards, "
                        "so classification never sees the raw shape). Test only: "
                        "default order classifies first and smooths only the "
                        "entries tagged shape.")
    p.add_argument("--flavtag-classic", action="store_true",
                   help="reproduce analysis/smoothing.py's exact flavTag forceSymm "
                        "logic (classic_forcesymm_smooth(): fits a scale from EACH "
                        "leg against its own raw data, prefers up if positive else "
                        "falls back to down; absolute-space, raw-Sumw2-weighted) "
                        "instead of the deterministic mirror-up "
                        "(configs.model.is_mirror_up). Test only, not the default.")

    p.add_argument("--plots", nargs="?", const="AUTO", default=None,
                   metavar="DIR",
                   help="write one diagnostic PNG per (category, process, systematic) "
                        "entry under DIR/<decision>/<category>/. With no value the "
                        "directory is <output>/<outname>_plots")
    p.add_argument("--plot-what", choices=["kept", "all"], default="all",
                   help="'kept' skips the entries that were dropped")
    p.add_argument("--plot-max", type=int, default=0, help="0 = no limit")
    p.add_argument("--plot-jobs", type=int, default=os.cpu_count(),
                   help="parallel plot workers")
    p.add_argument("--plot-channels", default=None,
                   help="regex; only plot categories that match (plotting all "
                        "8 categories is ~16000 figures)")
    return p


def main():
    args = make_parser().parse_args()
    if args.no_dps and args.dps_scale != 1.0:
        raise SystemExit("--no-dps and --dps-scale are mutually exclusive: there "
                          "is no DPS process left to scale once it is dropped")
    if args.no_dps:
        M.disable_dps()

    report = {"input": os.path.abspath(args.shapes), "settings": vars(args)}
    dec = defaultdict(dict)

    print(f"reading {args.shapes}")
    edges, data, nominal, variations, unknown = read_shapes(args.shapes)
    print(f"  {len(edges)} categories, {len(nominal)} nominal templates, {len(variations)} variation entries")

    print("[0.2] validate")
    validate(edges, data, nominal, variations, unknown, args, report)

    print("[0.3] process yield scaling")
    scale_processes(nominal, variations, args, report)

    print("[0] no-ops")
    drop_noops(nominal, variations, report, dec)

    # pre-smoothing snapshot, so the plots can show raw vs written
    raw_all = ({k: {leg: v[0].copy() for leg, v in legs.items()}
                for k, legs in variations.items()} if args.plots else {})

    print("[0.5] leg preprocess (flavTag mirror-up / two-point one-sided)")
    preprocess_legs(nominal, variations, args, report, dec)

    if args.flavtag_classic:
        print("[0.6] flavTag classic forceSymm (test only, replaces mirror-up)")
        smooth_flavtag_classic(edges, nominal, variations, args, report, dec)

    totals = {}
    for cat in edges:
        t = np.zeros(len(edges[cat]) - 1)
        for proc in M.ALL_PROCESSES:
            if (cat, proc) in nominal:
                t = t + nominal[(cat, proc)][0]
        totals[cat] = t

    if args.legacy_order and not args.no_smoothing:
        print("[3] smoothing ALL entries first (--legacy-order: smooth before "
              "classify, matching the old chain's order)")
        all_shape = _smooth_tags({k: "shape" for k in variations}, args)
        smooth(edges, nominal, variations, all_shape, args, report, dec)
        print("[3.5] re-validate")
        revalidate(nominal, variations, report)
        print("[1+2] significance tests (on already-smoothed shapes)")
        tags, pvals = classify(nominal, variations, args, report, dec, totals)
    else:
        print("[1+2] significance tests")
        tags, pvals = classify(nominal, variations, args, report, dec, totals)
        if args.no_smoothing:
            print("[3] smoothing skipped")
            report["smoothing"] = "skipped"
        else:
            print("[3] smoothing")
            smooth(edges, nominal, variations, _smooth_tags(tags, args), args, report, dec)
            print("[3.5] re-validate")
            revalidate(nominal, variations, report)

    print("[4] drops")
    apply_drops(edges, nominal, variations, tags, args, report, dec)

    print("[5] symmetrisation")
    policy, syst_rows = symmetrisation(nominal, variations, report, dec,
                                       args.no_symmetrisation, args.flavtag_mirror)

    print("[6] write")
    out = write_tensor(edges, data, nominal, variations, tags, policy, args, report, dec)

    print("[7] tables")
    rows, census = write_tables(dec, syst_rows, args)
    report["decision_census"] = census

    if args.plots:
        outdir = (os.path.join(args.output, args.outname + "_plots")
                  if args.plots == "AUTO" else args.plots)
        print(f"[8] plots -> {outdir}")
        from analysis.templatePlots import plot_all
        nan = float("nan")
        info = {}
        for k in dec:
            row = {c: dec[k].get(c, "") for c in COLUMNS}
            row.update(category=k[0], process=k[1], systematic=k[2])
            row["p_shift"] = row["p_shift"] if row["p_shift"] != "" else nan
            row["p_shape"] = row["p_shape"] if row["p_shape"] != "" else nan
            info[k] = row
        report["plots"] = plot_all(outdir, edges, nominal, raw_all, variations,
                                   info, M.CATEGORIES,
                                   only_kept=(args.plot_what == "kept"),
                                   limit=args.plot_max, jobs=args.plot_jobs,
                                   channels=args.plot_channels)

    report_path = os.path.join(args.output, args.outname + "_report.json")
    with open(report_path, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"\ntensor : {out}\nreport : {report_path}")


if __name__ == "__main__":
    main()
