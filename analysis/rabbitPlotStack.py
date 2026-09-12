#!/usr/bin/env python3
"""Stacked pre/post-fit plot with every category unrolled into one histogram.

rabbit_plot_hists.py draws one figure per channel. This draws the whole fit at
once: all categories concatenated along x, every process stacked, data on top,
category boundaries marked.

    analysis/rabbitPlotStack.py fit.hdf5 -o plots/
    analysis/rabbitPlotStack.py fit.hdf5 -o plots/ --logy

Masked channels have no data, so only the prediction is drawn there.
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

from rabbit import io_tools

from configs import model as M
from analysis.rabbitPlotStyle import cms_label, group_stack, region_label

# suggested factor for --dps-scale; scaling is off (1.0) unless requested.
DPS_SCALE_DEFAULT = 4.52

def collect(result, fittype, channels, show_signal_data=False, dps_scale=1.0):
    """Concatenate every channel: per-process yields, total, its error, data.

    The SR is blinded by default: an Asimov fit run with the SR unmasked
    still writes an (expected) hist_data_obs there, but the plot must not
    show a data marker in the blinded signal region unless the caller
    explicitly opts in via show_signal_data=True (a genuinely unblinded fit).
    """
    stack, total, err, data, edges, bounds = {}, [], [], [], 0, []
    for ch in channels:
        c = result[ch]
        inc = c[f"hist_{fittype}_inclusive"].get()
        n = inc.values().flatten()
        h = c[f"hist_{fittype}"].get()
        # the inclusive total doesn't know about the DPS scaling applied to the
        # stack below, so shift it by the same amount to keep the ratio panel
        # (data / total) consistent with what's actually drawn in the stack.
        dps_delta = np.zeros_like(n)
        for p in [str(x) for x in np.array(h.axes["processes"])]:
            vals = h[{"processes": p}].values().flatten()
            if "dps" in p.lower() and dps_scale != 1.0:
                dps_delta = dps_delta + (dps_scale - 1) * vals
                vals = vals * dps_scale
            stack.setdefault(p, []).append(vals)
        total.append(n + dps_delta)
        err.append(np.sqrt(np.abs(inc.variances().flatten())))
        d = c["hist_data_obs"].get() if "hist_data_obs" in c.keys() else None
        vals = d.values().flatten() if d is not None else np.full(len(n), np.nan)
        if ch in M.SIGNAL_REGIONS and not show_signal_data:
            vals = np.full(len(n), np.nan)
        data.append(vals)
        edges += len(n)
        bounds.append(edges)
    return ({p: np.concatenate(v) for p, v in stack.items()},
            np.concatenate(total), np.concatenate(err), np.concatenate(data), bounds)


def draw(path, stack, total, err, data, bounds, channels, fittype, logy,
         asimov, prefit_total=None, prefit_err=None):
    labels, colors, arrays = group_stack(stack)
    n = len(total)
    x = np.arange(n) + 0.5
    edges = np.arange(n + 1)

    fig, (a, r) = plt.subplots(
        2, 1, figsize=(max(12, n * 0.24), 8), dpi=160, sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0.05})

    a.stackplot(edges, *[np.append(v, v[-1]) for v in arrays], step="post",
                labels=labels, colors=colors, edgecolor="white", linewidth=0.2)
    a.fill_between(edges, np.append(total - err, (total - err)[-1]),
                   np.append(total + err, (total + err)[-1]), step="post",
                   facecolor="none", edgecolor="k", hatch="////", linewidth=0,
                   alpha=0.6, label="Total uncertainty")
    ok = np.isfinite(data)
    dlabel = "Asimov" if asimov else "Data"
    a.errorbar(x[ok], data[ok], yerr=np.sqrt(np.maximum(data[ok], 0)), fmt="ko",
               ms=4, lw=1.2, label=dlabel)

    safe = np.where(total > 0, total, np.nan)
    r.axhline(1.0, color="k", lw=1)
    r.fill_between(edges, np.append(1 - err / safe, (1 - err / safe)[-1]),
                   np.append(1 + err / safe, (1 + err / safe)[-1]), step="post",
                   facecolor="none", edgecolor="k", hatch="////", linewidth=0,
                   alpha=0.6)
    r.errorbar(x[ok], (data / safe)[ok], yerr=(np.sqrt(np.maximum(data, 0)) / safe)[ok],
               fmt="ko", ms=4, lw=1.2, zorder=3)
    if prefit_total is not None:
        # red line: prefit data/MC ratio (MC = prefit sum of processes)
        # red band: prefit MC relative stat. uncertainty, centered on 1
        pre_safe = np.where(prefit_total > 0, prefit_total, np.nan)
        pre_ratio = data / pre_safe
        pre_band = prefit_err / pre_safe
        r.fill_between(edges, np.append(1 - pre_band, (1 - pre_band)[-1]),
                       np.append(1 + pre_band, (1 + pre_band)[-1]),
                       step="post", facecolor="#cc2222", edgecolor="none",
                       alpha=0.18, zorder=1)
        r.step(edges, np.append(pre_ratio, pre_ratio[-1]), where="post",
               color="#cc2222", lw=1.3, zorder=2, label="Pre-fit data / pred.")

    for b in bounds[:-1]:
        for ax in (a, r):
            ax.axvline(b, color="k", lw=0.8, ls="--", alpha=0.5)
    lo = 0
    for b, ch in zip(bounds, channels):
        a.text((lo + b) / 2, 0.74, region_label(ch), transform=
               a.get_xaxis_transform(), ha="center", va="top", fontsize=9,
               rotation=90, color="#333333")
        lo = b

    if logy:
        a.set_yscale("log")
        a.set_ylim(max(1e-2, np.nanmin(total[total > 0]) * 0.2), np.nanmax(total) * 60)
    else:
        a.set_ylim(0, np.nanmax(np.where(ok, data, total)) * 1.6)
    # zoom the ratio panel to the actual spread: postfit bands/residuals are
    # much smaller than prefit, a fixed +-40% window wastes most of the panel.
    resid = np.abs((data / safe)[ok] - 1.0) if ok.any() else np.array([])
    band = (err / safe)[np.isfinite(err / safe)]
    extra = []
    if prefit_total is not None:
        pre_safe = np.where(prefit_total > 0, prefit_total, np.nan)
        pre_resid = np.abs((data / pre_safe)[ok] - 1.0) if ok.any() else np.array([])
        pre_rel = prefit_err / pre_safe
        extra.extend([pre_resid[np.isfinite(pre_resid)],
                      pre_rel[np.isfinite(pre_rel)]])
    spread = max([0.0] + [np.nanmax(a_) for a_ in (resid, band, *extra) if a_.size])
    half = max(spread * 1.35, 0.05)
    r.set_ylim(1 - half, 1 + half)
    r.set_ylabel("Data / pred.", fontsize=14)
    r.set_xlabel("Unrolled bin (category boundaries dashed)", fontsize=14)
    a.set_ylabel("Events", fontsize=14)
    a.set_xlim(0, n)

    handles, hlabels = a.get_legend_handles_labels()
    # stackplot fills bottom-up; show the legend top-down (signal first).
    n_groups = len(labels)
    handles = handles[:n_groups][::-1] + handles[n_groups:]
    hlabels = hlabels[:n_groups][::-1] + hlabels[n_groups:]
    if prefit_total is not None:
        handles.append(Line2D([0], [0], color="#cc2222", lw=1.5))
        hlabels.append("Pre-fit data / pred.")
    a.legend(handles, hlabels, fontsize=10, ncol=5, frameon=False, loc="lower left",
             bbox_to_anchor=(0, 1.10, 1, 0.2), mode="expand", borderaxespad=0)

    cms_label(a, data=not asimov, loc=0)
    a.text(0.98, 0.96, "Post-fit" if fittype == "postfit" else "Pre-fit",
           transform=a.transAxes, ha="right", va="top", fontsize=13,
           fontweight="bold", color="#222222")

    chi2 = float(np.nansum((data[ok] - total[ok]) ** 2 / np.maximum(total[ok], 1e-9)))
    if ok.any():
        r.text(0.995, 0.06, f"$\\chi^2$/bin = {chi2 / ok.sum():.2f}",
               transform=r.transAxes, ha="right", fontsize=9, color="#444444")

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.splitext(path)[0] + "." + ext, bbox_inches="tight")
    plt.close(fig)
    print(f"  {os.path.basename(path)}")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fitresult")
    p.add_argument("-o", "--outdir", default="./")
    p.add_argument("--postfix", default=None)
    p.add_argument("--result", default=None, help="e.g. asimov")
    p.add_argument("--logy", action="store_true")
    p.add_argument("-m", "--mapping", default="BaseMapping")
    p.add_argument("--asimov", action="store_true",
                   help="fit ran on Asimov/expected data (-t -1 or -t >= 1): "
                        "CMS label reads 'Simulation', points labelled 'Asimov'")
    p.add_argument("--show-signal-data", action="store_true",
                   help="draw the data marker in signal-region channels too. "
                        "Off by default: the SR is blinded here regardless of the "
                        "fit's own blinding state, so this is an explicit opt-in for "
                        "a fit whose SR data is meant to be shown.")
    p.add_argument("--dps-scale", type=float, default=1.0, nargs="?",
                   const=DPS_SCALE_DEFAULT,
                   help="scale processes with 'dps' in their name (e.g. tt2b-dps) "
                        f"by this factor. Off by default (1.0); bare flag uses "
                        f"{DPS_SCALE_DEFAULT}.")
    args = p.parse_args()

    fitresult = io_tools.get_fitresult(args.fitresult, result=args.result)
    result = fitresult["mappings"][args.mapping]["channels"]
    channels = [c for c in M.CATEGORIES if c in result.keys()]
    tag = args.postfix or os.path.splitext(os.path.basename(args.fitresult))[0]

    prefit_cache = None
    for fittype in ("prefit", "postfit"):
        if f"hist_{fittype}_inclusive" not in result[channels[0]].keys():
            print(f"  no {fittype} histograms in {args.fitresult}")
            continue
        stack, total, err, data, bounds = collect(result, fittype, channels,
                                                   show_signal_data=args.show_signal_data,
                                                   dps_scale=args.dps_scale)
        pre_kwargs = {}
        if fittype == "postfit" and prefit_cache is not None:
            pre_kwargs = dict(prefit_total=prefit_cache[0], prefit_err=prefit_cache[1])
        elif fittype == "prefit":
            prefit_cache = (total, err)
        for logy in ({False, True} if args.logy else {False}):
            name = f"{tag}_{fittype}" + ("_log" if logy else "") + ".png"
            draw(os.path.join(args.outdir, name), stack, total, err, data,
                 bounds, channels, fittype, logy, args.asimov, **pre_kwargs)


if __name__ == "__main__":
    main()

