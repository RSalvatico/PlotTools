#!/usr/bin/env python3
"""Likelihood scan plot for every scanned parameter.

    analysis/rabbitPlotScan.py fit.hdf5 -o plots/scans

One figure per `nll_scan_<param>` present in the fit result. The header carries
the profile-likelihood interval read off the 2*deltaNLL = 1 crossings and, for
comparison, the parabolic uncertainty from the inverse Hessian.

Two conventions have to be undone:
  * rabbit stores the scan as deltaNLL, so it is doubled here to make the 1 sigma
    crossing sit at 1;
  * a POI is scanned in the internal variable x = sqrt(mu), so the axis is
    squared back to mu (which is what makes the interval asymmetric).
"""

import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rabbit import io_tools

from analysis.rabbitPlotStyle import cms_label


def crossing(x, y, best, level):
    """(low, high) where 2*deltaNLL crosses `level`, nan if outside the range."""
    lo_m, hi_m = x < best, x > best
    lo = np.interp(level, y[lo_m][::-1], x[lo_m][::-1],
                   left=np.nan, right=np.nan) if lo_m.any() else np.nan
    hi = np.interp(level, y[hi_m], x[hi_m],
                   left=np.nan, right=np.nan) if hi_m.any() else np.nan
    if lo_m.any() and y[lo_m].max() < level:
        lo = np.nan
    if hi_m.any() and y[hi_m].max() < level:
        hi = np.nan
    return lo, hi


def read_scan(fr, name, pois):
    """(x, 2*deltaNLL) for `name` from a rabbit fit result, axis squared back
    to mu for a POI (same conventions as the main curve)."""
    h = fr[f"nll_scan_{name}"].get()
    x = np.array([float(v) for v in np.array(h.axes[0])])
    y = 2.0 * h.values()                            # rabbit stores deltaNLL
    if name in pois:                                # internal x = sqrt(mu)
        x = x ** 2
    return x, y


def read_combine(path, param):
    """(x, 2*deltaNLL) from a combine MultiDimFit --algo grid output."""
    import uproot
    t = uproot.open(path)["limit"]
    if param not in t.keys():
        return None
    x = np.asarray(t[param].array(library="np"), dtype=float)
    d = np.asarray(t["deltaNLL"].array(library="np"), dtype=float)
    ok = np.isfinite(x) & np.isfinite(d)
    x, d = x[ok], d[ok]
    o = np.argsort(x)
    return x[o], 2.0 * d[o]


def draw(path, name, x, y, hesse, title, combine=None, statonly=None, asimov=False):
    best = x[np.argmin(y)]
    lo1, hi1 = crossing(x, y, best, 1.0)
    lo2, hi2 = crossing(x, y, best, 4.0)

    fig, a = plt.subplots(figsize=(10, 9), dpi=160)
    a.grid(True, lw=0.5, alpha=0.35, zorder=0)
    lbl = "Rabbit (stat+syst)" if statonly is not None else "Rabbit"
    a.plot(x, y, color="#1f4fd8", lw=2.2, label=lbl, zorder=3)
    if statonly is not None:
        sx, sy = statonly
        a.plot(sx, sy, color="#e08a00", lw=2.0, ls="--", label="Rabbit (stat only)",
               zorder=3)
        sbest = sx[np.argmin(sy)]
        slo1, shi1 = crossing(sx, sy, sbest, 1.0)
    if combine is not None:
        cx, cy = combine
        a.plot(cx, cy, color="#cc2222", lw=1.8, ls="--", label="Combine", zorder=3)
        cbest = cx[np.argmin(cy)]
        clo, chi = crossing(cx, cy, cbest, 1.0)
    a.axvline(best, color="k", lw=1.1, ls="--", zorder=2)
    for lvl, lab, c in ((1.0, r"1$\,\sigma$", "#888888"), (4.0, r"2$\,\sigma$", "#bbbbbb")):
        a.axhline(lvl, color=c, lw=1.0, zorder=1)
        a.text(x.max(), lvl, f" {lab}", va="center", ha="left", fontsize=13, color=c)
    for v in (lo1, hi1):
        if np.isfinite(v):
            a.plot([v, v], [0, 1.0], color="#888888", lw=1.0, zorder=1)

    a.set_xlabel(name, fontsize=18)
    a.set_ylabel(r"$2\,\Delta\mathrm{NLL}$", fontsize=18)
    a.set_ylim(0, min(9.0, float(np.nanmax(y)) * 1.05))
    a.set_xlim(x.min(), x.max())

    up = hi1 - best if np.isfinite(hi1) else np.nan
    dn = best - lo1 if np.isfinite(lo1) else np.nan
    head = f"{name} = {best:.4f}  $-{dn:.4f}/+{up:.4f}$" if np.isfinite(up + dn) \
        else f"{name} = {best:.4f}  (1 sigma outside the scan range)"
    a.text(0.04, 0.66, head, transform=a.transAxes, ha="left", va="top", fontsize=14,
           fontweight="bold")
    sub = f"Hessian $\\pm${hesse:.4f}" if np.isfinite(hesse) else ""
    if np.isfinite(lo2) and np.isfinite(hi2):
        sub += f"     2$\\sigma$ [{lo2:.4f}, {hi2:.4f}]"
    if sub:
        a.text(0.04, 0.60, sub, transform=a.transAxes, ha="left", va="top",
               fontsize=11, color="#555555")
    if statonly is not None:
        su = shi1 - sbest if np.isfinite(shi1) else np.nan
        sd = sbest - slo1 if np.isfinite(slo1) else np.nan
        stxt = f"stat. only  $-{sd:.4f}/+{su:.4f}$" if np.isfinite(su + sd) \
            else "stat. only  (1 sigma outside the scan range)"
        a.text(0.04, 0.54, stxt, transform=a.transAxes, ha="left", va="top",
               fontsize=11, color="#e08a00")
    if title:
        a.text(0.04, 0.48 if statonly is not None else 0.54, title,
               transform=a.transAxes, ha="left", va="top", fontsize=11, color="#555555")
    if combine is not None:
        cu = chi - cbest if np.isfinite(chi) else np.nan
        cd = cbest - clo if np.isfinite(clo) else np.nan
        a.text(0.5, -0.14, f"Combine: {cbest:.4f}  $-{cd:.4f}/+{cu:.4f}$",
               transform=a.transAxes, ha="center", fontsize=11, color="#cc2222")
    a.legend(fontsize=12, frameon=False, loc="upper right")
    cms_label(a, data=not asimov)

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
    p.add_argument("--result", default=None)
    p.add_argument("--params", nargs="*", default=None,
                   help="default: everything that was scanned")
    p.add_argument("--combine", default=None,
                   help="combine MultiDimFit --algo grid output to overlay. Use "
                        "{param} in the path for a per-parameter file, e.g. "
                        "'scans/higgsCombine{param}.MultiDimFit.mH120.root'")
    p.add_argument("--statonly", default=None,
                   help="second rabbit scan fit result (systematics frozen) to "
                        "overlay as a 'stat only' curve")
    p.add_argument("--asimov", action="store_true",
                   help="fit ran on Asimov/expected data: CMS label reads 'Simulation'")
    args = p.parse_args()

    fr, meta = io_tools.get_fitresult(args.fitresult, result=args.result, meta=True)
    pois = {str(v) for v in np.asarray(meta["pois"]).astype(str)}
    tag = args.postfix or os.path.splitext(os.path.basename(args.fitresult))[0]

    scanned = [k[len("nll_scan_"):] for k in fr.keys() if k.startswith("nll_scan_")]
    names = [n for n in scanned if args.params is None or n in args.params]
    if not names:
        print(f"  no scans in {args.fitresult}")
        return

    fr_stat = None
    if args.statonly:
        fr_stat = io_tools.get_fitresult(args.statonly, result=args.result)

    parms = fr["parms"].get()
    for name in sorted(names):
        x, y = read_scan(fr, name, pois)
        sx = float(np.sqrt(parms[{"parms": name}].variance))
        if name in pois:                           # internal x = sqrt(mu)
            sx = 2.0 * abs(float(parms[{"parms": name}].value)) * sx
        stat = None
        if fr_stat is not None:
            if f"nll_scan_{name}" in fr_stat.keys():
                stat = read_scan(fr_stat, name, pois)
            else:
                print(f"  (no stat-only scan for {name} in {args.statonly})")
        comb = None
        if args.combine:
            cpath = args.combine.replace("{param}", name)
            try:
                comb = read_combine(cpath, name)
            except Exception as exc:
                print(f"  (no combine overlay for {name}: {exc})")
        draw(os.path.join(args.outdir, f"scan_{name}_{tag}.png"), name, x, y, sx,
             tag, combine=comb, statonly=stat, asimov=args.asimov)


if __name__ == "__main__":
    main()
