# Rabbit backend

TensorFlow-based alternative to the classic Combine workflow described in
`README.md`, sharing the same input: `<card>_shapes.root` from
`prepareDatacards.py`. The model (categories, processes, systematics metadata)
is defined once in `configs/model.py`.

## Setup

```
cmsenv
git submodule update --init rabbit   # first time only, fetches the fitter
./setup_rabbit_env.sh                # builds rabbit_env/, first time only
source setup_rabbit.sh               # activate; `rabbit_off` to leave
```

## Run

```
./runRabbitFits.sh
```
builds tensors from a shapes file (`analysis/prepareTensor.py`), fits them
(`rabbit_fit.py`, from the `rabbit/` submodule), and plots pre/post-fit
stacks, likelihood scans, and impacts (`analysis/rabbitPlot{Stack,Scan,Impacts}.py`,
all CMS-styled, grouped legends, png+pdf). Edit the settings block at the top
of the script to point at a different shapes file.

Cross-check against classic Combine on one datacard:
```
./validateCombineRabbit.sh <carddir> <outdir>
```

## Pre-unblinding study: CR + SR, expected/observed, impacts, prepostfit, GoF

`GoFStudy_orig_noFTS_lowess/README.md` -- setup, local run, HTCondor run for
large toy counts, blinding discipline, output layout.

## Defaults (lowess + flavTag mirror-up off)

`analysis/prepareTensor.py` defaults to `--smoothing-method lowess` and keeps
the stored flavTag down leg as-is (no `--flavtag-mirror`). To get the OLD
defaults (spline smoothing, rebuilt flavTag down leg), pass
`--smoothing-method spline --flavtag-mirror` explicitly.

`--dps-scale FACTOR` (default `1.0`, a no-op) multiplies the yield of the tt+bb
DPS processes `ttbb-dps, ttbj-dps, tt2b-dps` -- nominal and every systematic
variation leg -- by `FACTOR`, emulating a change to their simulated cross
section; bin variances scale as `FACTOR**2`, so every relative uncertainty on
those processes is preserved. See
`GoFStudy_orig_noFTS_lowess/run_dps_scaled.sh` for the CR+SR observed fits at
`DPS_SCALE=4.5`.

## Manual: one fit, then impacts / pre-postfit / scan / GoF toys

Build a tensor and fit it (SR unmasked here; add `--mask '_SR$'` for a
CR-only tensor):
```
python3 analysis/prepareTensor.py <shapes.root> -o out --outname SR
PM="--paramModel Mu --paramModel analysis.rabbit_models.FreeNorm ttbb,ttbj,tt2b,ttcc,ttcj,tt2c,ttLF"
rabbit_fit.py out/SR.hdf5 -o out --outname SR_fit -t -1 ${PM} \
    --doImpacts --globalImpacts --saveHists --saveHistsPerProcess --computeHistErrors \
    --scan tt-vcb --scanPoints 31
```
`-t -1` Asimov, `-t 0` real data, `-t N` (N>=1) generates N toys in one job.

**Impacts** (one call per POI, traditional or global):
```
python3 analysis/rabbitPlotImpacts.py out/SR_fit.hdf5 --poi tt-vcb \
    -o out/plots --postfix SR --impact-type traditional
```
Add `--redact tt-vcb` (or any parameter name) to hide that parameter's own
central value/uncertainty everywhere it appears as a pull row in this plot --
including when it merely correlates with `--poi` rather than being `--poi`
itself -- while every other parameter's impact still shows in full. Pass it
for every `--poi` call whenever a blinded parameter might leak by
correlation, not only the calls where it is the plotted POI. If it equals
`--poi`, the "total" subtitle is hidden too.

**Pre/postfit stack** (SR data marker is blinded by default):
```
python3 analysis/rabbitPlotStack.py out/SR_fit.hdf5 -o out/plots --postfix SR --logy
```
Add `--show-signal-data` to also draw the SR data marker -- only for a fit
whose SR is genuinely meant to be shown; it does not check blinding for you.

**Scan**: already produced above via `--scan`; read/plot with
`analysis/rabbitPlotScan.py out/SR_fit.hdf5 -o out/plots --params tt-vcb`.

**GoF toys**, storage-safe (keeps ~KB/config, not ~GB -- see
`GoFStudy_orig_noFTS_lowess/run_gof_toys_lean.sh` for the reasoning):
```
cd GoFStudy_orig_noFTS_lowess
./run_gof_toys_lean.sh <name> <tensor.hdf5> <expected|observed> "<extra rabbit_fit.py flags>"
python3 ../analysis/plotToyGoF.py --npz rabbit/toyGoF_<name>_summary.npz \
    --observed rabbit/<name>.hdf5 -o plots/toyGoF_<name>.png
```
`expected` mode is the correct choice for an Asimov config (the toy mean IS
the null hypothesis already). `observed` mode needs a postfit-conditioned
tensor first (Poisson-fluctuating around raw `data_obs` double-counts noise,
verified: inflates the toy mean ~2x) -- see `run_all.sh`'s
`patch_postfit_mean` step for the recipe.
GoF itself (`nllvalreduced`/`ndfsat`) needs no extra flags, always computed.

## Layout

- `analysis/` -- tensor building, fit glue, and all plotting
  (`rabbitPlot{Stack,Scan,Impacts,Style}.py`, `prepareTensor.py`,
  `smoothing.py`, `validation*.py`, ...).
- `configs/model.py` -- single source of truth for categories, processes, and
  per-systematic statistical metadata (correlation class, symmetrisation).
- `rabbit/` -- the vendored fitter itself: a git submodule with its own
  history and README, not part of this repo's commits.
- `COMBINE_FIXES.md` -- Combine-side bugs found while building this backend.
- `RABBIT_TENSOR_GUIDE.md` -- tensor pipeline internals, decision policy.
