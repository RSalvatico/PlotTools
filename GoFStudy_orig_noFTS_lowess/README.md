# Pre-unblinding GoF study: orig shapes, lowess smoothing, no flavTag mirror-up

Deliverable: for `Datacards_250826_preUnblinding/orig` (orig shapes, `analysis/prepareTensor.py`'s
lowess/no-flavtag-mirror pipeline defaults) --

- 4 fits: `CR_expected`, `CR_observed`, `SR_expected` (+ tt-vcb scan), `SR_observed`
- prepostfit + impacts (8 POIs x traditional/global) for all 4
- goodness-of-fit toys for the two **observed** configs only (`CR_observed`, `SR_observed`)

See `../README_rabbit.md` for the rabbit backend in general (setup, manual
one-off commands, flag reference). This file is only about reproducing this
one study.

## Setup (once)

```
cmsenv
source ../setup.sh
source ../setup_rabbit.sh
```

## Run locally

```
cd GoFStudy_orig_noFTS_lowess
./run_all.sh
```

One idempotent script, 7 steps (tensors -> fits -> prepostfit -> impacts ->
postfit-conditioned tensors -> GoF toys -> toy plots). Every step checks its
own output and skips if already there -- safe to Ctrl-C and rerun, safe to
rerun after a crash. `JOBS=16 TOYS=5000 ./run_all.sh` overrides local
parallelism/toy count for step 6; each rabbit fit is single-core, so `JOBS`
should not exceed the machine's core count. At `TOYS=5000` the two observed
configs' toy batches are the only real compute, ~2-2.5h each at 16-way local
parallelism -- for anything bigger (this study's actual target is 50k toys
per config), use the HTCondor path below instead of scaling `JOBS` up
locally.

## DPS-scaled cross-check (togglable)

`./run_dps_scaled.sh` repeats only the two **observed** fits (`CR_observed`,
`SR_observed`) -- plus their prepostfit and impacts -- with the tt+bb DPS
processes `ttbb-dps`, `ttbj-dps`, `tt2b-dps` and every one of their systematic
variations scaled up in the simulation. The scale is one knob, `DPS_SCALE`
(default `4.5`), passed straight to `prepareTensor.py --dps-scale`; the output
name carries a matching tag so nominal and scaled results never collide:

```
cd GoFStudy_orig_noFTS_lowess
./run_dps_scaled.sh                 # DPS_SCALE=4.5 -> rabbit/{CR,SR}_observed_dps4p5.hdf5
DPS_SCALE=1 ./run_dps_scaled.sh     # no-op scale; cross-check against run_all.sh
```

`--dps-scale` multiplies the nominal template and every variation leg of those
three processes by the factor and their bin variances by the factor squared, so
the relative MC-stat error and every relative systematic on them are unchanged
-- only their absolute normalisation moves. The related `norm_ttbb-dps` lnN is a
relative kappa and is left as is. Factor `1.0` is a genuine no-op: the tensor is
bit-for-bit what `run_all.sh` builds. `run_all.sh` itself is untouched and stays
the nominal reference; `run_dps_scaled.sh` writes none of its files. Same
blinding discipline as `run_all.sh` (SR unblinds tt-vcb in the fit, impacts
`--redact tt-vcb`, no SR data marker, fitted value never printed).

## Run at scale (50k toys) via HTCondor

`run_all.sh` still does everything through the postfit-conditioned tensors
(steps 1-5); only step 6 (the toy batches) moves to the batch system:

```
cd GoFStudy_orig_noFTS_lowess
./run_all.sh                       # steps 1-5: tensors, fits, prepostfit, impacts (fast, local)

./submit_toys_condor.sh CR_observed rabbit/ourCR_postfitmean.hdf5 observed \
    "--freezeParameters tt-vcb --unblind xsec_.*"
./submit_toys_condor.sh SR_observed rabbit/ourSR_postfitmean.hdf5 observed \
    "--unblind tt-vcb xsec_.*"

condor_q -batch                    # watch the two clusters drain (~1-1.5h each)

./submit_toys_condor.sh --collect CR_observed
./submit_toys_condor.sh --collect SR_observed

python3 ../analysis/plotToyGoF.py --npz rabbit/toyGoF_CR_observed_summary.npz \
    --observed rabbit/CR_observed.hdf5 -o plots/toyGoF_CR_observed.png --label "Lepton + jets"
python3 ../analysis/plotToyGoF.py --npz rabbit/toyGoF_SR_observed_summary.npz \
    --observed rabbit/SR_observed.hdf5 -o plots/toyGoF_SR_observed.png --label "Lepton + jets"
```

`TOYS=50000 JOBS=250 ./submit_toys_condor.sh ...` are the defaults (250 jobs
x 200 toys, ~85 min/job at ~25s/toy, `+MaxRuntime = 21600` = 6h cap -- edit
the script's settings block for a different split). Never pass
`--paramModel` in the `<extra flags>` string: the script already adds it,
doubling it crashes the fitter with `Duplicate parameter names`.

Each Condor job compacts its own toy batch (`nllvalreduced`/`ndfsat` only,
~KB) on the **worker's own local scratch** and transfers back only that tiny
`.npz` + a log -- the raw per-toy fit-result hdf5 (~1.1 MB/toy, hundreds of
MB per job at these batch sizes) never touches AFS. See
`run_gof_toys_lean.sh`'s header for the AFS-quota crash this design avoids,
and its trailing comment for why: `--toysDataMode observed` on a
*postfit-conditioned* tensor (built in step 5, `ourCR_postfitmean.hdf5` /
`ourSR_postfitmean.hdf5`) rather than the raw one -- Poisson-fluctuating
around raw `data_obs` double-counts noise already in the real data (toy mean
came out ~2x the asymptotic expectation when tried on the raw tensor).

## Blinding

- `SR_observed` genuinely unblinds tt-vcb in the fit (real answer, real
  postfit shapes, real GoF) -- the discipline is on reporting, not
  computation. **Never print or state its fitted value or uncertainty.**
- Its impacts plots pass `--redact tt-vcb` to `rabbitPlotImpacts.py` for
  *every* POI, not only when tt-vcb is the one being plotted: tt-vcb is a
  free parameter correlated with the tt+X norms, so its real pull can leak
  as printed text into *another* POI's impact plot too (this happened once
  during development -- see git log for the fix and how it was caught).
- Its prepostfit plot never draws an SR data marker by default, regardless
  of the fit's own blinding state (`rabbitPlotStack.py --show-signal-data`
  is an explicit, separate opt-in -- `run_all.sh` makes a second
  `SR_observed_withSRdata` copy with it, alongside the always-blinded
  default).
- `CR_expected` / `CR_observed` freeze tt-vcb (`--freezeParameters tt-vcb`);
  it is not a real measurement there, nothing to redact.
- `SR_expected` is Asimov (`-t -1`); no real data involved anywhere in it.

## Output layout

```
rabbit/                          fit results, tensors, toy batches (gitignored)
  ourSR.hdf5, ourCR.hdf5           tensors (lowess + no flavtag mirror-up, prepareTensor.py defaults)
  ourSR_postfitmean.hdf5,          copies with hdata_obs patched to the observed
  ourCR_postfitmean.hdf5           config's own postfit expectation (GoF toy conditioning)
  CR_expected/observed.hdf5,       the 4 fit results
  SR_expected/observed.hdf5
  toyGoF_<name>_batch*.hdf5        raw per-batch toy fits (local run only; deleted right
                                    after compaction -- never accumulates on disk)
  toyGoF_<name>_batch*_summary.npz per-batch compacted (q, ndf) arrays, ~KB each
  toyGoF_<name>_summary.npz        combined summary, consumed by plotToyGoF.py
  condor_<name>/                   HTCondor submit file + executable + job logs (submit_toys_condor.sh)
plots/
  prepostfit/                      5 = 4 configs + SR_observed_withSRdata, pre+postfit, linear+log, png+pdf
  impacts/                         4 configs x 8 POIs x traditional/global, png+pdf (128 files)
  toyGoF_CR_observed.png/.pdf      toy q-histogram + observed line + chi2 fit + p-value
  toyGoF_SR_observed.png/.pdf
```
