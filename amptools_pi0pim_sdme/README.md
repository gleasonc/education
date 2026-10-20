# AmpTools SDME fits for γd → ρ⁻p(p), ρ⁻ → π⁻π⁰

These scripts fit the ρ spin-density matrix elements (SDMEs) in bins of −t
with halld_sim's `TwoPiAngles`, for the four beam-polarization orientations
(0°, 45°, 90°, 135°). One command splits the trees, runs the fits, collects
the results and makes the plots, including data-vs-fit plots for every fit.

It is the SDME version of `../amptools_pippimetapr` (the η′π⁻ wave fits).
That folder is unchanged.

## 1. Input files

AmpTools trees from the v2 DSelector (`DSelector_gd_pi0pimprotmissprot__B4_F4_T2_S2_v2`):
tree `kin`, final state **proton, π⁻, π⁰** (TwoPiAngles: recoil, p1, p2), a
`Weight` branch (accidental weights, read by `ROOTDataReader`) and an `Int_t`
branch `pol` with the orientation (0/45/90/135, −1 = AMO).

You do **not** need to run `SplitByPol.C` first: the `bin` step splits the data
by `pol` itself.

## 2. Running it

```tcsh
cp example_sdme.ini my_sdme.ini      # set DATAFILE, ACCMCFILE, GENMCFILE and the binning
./run_all.csh my_sdme.ini            # or: python3 amptools_sdme.py all -c my_sdme.ini
```

Needs `root` and the halld_sim `fit` program in your path (the GlueX
environment), and the Python packages `uproot` and `matplotlib`
(`run_all.csh` installs them if missing). Nothing has to be compiled.

The steps can also run one at a time:

```tcsh
python3 amptools_sdme.py bin    -c my_sdme.ini   # split the files
python3 amptools_sdme.py fit    -c my_sdme.ini -r 5
python3 amptools_sdme.py gather -c my_sdme.ini
python3 amptools_sdme.py plot   -c my_sdme.ini
python3 amptools_sdme.py plotfit -c my_sdme.ini  # data vs. fit
```

Options on the command line override the INI. `--dry-run` prints the commands
without running them. `fit` and `plotfit` take `--bin NAME ...` to run only
some bins; `plotfit` skips fit plots that are newer than their `.fit` file
unless `--redo-plots` is given.

### On the JLab farm (swif2)

`submit_sdme_jobs.py` submits one swif2 job per bin; each job runs the five
fits of its bin (combined + one per orientation) and their fit plots. The job
settings (workflow, account, partition, `os = el9`, cores, ram, disk, time and
the tcsh lines that set up the GlueX environment) are in the `[farm]` section
of the INI.

```tcsh
python3 amptools_sdme.py bin -c my_sdme.ini            # on ifarm, once
python3 submit_sdme_jobs.py -c my_sdme.ini --dry-run   # check the job scripts
python3 submit_sdme_jobs.py -c my_sdme.ini             # submit (add --fit-args "-r 5")
swif2 status -workflow sdme_pi0pim                     # wait until all jobs are done
python3 amptools_sdme.py gather  -c my_sdme.ini
python3 amptools_sdme.py plot    -c my_sdme.ini
python3 amptools_sdme.py plotfit -c my_sdme.ini        # only merges the per-bin pages
```

Job scripts go to `<workdir>/jobs/<bin>.csh`, logs to `<workdir>/jobs/log/`.
Failed jobs can be resubmitted with `swif2 retry-jobs` or by running
`submit_sdme_jobs.py --bin NAME`.

## 3. What each step does

**bin.** `split_bins.C` splits every sample into M(π⁻π⁰) × −t bins
(−t = −(beam − π⁻π⁰)², the default is one mass bin 0.6–0.95 GeV and ten −t bins
from 0.1 to 1.1 GeV²). The samples listed under `split` (default `DATAFILE`)
are also split by orientation using the `pol` branch, giving
`datafile_00.root … datafile_135.root` in each bin. Samples that are not split
(the flat MC) are used whole for every orientation. If your MC also has the
right `pol` values, add it to `split`.

**fit.** Each bin gets two kinds of fit, each in its own folder:

| folder | template | what it is |
|---|---|---|
| `all/` | `fit_pi0pim_SDME_4pol_TEMPLATE.cfg` | one fit of all four orientations with shared SDMEs, a polarization fraction and angle per orientation, and a free scale per orientation (`parScale00` fixed to 1) |
| `pol00/` … `pol135/` | `fit_pi0pim_SDME_1pol_TEMPLATE.cfg` | one fit per orientation, for the per-orientation comparison |

Leave `cfg` or `cfg_single` empty in the INI to skip one kind. The polarization
fractions and angles are the `define polVal_..` / `polAngle_..` lines at the
top of both templates. With `fit -r N`, the `parSDME` keyword lines make
halld_sim's `fit` also start each try from random SDMEs that form a valid
density matrix; the best try is kept as `<name>.fit`. Bins with no events in a
needed file are skipped.

**gather.** Reads every `.fit` file and writes `results.csv`, one row per bin
and fit (`fit` = `all`, `pol00`, …), with the nine SDMEs, the derived
quantities below and the scale parameters, each with `_err`. Errors of the
derived quantities come from the fit's full covariance matrix.

**plot.** Writes three PDFs to `plots/`, one page per mass bin, one colour per
fit (pol=0 black, 45 green, 90 red, 135 blue, combined fit purple):

* `sdmes.pdf`: the nine SDMEs vs −t in a 3 × 3 grid
* `parity_components.pdf`: natural (top) and unnatural (bottom) components,
  laid out as Fig. 7 of [arXiv:2305.09047](https://arxiv.org/abs/2305.09047)
* `asymmetries.pdf`: the parity asymmetry P_σ and the beam asymmetry Σ

Each also comes as `sdmes_total.pdf`, `parity_components_total.pdf` and
`asymmetries_total.pdf` with only the combined fit (the total result, no
per-orientation points).

Dashed lines mark s-channel helicity conservation with natural-parity exchange
(ρᴺ₁₁ = 0.5, P_σ = Σ = 1). `--fits all` plots only the combined fit;
`--format png` also writes PNG files.

**plotfit.** Checks each fit by comparing the data with the accepted MC
weighted by the fitted intensity, as in Fig. 5 of the paper, in cos ϑ_H, φ_H,
Φ, Φ − φ_H, M(π⁻π⁰) and −t. Every fit folder gets `fitplots.pdf` (for the
combined fit: a page for the sum of the orientations, then one per
orientation), and `plots/fit_quality_all.pdf`, `plots/fit_quality_pol00.pdf`, …
collect all bins (needs `pypdf`, `pdfunite` or `gs` to merge). Data events use
their `Weight`; the fit histogram is normalized to the number of events the fit
predicts, so it is an absolute comparison. `sdme_fit_plotter.py` also runs on
its own inside a fit folder:

```tcsh
cd SDME_fit/bin_m000_t03/all
python3 ../../../sdme_fit_plotter.py all.fit --check
```

`--check` compares the predicted yield with the fit's normalization integrals
and should print `OK`. Leave the step out of `all` with `--no-plotfit` or
`plotfit = false`.

## 4. Derived quantities

From arXiv:2305.09047 Eq. (22), ρᴺ'ᵁᵢₖ = ½[ρ⁰ᵢₖ ∓ (−1)ⁱ ρ¹₋ᵢₖ], written with
the TwoPiAngles parameters:

| quantity | formula | parameters |
|---|---|---|
| ρᴺ₀₀, ρᵁ₀₀ | ½(ρ⁰₀₀ ∓ ρ¹₀₀) | `rho000`, `rho001` |
| ρᴺ₁₁, ρᵁ₁₁ | ½(ρ⁰₁₁ ± ρ¹₁₋₁), ρ⁰₁₁ = (1 − ρ⁰₀₀)/2 | `rho000`, `rho1m11` |
| ρᴺ₁₋₁, ρᵁ₁₋₁ | ½(ρ⁰₁₋₁ ± ρ¹₁₁) | `rho1m10`, `rho111` |
| Re ρᴺ₁₀, Re ρᵁ₁₀ | ½(Re ρ⁰₁₀ ∓ Re ρ¹₁₀) | `rho100`, `rho101` |
| P_σ (Eq. 23) | 2ρ¹₁₋₁ − ρ¹₀₀ | `rho1m11`, `rho001` |
| Σ | 2ρ¹₁₁ + ρ¹₀₀ | `rho111`, `rho001` |

Σ is the beam asymmetry of the decay-angle-integrated cross section,
dσ/dΦ ∝ 1 − P_γ Σ cos 2Φ, in the same sign convention as TwoPiAngles
(Schilling et al.).

## 5. Files

| file | what it is |
|---|---|
| `run_all.csh` / `run_all.sh` | one command for everything (tcsh / bash) |
| `amptools_sdme.py` | the driver: bin → fit → gather → plot → plotfit |
| `example_sdme.ini` | example settings; copy and edit |
| `fit_pi0pim_SDME_4pol_TEMPLATE.cfg` | combined fit of the four orientations |
| `fit_pi0pim_SDME_1pol_TEMPLATE.cfg` | fit of one orientation |
| `split_bins.C` | ROOT macro: splits trees into M × −t bins, optionally one orientation |
| `sdme_fit_plotter.py` | data vs. fit plots for one .fit file |
| `submit_sdme_jobs.py` | submits one swif2 farm job per bin |

## 6. Test

The chain was run end to end with AmpTools `fit` built from the AmpTools and
halld_sim sources (TwoPiAngles, ROOTDataReader) on a toy sample: 80k
γp → pπ⁻π⁰ events generated with known, −t-dependent SDMEs and random
orientations, with 200k phase-space events as MC. The combined fits returned
the input SDMEs in every bin (pulls with mean 0.0 and RMS 1.1). The
single-orientation fits agreed as well, except one fit at low −t whose
SDMEs make the angular distribution W negative in part of the angular range
(TwoPiAngles uses |W|), so it is not a physical solution; use `-r` for
several starts and compare with the combined fit.
