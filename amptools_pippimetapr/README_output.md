# What the output files mean

Everything goes into the `workdir` from your INI file (default
`EtaPrPi_fit/`), in the folder where you ran the command. Bins are named
`bin_m<mass index>_t<t index>`. For example, `bin_m004_t01` is the 5th
M(η′π⁻) bin and the 2nd −t bin, counting from 0.

## Folder layout

Without wave sets (no `[wavesets]` section):

```
EtaPrPi_fit/
├── bins.json                 bin edges and file names
├── bin_m000_t00/             one folder per bin: split ROOT files, fit, plots
├── bin_m001_t00/
├── ...
├── fits_raw.csv              intermediate: one row per .fit file
├── data_raw.csv              intermediate: one row per data file
├── results.csv               ★ the main table: one row per bin
└── plots/                    ★ every result vs. M(η′π⁻), plus fit_quality.pdf
```

With wave sets, the split ROOT files are made once at the top. Each wave set
then gets its own copy of the layout above, and the sets are compared at the end:

```
EtaPrPi_fit/
├── bins.json
├── bin_m000_t00/             split ROOT files only (shared by all wave sets)
├── ...
├── P1D1/                     one folder per wave set, laid out as above
│   ├── bin_m000_t00/         fit + plots for this set (ROOT files are links)
│   ├── results.csv
│   └── plots/
├── SPD/
├── results_wavesets.csv      ★ all wave sets in one table
└── plots_wavesets/           ★ each result vs. mass, one colour per wave set
```

## Inside one bin folder

| file | what it is |
|---|---|
| `datafile.root`, `accmcfile.root`, `genmcfile.root` (+ `bkgndfile.root` if used) | your input samples cut down to this bin. The names are the INI tokens in lower case. In a wave-set folder they are links to the shared copies. |
| `fit.cfg` | the AmpTools config actually used here: the template with this bin's file names filled in, and with the wave set's waves floating and all others fixed to 0. |
| `fit.log` | full screen output of `fit`. Look here first when a bin fails. |
| `bin_mXXX_tYY.fit` | **the fit result** (AmpTools format): best fit over all random starts. This is the file the next steps read. |
| `bin_mXXX_tYY_0.fit`, `_1.fit`, … | one result per random start (only with `rand`/`-r`). |
| `fit_ranking.txt` | the random starts ranked by likelihood, with Minuit status and error-matrix status. Use it to check that the best fit is not a fluke. |
| `bin_mXXX_tYY.ni` | normalization integrals (acceptance/phase-space integrals of each amplitude pair). AmpTools can reuse these. |
| `fitplots.pdf` | **data vs. fit** for this bin: M(η′π⁻), **M(pπ⁺)**, cosθ and φ in the GJ frame, −t, cosθ and φ in the helicity frame, and Φ (angle to the polarization). Black points are data. The blue histogram is accepted MC weighted by the fitted intensity, which is what the fit predicts you should see. |
| `fitplots.root` | the same histograms in ROOT: `<variable>dat` (data) and `<variable>acc` (fit prediction), e.g. `DeltaPPMassdat` / `DeltaPPMassacc` for M(pπ⁺). With `with_gen` there is also `<variable>gen` (acceptance corrected). With `per_sum` there are also per-sum versions such as `tacc_PositiveRe`. |
| `plotfit.log` | output of the plotter. |

## `bins.json`

The list of bins: name, mass and −t edges, and which sample file sits in each
bin folder. The later steps read it, so leave it alone.

## `results.csv`: one row per bin

The columns in order:

**Bin information**

| column | meaning |
|---|---|
| `bin`, `i_mass`, `i_t` | bin name and indices |
| `mass_lo`, `mass_hi`, `mass` | M(η′π⁻) edges and centre [GeV] |
| `t_lo`, `t_hi` | −t edges [GeV²] (−t to the pπ⁺ system); empty without t bins |
| `fit_ok` | 1 if a .fit file was found and read, 0 if the bin failed |
| `file` | the .fit file used |

**Fit quality**

| column | meaning |
|---|---|
| `likelihood` | −2 ln L as AmpTools reports it (lower is better). Compare it only between fits of the **same bin**, e.g. between wave sets. |
| `minuit_status` | status of the last Minuit command. 0 = OK. |
| `ematrix_status` | error-matrix quality. 3 = full and accurate; less than 3 means the errors can't be trusted. |

**Per wave** (`S0+`, `Pn1+`, `P0+`, `P1+`, `Dn2+`, `Dn1+`, `D0+`, `D1+`, `D2+` and the same with `-`)

Names are `<L><m><reflectivity>`, where `n` means negative m. For example,
`Dn1-` is L = 2, m = −1, negative reflectivity.

| column | meaning |
|---|---|
| `P1+` | intensity of the wave: number of events it contributes. By default this is the **accepted** yield, i.e. what the detector would see. Set `acceptance_corrected = true` in `[gather]` (or pass `-a`) to get the produced (efficiency-corrected) yield instead. |
| `P1+_err` | its uncertainty |
| `P1+_re`, `P1+_im` | the fitted complex production amplitude. The reference wave has `_im` = 0. Waves fixed to 0 show 0. |

**Sums**

| column | meaning |
|---|---|
| `pos_refl`, `neg_refl` | total intensity of all positive / all negative reflectivity waves |
| `all` | total intensity (should be close to the number of data events when not acceptance corrected) |
| `Ssum+`, `Psum+`, `Dsum+`, `Ssum-`, … | coherent sum of all m for one L in one reflectivity |
| `Ssum`, `Psum`, `Dsum` | the same over both reflectivities |
| every `_err` column | the uncertainty of the column before it |

**Phases**

| column | meaning |
|---|---|
| `phase_P1+_D1+` (and `_err`) | phase difference between two waves of the same reflectivity, in **radians**. By default all pairs of floating waves are included; choose them with `convert_args = --phases P1+:D1+`. |

**Data summary** (from `data_info = DATAFILE`)

| column | meaning |
|---|---|
| `data_entries` | number of data events in the bin |
| `data_weighted_events` | sum of the `Weight` branch (equals `data_entries` without weights) |
| `data_mass_mean`, `data_mass_rms`, `data_mass_min`, `data_mass_max` | M(η′π⁻) of the data in this bin. Use `data_info_args = --indices 0,3` to get M(pπ⁺) instead. |

## `fits_raw.csv`, `data_raw.csv`

Intermediate files, written by `convert_to_csv.py` before the bin
information is added. `results.csv` already contains everything in them.

## `plots/`

* `<column>.pdf`: that `results.csv` column vs. M(η′π⁻), with error bars, one
  colour per −t bin. Examples: `P1+.pdf`, `D1+_re.pdf`, `all.pdf`,
  `phase_P1+_D1+.pdf`, `likelihood.pdf`, `data_entries.pdf`.
* `fit_quality.pdf`: the `fitplots.pdf` page of every bin in one file (needs
  `pypdf`, `pdfunite` or `gs`; otherwise the pages stay in the bin folders).

Change the format with `format = png` in `[plot]`, and plot only some
columns with `only = P1+ D1+ likelihood`.

## `results_wavesets.csv` and `plots_wavesets/` (wave sets only)

* `results_wavesets.csv`: every wave set's `results.csv` stacked together,
  with an extra first column `waveset`.
* `plots_wavesets/<column>.pdf`: that column vs. M(η′π⁻) with one colour per
  wave set and one panel per −t bin. The points are shifted slightly
  sideways so they don't overlap. `likelihood.pdf` is the quickest way to see
  which wave set describes the data better in each bin.

## Optional extra files

| file | made when |
|---|---|
| `fits_raw_cov.csv`, `fits_raw_corr.csv` | `convert_args = --covariance` / `--correlation`: the covariance / correlation matrix of the fit parameters, one row per (bin, parameter i, parameter j) |
| `_split/` | temporary folder used during binning. It is removed when binning finishes, so if you see it, binning was interrupted. |
