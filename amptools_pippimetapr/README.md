# AmpTools fit automation for γp → η′π⁻π⁺p

These scripts run a binned (mass-independent) AmpTools fit of η′π⁻ waves
recoiling against a pπ⁺ (Δ⁺⁺) system, from splitting the ROOT files to
plots of the results. One command does the whole chain.

This is the η′π⁻π⁺p version of `simtest_oct2026/amptools_automation`.
That folder is unchanged.

For what every output file means, see [README_output.md](README_output.md).

---

## 1. What the input files must look like

ROOT files in the usual AmpTools `ROOTDataReader` format (tree `kin`), with
**four** final-state particles in this order, like
`pippimetapr/Simulation/Zlm/gen_pol0.root`:

| index | particle |
|---|---|
| 0 | proton |
| 1 | η′ |
| 2 | π⁻ |
| 3 | π⁺ |

The config template uses `reaction EtaPi Beam Proton EtaPrime Pi- Pi+`.
If your trees use another order, change that line and the
`BreitWigner deltapp 1 1 4` lines in the template (1 = proton, 4 = π⁺,
counting the beam as 0).

## 2. What is different from the η π⁻ version

* **Binning.** `split_mass` adds up every particle after the recoil, which
  here would give M(η′π⁻π⁺). The new ROOT macro `split_bins.C` bins in
  **M(η′π⁻)** (particles 1,2) and **−t = −(beam − η′π⁻)²**, the momentum
  transfer to the pπ⁺ system, in a single pass. It is used by default.
* **Fit template.** The template has the new reaction line, and every wave is
  multiplied by the same pπ⁺ Breit-Wigner (`BreitWigner deltapp 1 1 4`, the
  Δ⁺⁺ from your generator config). This lets the fit describe M(pπ⁺) too.
  Delete those lines if you don't want it.
* **Fit-quality plots.** `etapi_fit_plotter.py` now shows **M(pπ⁺)**. It also
  understands the BreitWigner factor and uses pπ⁺ as the recoil for −t and
  the GJ angles. The Fermi-momentum plot is gone because this is a proton target.

## 3. One-time setup (compiling)

Log in to the farm and set up your GlueX environment the way you usually do,
so that `fit`, `root`, `root-config` and `$AMPTOOLS` work. Then go to this folder:

```tcsh
cd /path/to/pippimetapr/amptools_automation
```

**a) `extract_fit_results_etapi`** (required; it reads the .fit files for
the results table). `run_all.csh` builds it automatically the first time.
To build it by hand:

```tcsh
g++ -O2 extract_fit_results_etapi.cc -I$AMPTOOLS -L$AMPTOOLS/lib -lAmpTools \
    `root-config --cflags --libs` -o extract_fit_results_etapi
```

**b) Python packages** for the plots (also installed by `run_all.csh`):

```tcsh
python3 -m pip install --user uproot awkward matplotlib pypdf
```

(`pypdf` is optional. It merges the per-bin PDFs into one file.)

**c) Optional: the C++ plotter `etapi_delta_plotter`.** The Python plotter
does the same job, so you only need this to cross-check. It uses halld_sim's
`EtaPiDeltaPlotGenerator`. The library folder name depends on your build, so
check `ls $HALLD_SIM_HOME`:

```tcsh
g++ -O2 etapi_delta_plotter.cc -I$AMPTOOLS -I$HALLD_SIM_HOME/src/libraries \
    -L$HALLD_SIM_HOME/$BMS_OSNAME/lib -lAMPTOOLS_DATAIO -lAMPTOOLS_AMPS \
    -L$AMPTOOLS/lib -lAmpTools `root-config --cflags --libs` -o etapi_delta_plotter
```

`split_bins.C` and `extract_bin_info_etapi.C` are ROOT macros. They run
through `root` and need no compiling.

## 4. Running it

1. Copy the example settings and edit the file paths:

   ```tcsh
   cp example_run.ini my_run.ini
   ```

   In `my_run.ini`, at minimum set `DATAFILE`, `ACCMCFILE` and `GENMCFILE`,
   plus the `mass` and `t` binning.

2. Run everything:

   ```tcsh
   ./run_all.csh my_run.ini
   ```

   That is the same as `python3 amptools_auto.py all -c my_run.ini`, plus the
   one-time build and package check. It splits the files, fits every bin,
   collects the results into one CSV and makes the plots.

Results go to the `workdir` folder from the INI (default `EtaPrPi_fit/`),
created where you run the command.

## 5. Use-case examples

**Generate a test sample with gen_amp.** `gen_EtaPiMinusDeltaLowerVertex_pol0.cfg`
makes η′π⁻Δ⁺⁺ events with a P1+ wave (π₁ Breit-Wigner) and a D1+ wave (a₂(1320)),
and with Δ⁺⁺ → pπ⁺ at the lower vertex. The trees come out in the order this
fit expects. The `-f` option generates the matching flat sample (no physics) for MC:

```tcsh
gen_amp -c gen_EtaPiMinusDeltaLowerVertex_pol0.cfg -o gen_pol0.root  -l 1.1 -u 2.6 -n 100000 -a 8.0 -b 9.0
gen_amp -c gen_EtaPiMinusDeltaLowerVertex_pol0.cfg -o gen_flat.root -l 1.1 -u 2.6 -n 500000 -a 8.0 -b 9.0 -f
```

`-l`/`-u` set the M(η′π⁻) range, `-n` the number of events and `-a`/`-b` the
beam-energy range. To change the physics, edit the `initialize` lines; for
example, set the negative-reflectivity waves to non-zero values. Fitting
`gen_pol0.root` with the wave set `P1D1 = D1+ P1+` should then give back the
input waves.


**Fit the simulated `gen_pol0.root` as "data".** Here the flat sample serves
as both generated and accepted MC, so the fit is a closure test:

```ini
[samples]
DATAFILE  = /path/to/pippimetapr/Simulation/Zlm/gen_pol0.root
ACCMCFILE = /path/to/pippimetapr/Simulation/Zlm/gen_flat.root
GENMCFILE = /path/to/pippimetapr/Simulation/Zlm/gen_flat.root
```

**Change the binning** (M(η′π⁻) from 1.1 to 2.6 GeV in 30 bins, three −t bins):

```ini
[bin]
mass = 1.1 2.6 30
t = 0.1 1.0 3
```

Leave `t =` empty for mass bins only. Set `t_log = true` for logarithmic −t bins.

**Compare wave sets.** Uncomment or add lines in `[wavesets]`. The listed
waves float, all others are fixed to 0, and the first wave of each
reflectivity is the real reference:

```ini
[wavesets]
P1D1 = D1+ P1+          # what gen_pol0.root was generated with
SPD  = D1 S0 P1 D0 D2   # D1 = D1+ and D1-
ALLp = D1+ S0+ P*+ D*+  # every positive-reflectivity wave
```

Every set gets its own folder, and `plots_wavesets/` overlays them.

**Run only one step or one wave set:**

```tcsh
python3 amptools_auto.py bin     -c my_run.ini      # just split the files
python3 amptools_auto.py fit     -c my_run.ini -r 20 --waveset P1D1
python3 amptools_auto.py gather  -c my_run.ini
python3 amptools_auto.py plot    -c my_run.ini
python3 amptools_auto.py plotfit -c my_run.ini      # data vs. fit, incl. M(p pi+)
```

Options on the command line override the INI. Add `--dry-run` to print the
commands without running them. `python3 amptools_auto.py all -h` lists every option.

**Plot data vs. fit for a single bin by hand**, with a check that the plotter
reproduces the fit's own normalization integrals:

```tcsh
cd EtaPrPi_fit/P1D1/bin_m010_t00
python3 ../../../etapi_fit_plotter.py bin_m010_t00.fit --check
```

This should print `OK`. A `MISMATCH` means the plotter's amplitudes or files
differ from what the fit used.

**Use the fit results in a spreadsheet or your own script.** Open
`EtaPrPi_fit/results.csv` (or `results_wavesets.csv`). It has one row per bin.

**Split a file on its own** (for example to look at one bin in ROOT):

```tcsh
root -l -b -q 'split_bins.C("gen_pol0.root","out/gen",1.1,2.6,30, 0.1,0.3,1, false)'
```

This writes `out/gen_m<i>_t<j>.root`.

## 6. The files in this folder

| file | what it is |
|---|---|
| `run_all.csh` / `run_all.sh` | one command for everything (tcsh / bash) |
| `amptools_auto.py` | the driver: bin → fit → gather → plot → plotfit |
| `example_run.ini` | example settings; copy and edit |
| `gen_EtaPiMinusDeltaLowerVertex_pol0.cfg` | gen_amp config for test samples (P1+ and D1+ waves, Δ⁺⁺ lower vertex) |
| `fit_etaprpim_SPD_TEMPLATE.cfg` | AmpTools config template (S, P, D waves, both reflectivities, Δ⁺⁺ Breit-Wigner) |
| `split_bins.C` | ROOT macro that splits trees into M(η′π⁻) × −t bins |
| `convert_to_csv.py` | turns .fit files (or ROOT files) into CSV |
| `extract_fit_results_etapi.cc` | C++ helper used by `convert_to_csv.py` (compile once) |
| `extract_bin_info_etapi.C` | ROOT macro for per-bin event counts and masses |
| `etapi_fit_plotter.py` | data vs. fit plots (Python, no compiling) |
| `etapi_delta_plotter.cc` | the same plots in C++ (optional) |

## 7. Notes

* The template keeps the reaction name `EtaPi` and the same wave names as the
  η π⁻ version, so `extract_fit_results_etapi` works unchanged.
* `Zlm` in halld_sim takes the angles from particles 1–3 (proton as the
  recoil, η′π⁻ as the resonance) and ignores the π⁺. Your generator config
  did the same, so the fit is consistent with how the samples were made.
* `polAngle` (0°) and `polVal` (0.3519) are set with `define` in the template.
  Change them there for other orientations.
* The whole chain was tested end to end with real AmpTools `fit` (built from
  AmpTools and halld_sim sources) on `gen_pol0.root` with `gen_flat.root` as MC.
  The Python plotter's `--check` matched the fit's normalization integrals to
  about 1e-15.
