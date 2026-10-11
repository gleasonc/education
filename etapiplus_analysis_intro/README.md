# GlueX analysis tutorial: γp → ηπ⁺n

A first hands-on analysis for new undergraduate researchers, in both **C++ (ROOT)** and **Python**. You will open a real GlueX analysis tree, find your way around it, make labelled plots, apply your own cuts, subtract accidental beam photons, look at 2D distributions, fit a meson peak, and measure decay angles in the Gottfried-Jackson frame.

It fills the Tuesday and Thursday slots of the bootcamp week in the *GlueX New Mentor Guides* (Guide 2): "open a sample ROOT file, fill and fit a histogram in Python and in ROOT", then "explain each cut". Plan on one morning for exercises 1 to 3 and one afternoon for 4 to 6. Exercises 7 and 8 (2D mass plots and Gottfried-Jackson angles) are the next step once those work.

**Done when you can:** make a labelled histogram from a tree you have not seen before, apply a cut of your own and say what it did, and fit a peak and explain what the numbers mean.

> **New to programming?** Start with the macros in [`cpp_beginner/`](cpp_beginner/README.md). They do the same analysis in the plainer style of the group's own macros (like `MakeFlat.C`): one variable per branch, one event loop, cuts written as `if (...) continue;`. The macros in `cpp/` are more compact and use a few more ROOT shortcuts.
>
> An earlier version of this tutorial, built on the histogram-only file `etapipmissn__U1.root`, is kept in `old_etapipmissn_U1/`.

---

## 0. The physics in one paragraph

A high-energy photon (γ) from the GlueX beam hits a proton (p) in the liquid-hydrogen target. In the reaction studied here, the products are an **η meson**, a **π⁺** and a **neutron**:

γ p → η π⁺ n,  with  η → γ γ

Short-lived mesons that decay to ηπ⁺ show up as **peaks in the ηπ⁺ invariant mass**. Two you will see:

| Meson | Mass (PDG) | Width (PDG) |
|---|---|---|
| a₀(980) | 980 ± 20 MeV | 50 to 100 MeV |
| a₂(1320) | 1318 MeV | about 107 MeV |

---

## 1. What is in the file

`flat_pipneta__B4_v6.root` (810 MB) holds **one TTree** called `pipneta__B4` with **6,091,343 entries**. A tree is a table: each **entry** (row) is one candidate event (one *combo*, see below), and each **branch** (column) is one number about it. "Flat" means every branch holds a single number per entry, so it is easy to read.

| Branch(es) | Type | What it is |
|---|---|---|
| `PiPlus_Px, _Py, _Pz, _E` | double | π⁺ four-momentum (GeV) |
| `Eta_Px, _Py, _Pz, _E` | double | η four-momentum (GeV) |
| `Neutron_Px, _Py, _Pz, _E` | double | neutron four-momentum (GeV) |
| `Beam_Px, _Py, _Pz, _E` | double | beam photon four-momentum (GeV); it travels along z, so Px = Py = 0 |
| `Neutron_Beta` | double | neutron speed as a fraction of the speed of light (0.03 to 0.90) |
| `Mandelstam_t` | double | momentum transfer t (GeV²), always negative here (−2.28 to 0) |
| `Chi2_NDF` | double | kinematic-fit χ² per degree of freedom (0.005 to 19,270) |
| `MissingMassSquared` | double | missing mass² (GeV²), already restricted to −0.1 to 0.1 |
| `accidweight` | double | +1 for in-time beam photons, about −0.17 to −0.28 for out-of-time ones (exercise 4) |
| `pol` | int | polarization angle of the beam in degrees: 0, 45, 90, 135, or −1 for no polarization |
| `deltaPhi` | double | an azimuthal angle difference in degrees (−360 to 360); ask your mentor exactly which |

The tree has no branch for the ηπ⁺ mass; you compute it yourself from the four-momenta:

M(ηπ⁺)² = (E_η + E_π)² − (px_η + px_π)² − (py_η + py_π)² − (pz_η + pz_π)²

A **combo** is one way of assigning the detected tracks, showers and beam photons to the particles in the reaction. One real event can give several combos (for example several beam photons close in time), so entries are combos, not events.

---

## 2. Set up

You need **one** of the following.

**Option A: JLab farm (ifarm).** ROOT and Python are available once you load the GlueX environment. Your mentor will give you the setup line for your group.

**Option B: your laptop, using conda** (Mac, Linux or Windows with WSL):

```bash
conda create -n gluex-tutorial -c conda-forge root uproot matplotlib scipy jupyter
conda activate gluex-tutorial
```

**Option C: Python only** (Jupyter, Google Colab, any Python 3). `uproot` reads ROOT files without ROOT installed:

```bash
pip install uproot matplotlib scipy
```

The Python scripts read the columns they need into memory, about 0.6 GB for the biggest exercises. That is fine on a laptop or on Colab.

### Folder layout

```
UndergraduateResearch/
├── flat_pipneta__B4_v6.root  <- the data file
└── tutorial/                 <- run Python from here, C++ from tutorial/cpp/
    ├── README.md             <- this guide
    ├── tutorial.ipynb        <- Python exercises as a notebook
    ├── cpp/ex1_explore.C ... ex8_gj_angles.C
    ├── cpp_beginner/         <- the same analysis in a plainer style (has its own README)
    ├── python/ex1_explore.py ... ex8_gj_angles.py
    ├── reference/            <- the plots and printout you should get
    └── old_etapipmissn_U1/   <- earlier version (histogram-only file)
```

Open a terminal and go to the tutorial folder (C++ users then go one level further, into `cpp/`):

```bash
cd UndergraduateResearch/tutorial
```

Run the **C++ macros from inside `tutorial/cpp/`** (and the beginner macros from inside `tutorial/cpp_beginner/`), and the **Python scripts from `tutorial/`**. Each script finds the data file relative to where it is run, as in the table below. If your copy is somewhere else, pass the path as an argument. Plots are written to a new `plots/` folder.

| | C++ (ROOT) | Python |
|---|---|---|
| Run from | `tutorial/cpp/` | `tutorial/` |
| Run exercise N | `root -l -b -q exN_name.C` | `python python/exN_name.py` |
| Use another file | `root -l -b -q 'exN_name.C("/path/file.root")'` | `python python/exN_name.py /path/file.root` |
| Data file it looks for | `../../flat_pipneta__B4_v6.root` | `../flat_pipneta__B4_v6.root` |
| Plots saved to | `tutorial/cpp/plots/` | `tutorial/plots/` |

`-l` skips the splash screen, `-b` runs without opening windows, and `-q` quits when the macro ends. Leave out `-b -q` to see the plots on screen. Each C++ macro defines one function with the same name as the file (`ex1_explore.C` defines `ex1_explore`), which is the function ROOT runs.

**Read the C++ macros, not only their output.** They are commented line by line for people who have not programmed before: what a pointer is, what `->` does, how a loop works, what a function argument is. Exercises 1 and 3 explain the most. If they still feel dense, use `cpp_beginner/` first.

Try the C++ and the Python version of each exercise and check that they print the same numbers. They do in `reference/expected_output.txt`.

---

## Exercise 1: What is inside the file?

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex1_explore.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex1_explore.py
```

The script lists the file contents, the 23 branches with their types, the first 5 entries of a few branches, and some counts.

**You should see:** 6,091,343 entries; 1,763,934 with `accidweight == 1` and 4,327,409 with `accidweight < 0`; 699,087 with `Chi2_NDF < 10`; and per `pol` value: −1: 507,440, 0: 1,390,581, 45: 1,365,837, 90: 1,439,007, 135: 1,388,478.

**Ideas to take away**
- A **TTree** is a table: entries are rows, branches are columns.
- `/D` in a branch title means a decimal number (double), `/I` a whole number (integer).
- `tree->Scan(...)` prints rows; `tree->GetEntries("condition")` counts rows where the condition is true.

**Questions**
1. Rows 3 and 4 have `Chi2_NDF` above 500. Do you think those combos are well described by γp → ηπ⁺n?
2. Why are there more out-of-time entries (`accidweight < 0`) than in-time ones?

---

## Exercise 2: Draw labelled histograms

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex2_draw.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex2_draw.py
```

Output: `plots/cpp_ex2_draw.png` and `plots/py_ex2_draw.png`. Left: the beam energy. Right: `Chi2_NDF` for all entries and for `Chi2_NDF < 10`, on a log y-axis. Of the 6,091,343 entries, 1,629,202 have `Chi2_NDF` between 0 and 50.

This uses **`TTree::Draw`**, the quickest way to histogram a branch:

```cpp
tree->Draw("Beam_E>>hBeam", "cut", "goff");   // fill histogram hBeam with Beam_E, only where cut is true
```

**Ideas to take away**
- Make the empty histogram first (`new TH1D(name, "title;x label;y label", bins, low, high)`), then fill it.
- A plot without axis labels and units cannot be read by anyone else, including you in a month.
- A log y-axis shows big and small numbers on the same plot.

**Questions**
1. The beam energy has a bump between 8.2 and 8.8 GeV. This is the *coherent peak*, where the beam is linearly polarized. Why does the plot start sharply at 8.2 GeV?
2. Most entries have `Chi2_NDF` far above 10. What does that tell you about most of the combos in the file?

---

## Exercise 3: Your own cuts and a cut flow

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex3_cutflow.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex3_cutflow.py
```

The C++ macro uses an **event loop**: it connects each branch to a variable with `SetBranchAddress`, then reads the entries one at a time with `GetEntry(i)`, computes M(ηπ⁺) with `TLorentzVector`, and applies three cuts. The Python version does the same with numpy arrays. The cut values are at the top of each file under **YOUR CUTS**.

**You should see:**

```
Step                           Entries    % kept      Weighted
All entries                    6091343    100.0%        930757
Chi2_NDF cut                    699087     11.5%        177950
-t cut                          433403      7.1%        109496
|MM2| cut                       362187      5.9%         98492

Weighted entries in the mass plot: 76108
```

The cuts are `Chi2_NDF < 10`, `−t < 1.0 GeV²` and `|MissingMassSquared| < 0.05 GeV²`. "Weighted" is the sum of `accidweight`, which is the number of entries after subtracting accidental beam photons (exercise 4). The mass plot has fewer weighted entries (76,108) than the last step (98,492) because it only shows 0.7 to 2.5 GeV, and M(ηπ⁺) goes up to 3.8 GeV.

**Ideas to take away**
- A **cut** keeps entries where a condition is true. Applying cuts one after another and counting after each one gives the **cut flow**.
- `continue` in a loop skips to the next entry, which is how a failed cut throws an entry away.
- After the cuts, two peaks appear near 1.0 and 1.3 GeV: the a₀(980) and a₂(1320).

**Questions**
1. Which cut removes the most? Change `maxChi2NDF` to 5 and to 20. What happens to the peaks and to the background?
2. Add a fourth cut of your own (ideas: `Beam_E` in the coherent peak 8.2 to 8.8 GeV, or one `pol` value). Put it in the cut flow.

---

## Exercise 4: What does a cut do? And what are accidentals?

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex4_compare.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex4_compare.py
```

Left: M(ηπ⁺) with no cut and with `Chi2_NDF < 10`, each scaled to area 1. Right: after the χ² cut, in-time entries, the accidental estimate, and their difference.

**You should see:**

```
Weighted sum, no cut:              688928
Weighted sum, chi2/ndf < 10:       119942
Fraction kept by the cut:          0.1741
Prompt entries (chi2/ndf < 10):    176715
Accidental estimate:               56772.6
Prompt - accidental:               119942
```

This macro introduces two `TTree::Draw` tricks: `tree->SetAlias("M_etapi", "sqrt(...)")` gives the mass formula a nickname, and the second argument of `Draw` can be a **weight** such as `"accidweight * (Chi2_NDF < 10)"` (a true/false condition counts as 1 or 0).

**Ideas to take away**
- **Normalizing** to the same area lets you compare *shapes* when the totals are very different.
- **Accidentals:** the beam photons arrive in bunches every 4 ns. Sometimes the photon matched to an event came from a different bunch than the one that made the event. The **in-time** (prompt) sample contains some of these. Combos with photons from **out-of-time** bunches estimate how many, so they get a negative weight; filling with `accidweight` subtracts them.
- The accidental estimate still shows small a₀ and a₂ peaks: those are real events paired with the wrong beam photon.

**Questions**
1. Which mass region does the χ² cut keep, and which does it remove?
2. Why is it fine for a histogram bin to have a negative weighted content, and when would that worry you?

---

## Exercise 5: A 2D plot and its projections

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex5_2d.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex5_2d.py
```

Left: a 2D histogram of M(ηπ⁺) (x) against −t (y), with the χ² and missing-mass cuts, weighted. Right: **projections** onto the mass axis for 0 < −t < 0.5 GeV² and 0.5 < −t < 1.0 GeV².

**You should see:** 104,480 weighted entries in the 2D plot, 32,873 in the low −t slice and 43,235 in the high −t slice.

**Ideas to take away**
- In `TTree::Draw`, a 2D plot is written `"y:x"` (y first).
- `ProjectionX(name, firstYbin, lastYbin)` adds up a range of y bins and gives a 1D histogram. In Python, sum the 2D array over those bins.
- The a₂(1320) is relatively stronger at small −t, the a₀(980) at larger −t. How a peak changes with −t tells you about how the meson was produced.

**Questions**
1. There are very few entries below −t ≈ 0.2 GeV² after the cuts, even though a third of all entries in the file have −t < 0.2. Which cut removes them? (Hint: rerun exercise 3 with the χ² cut loosened.) Why might the kinematic fit struggle when the neutron is slow?
2. Make the same 2D plot of M(ηπ⁺) against `Beam_E`. Do the peaks move?

---

## Exercise 6: Fit the a₀(980) peak

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex6_fit.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex6_fit.py
```

The script makes M(ηπ⁺) with the exercise 3 cuts in 10 MeV bins and fits 0.88 to 1.10 GeV with a **Gaussian** for the peak plus a **straight line** for the background.

**You should see:**

| Quantity | C++ (ROOT/MINUIT) | Python (scipy) |
|---|---|---|
| Mean | 985.7 ± 1.1 MeV | 985.7 ± 1.2 MeV |
| σ (Gaussian width) | 26.4 ± 1.6 MeV | 26.4 ± 1.3 MeV |
| χ²/ndf | 31.9 / 17 (p-value 0.016) | 31.9 / 17 (p-value 0.016) |
| Signal yield | about 3,320 weighted entries | about 3,320 weighted entries |

The fitted values agree. The uncertainty on σ differs slightly because MINUIT and scipy estimate uncertainties in different ways (MINUIT uses the full curvature of χ², scipy an approximation of it). Both are approximations; the difference is a reminder not to quote more digits than you trust.

**Ideas to take away**
- A fit needs **starting values** near the answer, or it can wander off. Look at the plot first and guess.
- **χ²/ndf** near 1 means the model describes the data about as well as the error bars allow. Here it is 1.9 (p-value 0.016): the model is missing something. The a₀(980) sits right at the K K̄ threshold (987 MeV), which gives it a sharp, non-Gaussian shape.
- With weighted entries, the error on a bin is √(sum of weights²), not √N. ROOT does this when you call `Sumw2()`.
- The **yield** is the area under the Gaussian divided by the bin width.
- Changing the fit range, background shape or signal shape moves the answer. How much it moves is a **systematic uncertainty**.

**Questions**
1. Is the fitted mean consistent with the PDG a₀(980) mass? How many standard deviations apart are they?
2. Change the fit range to 0.85 to 1.15 GeV (C++, from `cpp/`: `root -l -b -q 'ex6_fit.C("../../flat_pipneta__B4_v6.root", 0.85, 1.15)'`; Python: `python python/ex6_fit.py ../flat_pipneta__B4_v6.root 0.85 1.15`). How much do the mean, width and yield change?

---

## Exercise 7: Two-dimensional mass plots

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex7_masses2d.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex7_masses2d.py
```

With three particles in the final state there are three pairs, so three two-body masses: M(ηπ⁺), M(ηn) and M(π⁺n). After the standard cuts, the macro plots them against each other, plus a **Dalitz plot** (M²(ηπ⁺) against M²(π⁺n)), all weighted by `accidweight`. It then applies a **baryon veto**, keeping only M(π⁺n) > 2 GeV and M(ηn) > 2 GeV, and compares M(ηπ⁺) with and without it.

**You should see:**

```
Weighted entries in the Dalitz plot:          96168.3
M(eta pi+) 0.7-2.5 GeV, standard cuts:        76107.8
M(eta pi+) 0.7-2.5 GeV, plus baryon veto:     40348.1
Fraction kept by the veto:                    0.5301
```

**Ideas to take away**
- M(ηπ⁺) is on the y-axis, so a meson that decays to ηπ⁺ makes a **horizontal band**: the a₀(980) at 0.98 GeV and the a₂(1320) at 1.32 GeV. In the Dalitz plot these sit at M² ≈ 0.96 and 1.74 GeV².
- A baryon makes a **vertical band**: the Δ(1232) → π⁺n near M(π⁺n) = 1.23 GeV, and an ηn enhancement right at threshold, M(ηn) ≈ 1.5 GeV, where the N(1535) sits.
- Every band in one plot shows up as a curved band in the others. These **reflections** are why the M(ηπ⁺) background rises at high mass.
- The edges of each plot are set by energy and momentum conservation (the kinematic boundary).
- The veto keeps 53% of the weighted entries. It removes most of the rising high-mass background while keeping the a₀ and a₂ peaks.

**Questions**
1. Where does the Δ(1232) band cross the a₂(1320) band? What could happen to the a₂ peak there?
2. Try the veto with 1.5 GeV instead of 2 GeV. What changes in the M(ηπ⁺) plot?

---

## Exercise 8: Gottfried-Jackson angles

```bash
cd tutorial/cpp            # C++: run from inside cpp/
root -l -b -q ex8_gj_angles.C
```

```bash
cd tutorial                # Python: run from tutorial/
python python/ex8_gj_angles.py
```

A meson's **spin** shows up in the angle at which it decays. A spin-0 a₀(980) decays the same way in every direction, but a spin-2 a₂(1320) does not. To measure that angle, you look at the decay in the meson's own rest frame, with axes fixed by the reaction: the **Gottfried-Jackson (GJ) frame**. The steps are the same as in the group's `MakeFlat.C`:

1. **Lab → centre-of-mass (CM) frame** of beam + target (target = proton at rest, (0, 0, 0, 0.938272) GeV). A **boost** rewrites every four-vector as seen by an observer moving at some velocity. To reach the rest frame of something with momentum p and energy E, boost by −p/E.
2. **CM → ηπ⁺ rest frame**, boosting by −p/E of the ηπ⁺ system (in the CM frame).
3. **Axes:** z = beam direction in the ηπ⁺ rest frame; y = beam × ηπ⁺ in the CM frame (perpendicular to the production plane); x = y × z.
4. **Angles of the η:** cos θ_GJ is its direction along z, and φ_GJ is its angle around z measured from x.

As a **cross-check**, the script computes cos θ_GJ a second way: it boosts the η and the beam straight from the lab into the ηπ⁺ rest frame and takes the angle between them. Rotations don't change angles between two vectors, so this must give the same cos θ_GJ for every entry. It does not check φ_GJ, which depends on how the x and y axes were chosen.

**You should see:**

```
Entries used:                    362187
cos(theta_GJ) range:             -0.999507 to 0.999979
phi_GJ range (deg):              -179.998 to 179.999
Max |cos - cross-check|:         1.25e-12
Weighted entries, all masses:    98492.2
Weighted entries, a0 window:     8314.83  (mean cos = -0.0508)
Weighted entries, a2 window:     8535.99  (mean cos = -0.0746)
```

C++ and Python agree. The mass windows are 0.92 to 1.04 GeV (a₀) and 1.24 to 1.40 GeV (a₂).

**Ideas to take away**
- The big peaks at cos θ_GJ near ±1 come from high M(ηπ⁺) (see the 2D plot). Those are mostly the baryon reflections from exercise 7, where the η or π⁺ goes along with the neutron.
- The a₀ and a₂ windows look similar here. That is **not** evidence that both have spin 0. Each window also holds background: the exercise 6 fit gives about 3,300 a₀ entries, so more than half of the 8,315 in the a₀ window are background, and the detector doesn't see every angle equally well (its **acceptance**): see the drop in the first bin at cos θ = −1. Measuring spin properly needs acceptance-corrected fits (partial-wave analysis, for example with AmpTools, which is what the rest of `MakeFlat.C` prepares files for).
- φ_GJ is not flat either: there is a dip near 0° and a rise towards ±180°.

**Questions**
1. Add the baryon veto from exercise 7. How do the cos θ_GJ plots change?
2. Make cos θ_GJ for a sideband (for example 1.45 to 1.60 GeV) and subtract it from the a₂ window. Does the a₂ shape change?

---

## 7. Challenges

1. **Quadratic background.** C++: change `pol1(3)` to `pol2(3)` and add a starting value. Python: add a `b2 * x**2` term. Does the yield change by more than its statistical error?
2. **The a₂(1320).** Fit the second peak, 1.15 to 1.50 GeV. Compare the mass and width with the PDG.
3. **Breit-Wigner.** Replace the Gaussian with a Breit-Wigner (ROOT: `TMath::BreitWigner`).
4. **Other combinations.** Compute M(π⁺n) and M(ηn) from the four-momenta. Do you see the Δ⁺(1232) or N* resonances?
5. **Polarization.** Make the mass plot separately for each `pol` value. Should the peaks look different?
6. **deltaPhi.** Plot `deltaPhi`. Most entries sit in a narrow peak at 0. Try a cut such as |deltaPhi| < 20° and see what it does to the peaks and to the background. Before using it, find out what the angle is.
7. **PyROOT.** In Python, `import ROOT` gives the same classes as the C++ macros (`ROOT.TFile.Open`, `tree.Draw`, `h.Fit`, ...). Rewrite exercise 2 with PyROOT.

---

## 8. Check the file before you trust it

Real analysis files have rough edges. These are facts about this file worth asking your mentor about:

- **Fixed masses.** The η mass is exactly 0.5479 GeV and the neutron mass exactly 0.9396 GeV in every entry. The four-momenta are the kinematic-fit output with those masses constrained, so plotting the η or neutron mass tells you nothing.
- **Pre-applied cuts.** `MissingMassSquared` stops sharply at ±0.1 GeV², `Neutron_Beta` at 0.9, and `Beam_E` starts at 8.2 GeV. These cuts were applied when the tree was made.
- **Missing mass² peaks at 0**, not at the neutron mass² (0.883 GeV²). That suggests the neutron was measured (it also has a `Neutron_Beta`) and the missing mass is for the whole reaction. The older file in `old_etapipmissn_U1/` treated the neutron as missing.
- **`Chi2_NDF` is very large for most entries:** median 210, mean 486, and 16% above 1,000. Check whether this branch is really χ² per degree of freedom.
- **Accidental weights.** Out-of-time weights run from −1/6 (−0.167) to −0.283, not one fixed value, and there are 2.45 out-of-time entries per in-time one. Ask how many beam bunches were used (`B4` in the name) and whether a scale factor is included.
- **`deltaPhi` runs from −360° to +360°**, so it is the difference of two angles without wrapping into ±180°. Values near ±360° are actually close to 0°.
- **No run or event number.** You cannot tell which combos belong to the same event, so you cannot pick a best combo per event.

None of this stops the exercises, but note each in your log and raise it at your next meeting.

---

## 9. C++ and Python side by side

| Task | C++ (ROOT) | Python (uproot + numpy + matplotlib) |
|---|---|---|
| Open a file | `TFile* f = TFile::Open("file.root");` | `f = uproot.open("file.root")` |
| Get the tree | `TTree* t = (TTree*)f->Get("pipneta__B4");` | `t = f["pipneta__B4"]` |
| Number of entries | `t->GetEntries()` | `t.num_entries` |
| List branches | `t->Print()` | `t.keys()`, `t.typenames()` |
| Read branches | `t->SetBranchAddress("Beam_E", &beamE); t->GetEntry(i);` | `d = t.arrays(["Beam_E"], library="np")` |
| Quick histogram | `t->Draw("Beam_E>>h", "Chi2_NDF<10")` | `np.histogram(d["Beam_E"][d["Chi2_NDF"] < 10], bins=..., range=...)` |
| Weighted fill | `t->Draw("x>>h", "accidweight*(cut)")` or `h->Fill(x, w)` | `np.histogram(x, bins, weights=w)` |
| Weighted errors | `h->Sumw2();` | `np.sqrt(np.histogram(x, bins, weights=w**2)[0])` |
| Normalize | `h->Scale(1.0 / h->Integral());` | `counts / counts.sum()` |
| 2D plot | `t->Draw("y:x>>h2", ...); h2->Draw("colz");` | `np.histogram2d(x, y, ...)`; `ax.pcolormesh(xe, ye, h2.T)` |
| Project onto x | `h2->ProjectionX("p", ybinLo, ybinHi)` | `h2[:, ymask].sum(axis=1)` |
| Fit | `h->Fit(tf1, "RS");` | `scipy.optimize.curve_fit(model, x, y, sigma=err)` |
| Save the plot | `c->SaveAs("plot.png");` | `fig.savefig("plot.png")` |

---

## 10. What to bring to your mentor

- Your plots, with labelled axes, from at least one of the two languages (or from `cpp_beginner/`).
- Your cut flow with at least one cut you chose yourself, and what it did to the a₀ and a₂ peaks.
- Your fitted a₀(980) mean, width and yield, and how much they moved when you changed the fit range.
- Your answers to the questions, and one thing from section 8 you would like explained.

The `reference/` folder has every plot and the full printed output from a run of both versions, so you can check your own.
