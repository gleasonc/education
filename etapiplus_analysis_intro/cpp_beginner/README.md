# Beginner C++ macros

These macros do the same analysis as the ones in `../cpp/`, written in the plainer style of the group's own analysis macros (like `MakeFlat.C`). Start here if you have not programmed before.

Every macro has the same five parts, in the same order:

1. open the file and get the tree
2. make one variable per branch and connect it with `SetBranchAddress`
3. make empty histograms (`new TH1F(...)`, `new TH2F(...)`)
4. loop over every entry: read it, build four-vectors, apply cuts with `if (...) continue;`, fill histograms
5. draw the histograms on a canvas and save a picture

Once you can read one macro, you can read all of them. Each one adds a little to the one before.

## How to run

```bash
cd tutorial/cpp_beginner
root -l ex1_explore.C
```

Type `.q` to quit ROOT afterwards. Add `-b -q` (`root -l -b -q ex1_explore.C`) to run without opening windows and quit at the end. Pictures are saved in this folder as `.pdf` (for printing) and `.png` (for a quick look).

Each macro opens `../../flat_pipneta__B4_v6.root`, which is the data file two folders up. If your file is somewhere else, change that line near the top.

Exercises 2 to 8 loop over all 6 million entries, which takes about 5 to 15 seconds each.

## The macros

| Macro | What it teaches | What you should see |
|---|---|---|
| `ex1_explore.C` | Open the file, list the branches, print the first 5 entries, `tree->Draw("Beam_E")` | 6091343 entries; `ex1_BeamE.pdf` |
| `ex2_draw.C` | The 5-part layout; fill two histograms in a loop | `ex2_draw.pdf`: beam energy and χ²/NDF |
| `ex3_masses.C` | `TLorentzVector`: add four-vectors, take `.M()` | `ex3_masses.pdf`: M(ηπ⁺), M(ηn), M(π⁺n), no cuts |
| `ex4_cuts.C` | Cuts with `continue`, counters, before/after on one plot | 6091343 → 699087 → 433403 → 362187; `ex4_cuts.pdf` |
| `ex5_2d.C` | 2D histograms (`TH2F`, `colz`), `ProjectionX` | 32873.1 and 43234.7 in the two −t slices; `ex5_2d.pdf` |
| `ex6_fit.C` | Fit a Gaussian + straight line with `TF1` | Mean ≈ 0.9857 GeV, Sigma ≈ 0.0264 GeV, χ²/ndf = 31.86/17; `ex6_fit.pdf` |
| `ex7_masses2d.C` | All three mass-vs-mass plots, a Dalitz plot, a baryon veto | 76107.8 → 40348.1 with the veto; `ex7_masses2d.pdf`, `ex7_veto.pdf` |
| `ex8_gj_angles.C` | Boosts, the Gottfried-Jackson frame, cos θ_GJ and φ_GJ | 8314.83 (a₀ window) and 8536 (a₂ window); `ex8_gj_angles.pdf` |
| `MakePlots_EtaPiPlus.C` | Optional: exercises 2 to 5 in one macro, laid out like `MakeFlat.C` | `ChiSqNDF.pdf`, `NeutronBeta.pdf`, `Masses.pdf`, `EtaPiPlusMass_Cuts.pdf`, `Masses2D.pdf` |

The cuts used throughout are χ²/NDF < 10, −t < 1 GeV² and |MM²| < 0.05 GeV², and histograms are filled with `accidweight` to subtract accidental beam photons. `../README.md` explains the physics, what each plot means, and the questions to answer. The numbers above are the same as the ones from `../cpp/` and `../python/`.

The full printout of a run is in `../reference/beginner_output.txt`, and all the pictures are in `../reference/beginner/`.

## Words you will meet

| Word | Meaning |
|---|---|
| `TFile*`, `TTree*`, `TH1F*` | The `*` means a **pointer**: the location of an object in memory. You use it with `->`, as in `tree->GetEntries()`. |
| `&Beam_E` | "The location of the variable `Beam_E`". `SetBranchAddress` needs it so it knows where to put each value. |
| `for (Long64_t i_entry=0; i_entry<NEntries; i_entry++)` | Repeat the code in `{ }` for i_entry = 0, 1, 2, … up to NEntries − 1. |
| `continue` | Skip the rest of the loop for this entry and go on to the next one. This is how a cut throws an entry away. |
| `&&`, `\|\|`, `==` | "and", "or", "is equal to". |
| `Fill(x, accidweight)` | Add `accidweight` (instead of 1) to the bin that contains x. |
| Boost | Re-write a four-vector as seen by an observer moving at some velocity (exercise 8). |
