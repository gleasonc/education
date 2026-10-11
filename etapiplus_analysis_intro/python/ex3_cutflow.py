"""Exercise 3: Apply your own cuts and make a cut flow.

In C++ we looped over the entries one by one. With numpy we work on whole
columns at once: a comparison such as  chi2 < 10  gives an array of
True/False, one per entry, and we combine those with  &  ("and").
numpy does the loop for us (and much faster).

Run from the tutorial/ folder:
    python python/ex3_cutflow.py
Output: a cut-flow table in the terminal and plots/py_ex3_cutflow.png

Try it: change the numbers in YOUR CUTS and run again.
"""
import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

filename = sys.argv[1] if len(sys.argv) > 1 else "../flat_pipneta__B4_v6.root"
os.makedirs("plots", exist_ok=True)

# ---- YOUR CUTS ------------------------------------------------------
max_chi2ndf = 10.0    # keep Chi2_NDF below this
max_minus_t = 1.0     # keep -t below this (GeV^2)
max_abs_mm2 = 0.05    # keep |MissingMassSquared| below this (GeV^2)

tree = uproot.open(filename)["pipneta__B4"]
branches = ["PiPlus_Px", "PiPlus_Py", "PiPlus_Pz", "PiPlus_E",
            "Eta_Px", "Eta_Py", "Eta_Pz", "Eta_E",
            "Chi2_NDF", "Mandelstam_t", "MissingMassSquared", "accidweight"]
print(f"Reading {tree.num_entries} entries...")
d = tree.arrays(branches, library="np")

# Invariant mass of the eta pi+ pair: M^2 = (E1+E2)^2 - |p1+p2|^2
E  = d["PiPlus_E"]  + d["Eta_E"]
px = d["PiPlus_Px"] + d["Eta_Px"]
py = d["PiPlus_Py"] + d["Eta_Py"]
pz = d["PiPlus_Pz"] + d["Eta_Pz"]
mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
w = d["accidweight"]

# Each step keeps what the previous step kept AND passes one more cut.
step0 = np.ones(len(w), dtype=bool)                        # everything
step1 = step0 & (d["Chi2_NDF"] < max_chi2ndf)
step2 = step1 & (-d["Mandelstam_t"] < max_minus_t)
step3 = step2 & (np.abs(d["MissingMassSquared"]) < max_abs_mm2)
steps = [step0, step1, step2, step3]
names = ["All entries", "Chi2_NDF cut", "-t cut", "|MM2| cut"]
labels = ["All entries", rf"$\chi^2$/ndf < {max_chi2ndf:g}",
          rf"$-t$ < {max_minus_t:g} GeV$^2$", rf"|MM$^2$| < {max_abs_mm2:g} GeV$^2$"]

print(f"\n{'Step':26s}{'Entries':>12s}{'% kept':>10s}{'Weighted':>14s}")
weighted = []
for name, keep in zip(names, steps):
    n = keep.sum()
    sw = w[keep].sum()                     # sum of weights = accidental-subtracted count
    weighted.append(sw)
    print(f"{name:26s}{n:12d}{100 * n / len(w):9.1f}%{sw:14.0f}")

# Mass histogram after all cuts, weighted. The error on a weighted bin is
# sqrt(sum of weights squared) - the same thing ROOT's Sumw2() does.
edges = np.linspace(0.7, 2.5, 91)                     # 90 bins of 20 MeV
counts, _ = np.histogram(mass[step3], bins=edges, weights=w[step3])
sumw2, _ = np.histogram(mass[step3], bins=edges, weights=w[step3] ** 2)
errors = np.sqrt(sumw2)
centres = 0.5 * (edges[1:] + edges[:-1])
print(f"\nWeighted entries in the mass plot: {counts.sum():.0f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
ax1.bar(range(4), weighted, width=1.0, color="lightsteelblue", edgecolor="navy")
ax1.set_xticks(range(4), labels, rotation=20)
ax1.set_yscale("log")
ax1.set_ylim(1e4, None)
ax1.set_ylabel("Weighted entries")
ax1.set_title("Cut flow (accidental-subtracted)")

ax2.errorbar(centres, counts, yerr=errors, fmt="o", ms=3, color="black")
ax2.set_title(r"$\eta\pi^+$ mass after all cuts")
ax2.set_xlabel(r"M($\eta\pi^+$) (GeV/$c^2$)")
ax2.set_ylabel("Weighted entries / 20 MeV")
ax2.set_ylim(0, None)

fig.tight_layout()
fig.savefig("plots/py_ex3_cutflow.png", dpi=100)
print("Saved plots/py_ex3_cutflow.png")
