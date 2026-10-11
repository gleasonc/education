"""Exercise 2: Draw labelled histograms from the tree.

Left: the photon beam energy. Right: the kinematic-fit chi2/ndf for all
entries and for entries with Chi2_NDF < 10 (log y-axis).

Run from the tutorial/ folder:
    python python/ex2_draw.py
Output: a few counts in the terminal and plots/py_ex2_draw.png
"""
import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")            # save to file without opening a window
import matplotlib.pyplot as plt

filename = sys.argv[1] if len(sys.argv) > 1 else "../flat_pipneta__B4_v6.root"
os.makedirs("plots", exist_ok=True)

tree = uproot.open(filename)["pipneta__B4"]
d = tree.arrays(["Beam_E", "Chi2_NDF"], library="np")

# np.histogram counts entries in each bin: same bins as the C++ version.
beam_counts, beam_edges = np.histogram(d["Beam_E"], bins=200, range=(8.0, 12.0))
chi_all, chi_edges = np.histogram(d["Chi2_NDF"], bins=100, range=(0.0, 50.0))
cut = d["Chi2_NDF"] < 10                      # True/False for every entry
chi_cut, _ = np.histogram(d["Chi2_NDF"][cut], bins=100, range=(0.0, 50.0))

print("Entries in hBeam:   ", len(d["Beam_E"]))
print("Entries in hChi2All:", len(d["Chi2_NDF"]), f" (inside 0-50: {chi_all.sum()})")
print("Entries in hChi2Cut:", cut.sum())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.stairs(beam_counts, beam_edges, color="tab:blue", linewidth=1.5)
ax1.set_title("Photon beam energy")
ax1.set_xlabel(r"$E_\gamma$ (GeV)")
ax1.set_ylabel("Entries / 20 MeV")

ax2.stairs(chi_all, chi_edges, color="black", linewidth=1.5, label="All entries")
ax2.stairs(chi_cut, chi_edges, color="tab:red", fill=True, alpha=0.4,
           label=r"$\chi^2$/ndf < 10")
ax2.set_yscale("log")
ax2.set_title(r"Kinematic-fit $\chi^2$/ndf")
ax2.set_xlabel(r"$\chi^2$/ndf")
ax2.set_ylabel("Entries / 0.5")
ax2.legend()

fig.tight_layout()
fig.savefig("plots/py_ex2_draw.png", dpi=100)
print("Saved plots/py_ex2_draw.png")
