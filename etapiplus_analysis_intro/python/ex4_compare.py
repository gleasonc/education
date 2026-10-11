"""Exercise 4: What does a cut do? And what are accidentals?

Left:  eta pi+ mass with no cut and with Chi2_NDF < 10, each scaled to
       area 1 so the SHAPES can be compared.
Right: after the chi2 cut: in-time (prompt) photons, the accidental
       estimate from out-of-time photons, and prompt minus accidental.

Run from the etapiplus_analysis_intro/ folder:
    python python/ex4_compare.py
Output: sums in the terminal and plots/py_ex4_compare.png
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

tree = uproot.open(filename)["pipneta__B4"]
d = tree.arrays(["PiPlus_Px", "PiPlus_Py", "PiPlus_Pz", "PiPlus_E",
                 "Eta_Px", "Eta_Py", "Eta_Pz", "Eta_E",
                 "Chi2_NDF", "accidweight"], library="np")

E  = d["PiPlus_E"]  + d["Eta_E"]
px = d["PiPlus_Px"] + d["Eta_Px"]
py = d["PiPlus_Py"] + d["Eta_Py"]
pz = d["PiPlus_Pz"] + d["Eta_Pz"]
mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
w = d["accidweight"]
chi2cut = d["Chi2_NDF"] < 10
prompt = w == 1
accid = w < 0

edges = np.linspace(0.7, 2.5, 91)

def hist(mask, weights):
    """Weighted histogram of the mass for the entries where mask is True."""
    counts, _ = np.histogram(mass[mask], bins=edges, weights=weights[mask])
    return counts

everything = np.ones(len(w), dtype=bool)
h_nocut  = hist(everything, w)
h_cut    = hist(chi2cut, w)
h_prompt = hist(chi2cut & prompt, np.ones(len(w)))     # each prompt entry counts 1
h_accid  = hist(chi2cut & accid, -w)                   # -weight is positive
h_sub    = hist(chi2cut, w)                            # prompt minus accidental

print(f"Weighted sum, no cut:              {h_nocut.sum():.0f}")
print(f"Weighted sum, chi2/ndf < 10:       {h_cut.sum():.0f}")
print(f"Fraction kept by the cut:          {h_cut.sum() / h_nocut.sum():.4f}")
print(f"Prompt entries (chi2/ndf < 10):    {h_prompt.sum():.0f}")
print(f"Accidental estimate:               {h_accid.sum():.1f}")
print(f"Prompt - accidental:               {h_sub.sum():.0f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
xlabel = r"M($\eta\pi^+$) (GeV/$c^2$)"

ax1.stairs(h_nocut / h_nocut.sum(), edges, color="black", lw=1.5, label="No cut")
ax1.stairs(h_cut / h_cut.sum(), edges, color="tab:red", lw=1.5, label=r"$\chi^2$/ndf < 10")
ax1.set_title(r"Effect of the $\chi^2$/ndf cut")
ax1.set_xlabel(xlabel)
ax1.set_ylabel("Fraction / 20 MeV")
ax1.set_ylim(0, None)
ax1.legend(loc="lower right")

ax2.stairs(h_prompt, edges, color="tab:blue", lw=1.5, label="In-time (prompt)")
ax2.stairs(h_accid, edges, color="gray", fill=True, alpha=0.5, label="Accidental estimate")
ax2.stairs(h_sub, edges, color="tab:red", lw=1.5, label="Prompt − accidental")
ax2.set_title(r"Accidental subtraction ($\chi^2$/ndf < 10)")
ax2.set_xlabel(xlabel)
ax2.set_ylabel("Entries / 20 MeV")
ax2.set_ylim(0, 1.4 * h_prompt.max())
ax2.legend(loc="upper right")

fig.tight_layout()
fig.savefig("plots/py_ex4_compare.png", dpi=100)
print("Saved plots/py_ex4_compare.png")
