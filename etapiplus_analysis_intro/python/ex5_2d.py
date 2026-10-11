"""Exercise 5: A 2D plot and its projections.

Left:  2D histogram of the eta pi+ mass (x) against -t (y).
Right: projections onto the mass axis for 0 < -t < 0.5 and 0.5 < -t < 1.0.
Cuts: Chi2_NDF < 10 and |MissingMassSquared| < 0.05, weighted by accidweight.

Run from the etapiplus_analysis_intro/ folder:
    python python/ex5_2d.py
Output: sums in the terminal and plots/py_ex5_2d.png
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
                 "Chi2_NDF", "MissingMassSquared", "Mandelstam_t", "accidweight"],
                library="np")

E  = d["PiPlus_E"]  + d["Eta_E"]
px = d["PiPlus_Px"] + d["Eta_Px"]
py = d["PiPlus_Py"] + d["Eta_Py"]
pz = d["PiPlus_Pz"] + d["Eta_Pz"]
mass = np.sqrt(E**2 - px**2 - py**2 - pz**2)
minus_t = -d["Mandelstam_t"]
keep = (d["Chi2_NDF"] < 10) & (np.abs(d["MissingMassSquared"]) < 0.05)

x_edges = np.linspace(0.7, 2.5, 91)    # 90 mass bins
y_edges = np.linspace(0.0, 2.0, 41)    # 40 -t bins of 0.05 GeV^2
h2, _, _ = np.histogram2d(mass[keep], minus_t[keep], bins=[x_edges, y_edges],
                          weights=d["accidweight"][keep])
# h2[i, j] = content of mass bin i, -t bin j

# Projection = add up the -t bins inside a range (axis=1 is the -t direction)
y_centres = 0.5 * (y_edges[1:] + y_edges[:-1])
p_low  = h2[:, (y_centres > 0.0) & (y_centres < 0.5)].sum(axis=1)
p_high = h2[:, (y_centres > 0.5) & (y_centres < 1.0)].sum(axis=1)

print(f"Weighted entries in the 2D plot:     {h2.sum():.0f}")
print(f"Weighted entries with 0 < -t < 0.5:   {p_low.sum():.1f}")
print(f"Weighted entries with 0.5 < -t < 1.0: {p_high.sum():.1f}")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
xlabel = r"M($\eta\pi^+$) (GeV/$c^2$)"

# pcolormesh wants [y, x], so transpose with .T; vmin=0 hides negative bins
mesh = ax1.pcolormesh(x_edges, y_edges, h2.T, cmap="viridis", vmin=0)
fig.colorbar(mesh, ax=ax1, label="Weighted entries")
ax1.set_title(r"$\eta\pi^+$ mass vs $-t$")
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"$-t$ (GeV$^2$)")

ax2.stairs(p_low, x_edges, color="tab:blue", lw=1.5, label=r"0 < $-t$ < 0.5 GeV$^2$")
ax2.stairs(p_high, x_edges, color="tab:red", lw=1.5, label=r"0.5 < $-t$ < 1.0 GeV$^2$")
ax2.set_title("Projections onto the mass axis")
ax2.set_xlabel(xlabel)
ax2.set_ylabel("Weighted entries / 20 MeV")
ax2.set_ylim(0, 1.3 * max(p_low.max(), p_high.max()))
ax2.legend(loc="upper right")

fig.tight_layout()
fig.savefig("plots/py_ex5_2d.png", dpi=100)
print("Saved plots/py_ex5_2d.png")
