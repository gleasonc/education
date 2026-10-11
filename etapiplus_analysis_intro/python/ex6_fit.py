"""Exercise 6: Fit the a0(980) peak.

Makes the eta pi+ mass histogram with the cuts from exercise 3 (weighted by
accidweight, 10 MeV bins) and fits the peak near 1 GeV with a Gaussian plus
a straight-line background.

Run from the tutorial/ folder:
    python python/ex6_fit.py
    python python/ex6_fit.py ../flat_pipneta__B4_v6.root 0.85 1.15   # other fit range
Output: fit results in the terminal and plots/py_ex6_fit.png
"""
import os
import sys
import numpy as np
import uproot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chi2 as chi2_dist

filename = sys.argv[1] if len(sys.argv) > 1 else "../flat_pipneta__B4_v6.root"
fit_lo = float(sys.argv[2]) if len(sys.argv) > 2 else 0.88
fit_hi = float(sys.argv[3]) if len(sys.argv) > 3 else 1.10
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
w = d["accidweight"]
keep = ((d["Chi2_NDF"] < 10) & (-d["Mandelstam_t"] < 1.0)
        & (np.abs(d["MissingMassSquared"]) < 0.05))

edges = np.linspace(0.7, 2.5, 181)                    # 180 bins of 10 MeV
counts, _ = np.histogram(mass[keep], bins=edges, weights=w[keep])
sumw2, _ = np.histogram(mass[keep], bins=edges, weights=w[keep] ** 2)
errors = np.sqrt(sumw2)                               # error on a weighted bin
centres = 0.5 * (edges[1:] + edges[:-1])


def model(x, amp, mean, sigma, b0, b1):
    """Gaussian peak plus a straight-line background."""
    return amp * np.exp(-0.5 * ((x - mean) / sigma) ** 2) + b0 + b1 * x


sel = (centres > fit_lo) & (centres < fit_hi)         # bins inside the fit range

# Starting values: height ~500 above background at ~0.99 GeV, ~30 MeV wide.
start = [500, 0.99, 0.03, 450, 0]
params, cov = curve_fit(model, centres[sel], counts[sel], p0=start,
                        sigma=errors[sel], absolute_sigma=True)
uncert = np.sqrt(np.diag(cov))
amp, mean, sigma, b0, b1 = params

chi2 = np.sum(((counts[sel] - model(centres[sel], *params)) / errors[sel]) ** 2)
ndf = sel.sum() - len(params)

print("===== Summary =====")
print(f"Mean     = {1000 * mean:.2f} +/- {1000 * uncert[1]:.2f} MeV")
print(f"Sigma    = {1000 * abs(sigma):.2f} +/- {1000 * uncert[2]:.2f} MeV")
print(f"chi2/ndf = {chi2:.2f} / {ndf}  (p-value {chi2_dist.sf(chi2, ndf):.4f})")
# Area of a Gaussian = amp * sigma * sqrt(2 pi); divide by bin width for a count
yield_ = amp * abs(sigma) * np.sqrt(2 * np.pi) / (edges[1] - edges[0])
print(f"Signal yield ~ {yield_:.0f} (weighted entries)")

fig, ax = plt.subplots(figsize=(8, 6))
show = (centres > 0.7) & (centres < 1.6)
ax.errorbar(centres[show], counts[show], yerr=errors[show], fmt="o", ms=3,
            color="black", label="Data")
x = np.linspace(fit_lo, fit_hi, 400)
ax.plot(x, model(x, *params), color="tab:red", lw=2, label="Gaussian + line")
ax.plot(x, b0 + b1 * x, color="gray", ls="--", label="Background only")
ax.set_title(r"a$_0$(980) $\rightarrow\eta\pi^+$ fit")
ax.set_xlabel(r"M($\eta\pi^+$) (GeV/$c^2$)")
ax.set_ylabel("Weighted entries / 10 MeV")
ax.set_ylim(0, 1.6 * counts[show].max())
ax.legend(loc="upper left")
ax.text(0.97, 0.97,
        f"mean = {1000 * mean:.1f} ± {1000 * uncert[1]:.1f} MeV\n"
        f"σ = {1000 * abs(sigma):.1f} ± {1000 * uncert[2]:.1f} MeV\n"
        f"χ²/ndf = {chi2:.1f} / {ndf}",
        transform=ax.transAxes, ha="right", va="top",
        bbox=dict(facecolor="white", edgecolor="gray"))
fig.tight_layout()
fig.savefig("plots/py_ex6_fit.png", dpi=100)
print("Saved plots/py_ex6_fit.png")
