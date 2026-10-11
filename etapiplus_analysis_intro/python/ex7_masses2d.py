"""Exercise 7: Two-dimensional mass plots (and a Dalitz plot).

After the standard cuts (Chi2_NDF < 10, -t < 1, |MM2| < 0.05), weighted by
accidweight:
    pad 1: M(eta pi+) vs M(eta n)
    pad 2: M(eta pi+) vs M(pi+ n)
    pad 3: M(eta n)   vs M(pi+ n)
    pad 4: Dalitz plot M^2(eta pi+) vs M^2(pi+ n)
and M(eta pi+) with and without a baryon veto (M(pi+ n) > 2 and M(eta n) > 2 GeV).

Mesons decaying to eta pi+ make HORIZONTAL bands (M(eta pi+) is on y);
baryons decaying to pi+ n or eta n make VERTICAL bands.

Run from the etapiplus_analysis_intro/ folder:
    python python/ex7_masses2d.py
Output: sums in the terminal, plots/py_ex7_masses2d.png and plots/py_ex7_veto.png
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
                 "Neutron_Px", "Neutron_Py", "Neutron_Pz", "Neutron_E",
                 "Chi2_NDF", "Mandelstam_t", "MissingMassSquared", "accidweight"],
                library="np")

# Standard cuts, then keep only the entries that pass
keep = ((d["Chi2_NDF"] < 10) & (-d["Mandelstam_t"] < 1.0)
        & (np.abs(d["MissingMassSquared"]) < 0.05))
d = {k: v[keep] for k, v in d.items()}
w = d["accidweight"]


def p4(name):
    """(E, px, py, pz) of one particle as four numpy arrays."""
    return np.array([d[name + "_E"], d[name + "_Px"], d[name + "_Py"], d[name + "_Pz"]])


def mass2(a, b):
    """Invariant mass squared of two four-vectors: (E1+E2)^2 - |p1+p2|^2."""
    s = a + b
    return s[0] ** 2 - s[1] ** 2 - s[2] ** 2 - s[3] ** 2


pip, eta, neu = p4("PiPlus"), p4("Eta"), p4("Neutron")
m2_etapi, m2_etan, m2_pin = mass2(eta, pip), mass2(eta, neu), mass2(pip, neu)
m_etapi, m_etan, m_pin = np.sqrt(m2_etapi), np.sqrt(m2_etan), np.sqrt(m2_pin)

panels = [  # (x, y, x bins, y bins, x label, y label, title)
    (m_etan, m_etapi, np.linspace(1.4, 4.4, 101), np.linspace(0.7, 2.5, 91),
     r"M($\eta n$) (GeV/$c^2$)", r"M($\eta\pi^+$) (GeV/$c^2$)", r"M($\eta\pi^+$) vs M($\eta n$)"),
    (m_pin, m_etapi, np.linspace(1.0, 4.0, 101), np.linspace(0.7, 2.5, 91),
     r"M($\pi^+n$) (GeV/$c^2$)", r"M($\eta\pi^+$) (GeV/$c^2$)", r"M($\eta\pi^+$) vs M($\pi^+n$)"),
    (m_pin, m_etan, np.linspace(1.0, 4.0, 101), np.linspace(1.4, 4.4, 101),
     r"M($\pi^+n$) (GeV/$c^2$)", r"M($\eta n$) (GeV/$c^2$)", r"M($\eta n$) vs M($\pi^+n$)"),
    (m2_pin, m2_etapi, np.linspace(1.0, 16.0, 101), np.linspace(0.0, 9.0, 101),
     r"M$^2$($\pi^+n$) (GeV$^2$/$c^4$)", r"M$^2$($\eta\pi^+$) (GeV$^2$/$c^4$)", "Dalitz plot"),
]

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, (x, y, xe, ye, xl, yl, title) in zip(axes.flat, panels):
    h, _, _ = np.histogram2d(x, y, bins=[xe, ye], weights=w)
    h = np.where(h < 0.1, np.nan, h)    # leave empty / negative bins white, like SetMinimum(0.1)
    mesh = ax.pcolormesh(xe, ye, h.T, cmap="viridis")
    fig.colorbar(mesh, ax=ax)
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title)
    if title == "Dalitz plot":
        print(f"Weighted entries in the Dalitz plot:          {np.nansum(np.histogram2d(x, y, bins=[xe, ye], weights=w)[0]):.1f}")
fig.tight_layout()
fig.savefig("plots/py_ex7_masses2d.png", dpi=100)

# Baryon veto
edges = np.linspace(0.7, 2.5, 91)
veto = (m_pin > 2.0) & (m_etan > 2.0)
h_all, _ = np.histogram(m_etapi, bins=edges, weights=w)
h_veto, _ = np.histogram(m_etapi[veto], bins=edges, weights=w[veto])
print(f"M(eta pi+) 0.7-2.5 GeV, standard cuts:        {h_all.sum():.1f}")
print(f"M(eta pi+) 0.7-2.5 GeV, plus baryon veto:     {h_veto.sum():.1f}")
print(f"Fraction kept by the veto:                    {h_veto.sum() / h_all.sum():.4f}")

fig, ax = plt.subplots(figsize=(8, 6))
ax.stairs(h_all, edges, color="black", lw=1.5, label="Standard cuts")
ax.stairs(h_veto, edges, color="tab:red", lw=1.5, label=r"+ M($\pi^+n$) > 2 and M($\eta n$) > 2 GeV")
ax.set_title("Baryon veto")
ax.set_xlabel(r"M($\eta\pi^+$) (GeV/$c^2$)")
ax.set_ylabel("Weighted entries / 20 MeV")
ax.set_ylim(0, 1.3 * h_all.max())
ax.legend(loc="upper left")
fig.tight_layout()
fig.savefig("plots/py_ex7_veto.png", dpi=100)
print("Saved plots/py_ex7_masses2d.png and plots/py_ex7_veto.png")
