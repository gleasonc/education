"""Exercise 8: Decay angles in the Gottfried-Jackson (GJ) frame.

Same steps as cpp/ex8_gj_angles.C (and the group's MakeFlat.C):
  1. boost everything from the lab to the CM frame of beam + target,
  2. boost from the CM to the eta pi+ rest frame,
  3. axes: z = beam direction in the eta pi+ rest frame,
           y = (beam x eta pi+) in the CM frame, x = y x z,
  4. cos(theta_GJ) and phi_GJ are the polar and azimuthal angles of the eta
     in those axes.
Cross-check: cos(theta_GJ) again, boosting straight from the lab to the
eta pi+ rest frame and taking the angle between eta and beam.

Run from the etapiplus_analysis_intro/ folder:
    python python/ex8_gj_angles.py
Output: checks and sums in the terminal and plots/py_ex8_gj_angles.png
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
                 "Beam_Px", "Beam_Py", "Beam_Pz", "Beam_E",
                 "Chi2_NDF", "Mandelstam_t", "MissingMassSquared", "accidweight"],
                library="np")
keep = ((d["Chi2_NDF"] < 10) & (-d["Mandelstam_t"] < 1.0)
        & (np.abs(d["MissingMassSquared"]) < 0.05))
d = {k: v[keep] for k, v in d.items()}
w = d["accidweight"]


def p4(name):
    """(E, px, py, pz) of one particle as a 4 x N array."""
    return np.array([d[name + "_E"], d[name + "_Px"], d[name + "_Py"], d[name + "_Pz"]])


def boost(v, b):
    """Boost four-vectors v (4 x N) by velocity b (3 x N), same formula as TLorentzVector::Boost."""
    E, p = v[0], v[1:]
    b2 = np.sum(b * b, axis=0)
    gamma = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * p, axis=0)
    gamma2 = np.where(b2 > 0, (gamma - 1.0) / np.where(b2 > 0, b2, 1.0), 0.0)
    p_new = p + gamma2 * bp * b + gamma * b * E
    E_new = gamma * (E + bp)
    return np.vstack([E_new, p_new])


def to_rest_frame_of(obj):
    """Velocity that takes you to the rest frame of obj: -p/E."""
    return -obj[1:] / obj[0]


def unit(v):
    return v / np.linalg.norm(v, axis=0)


beam, pip, eta = p4("Beam"), p4("PiPlus"), p4("Eta")
etapi = eta + pip
m_etapi = np.sqrt(etapi[0] ** 2 - np.sum(etapi[1:] ** 2, axis=0))
target = np.zeros_like(beam)
target[0] = 0.938272                      # proton at rest

# Step 1: lab -> CM
b_cm = to_rest_frame_of(beam + target)
beam_cm, eta_cm, etapi_cm = boost(beam, b_cm), boost(eta, b_cm), boost(etapi, b_cm)

# Step 2: CM -> eta pi+ rest frame
b_gj = to_rest_frame_of(etapi_cm)
beam_gj, eta_gj = boost(beam_cm, b_gj), boost(eta_cm, b_gj)

# Step 3: axes
z = unit(beam_gj[1:])
y = unit(np.cross(beam_cm[1:], etapi_cm[1:], axis=0))
x = np.cross(y, z, axis=0)

# Step 4: eta direction in (x, y, z)
vx, vy, vz = (np.sum(eta_gj[1:] * a, axis=0) for a in (x, y, z))
cos_theta = vz / np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
phi = np.degrees(np.arctan2(vy, vx))

# Cross-check: straight from the lab to the eta pi+ rest frame
b_direct = to_rest_frame_of(etapi)
cos_check = np.sum(unit(boost(eta, b_direct)[1:]) * unit(boost(beam, b_direct)[1:]), axis=0)

a0 = (m_etapi > 0.92) & (m_etapi < 1.04)
a2 = (m_etapi > 1.24) & (m_etapi < 1.40)
print(f"Entries used:                    {len(w)}")
print(f"cos(theta_GJ) range:             {cos_theta.min():.6f} to {cos_theta.max():.6f}")
print(f"phi_GJ range (deg):              {phi.min():.3f} to {phi.max():.3f}")
print(f"Max |cos - cross-check|:         {np.max(np.abs(cos_theta - cos_check)):.3g}")
in_range = (cos_theta >= -1) & (cos_theta < 1)
print(f"Weighted entries, all masses:    {w[in_range].sum():.1f}")
print(f"Weighted entries, a0 window:     {w[a0].sum():.2f}  (mean cos = {np.average(cos_theta[a0], weights=w[a0]):.4f})")
print(f"Weighted entries, a2 window:     {w[a2].sum():.2f}  (mean cos = {np.average(cos_theta[a2], weights=w[a2]):.4f})")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
ce = np.linspace(-1, 1, 51)
pe = np.linspace(-180, 180, 61)

ax = axes[0, 0]
ax.stairs(np.histogram(cos_theta, ce, weights=w)[0], ce, lw=1.5)
ax.set_title(r"cos$\theta_{GJ}$ (all masses)")
ax.set_xlabel(r"cos$\theta_{GJ}$")
ax.set_ylabel("Weighted entries")
ax.set_ylim(0, None)

ax = axes[0, 1]
ax.stairs(np.histogram(phi, pe, weights=w)[0], pe, lw=1.5)
ax.set_title(r"$\phi_{GJ}$ (all masses)")
ax.set_xlabel(r"$\phi_{GJ}$ (deg)")
ax.set_ylabel("Weighted entries")
ax.set_ylim(0, None)

ax = axes[1, 0]
me = np.linspace(0.7, 2.5, 46)
ce2 = np.linspace(-1, 1, 41)
h, _, _ = np.histogram2d(m_etapi, cos_theta, bins=[me, ce2], weights=w)
mesh = ax.pcolormesh(me, ce2, np.where(h < 0.1, np.nan, h).T, cmap="viridis")
fig.colorbar(mesh, ax=ax)
ax.set_title(r"cos$\theta_{GJ}$ vs M($\eta\pi^+$)")
ax.set_xlabel(r"M($\eta\pi^+$) (GeV/$c^2$)")
ax.set_ylabel(r"cos$\theta_{GJ}$")

ax = axes[1, 1]
cw = np.linspace(-1, 1, 26)
h0 = np.histogram(cos_theta[a0], cw, weights=w[a0])[0]
h2 = np.histogram(cos_theta[a2], cw, weights=w[a2])[0]
ax.stairs(h0 / h0.sum(), cw, color="tab:blue", lw=1.5, label="a$_0$(980): 0.92-1.04 GeV")
ax.stairs(h2 / h2.sum(), cw, color="tab:red", lw=1.5, label="a$_2$(1320): 1.24-1.40 GeV")
ax.set_title(r"cos$\theta_{GJ}$ in mass windows")
ax.set_xlabel(r"cos$\theta_{GJ}$")
ax.set_ylabel("Fraction / 0.08")
ax.set_ylim(0, 1.5 * max(h0.max() / h0.sum(), h2.max() / h2.sum()))
ax.legend(loc="upper left")

fig.tight_layout()
fig.savefig("plots/py_ex8_gj_angles.png", dpi=100)
print("Saved plots/py_ex8_gj_angles.png")
