#!/usr/bin/env python3
"""
Standalone replacement for etapi_delta_plotter: draw data vs. fit for
AmpTools eta' pi- Delta++ fits (reaction EtaPi Beam Proton EtaPrime Pi- Pi+)
that use the Zlm amplitude, optionally times BreitWigner factors (e.g. the
p pi+ Delta++ line shape).  Needs only Python with numpy, uproot and
matplotlib -- no AmpTools or ROOT build.

Final-state order: 0 = Proton, 1 = EtaPrime, 2 = Pi-, 3 = Pi+.  Plotted:
M(eta' pi-), M(p pi+), the GJ and helicity angles, -t (with the p pi+
system as the recoil, as in EtaPiDeltaPlotGenerator), and Phi.

Everything is read from the .fit file: the amplitudes, the fitted production
parameters and the data / accepted-MC / generated-MC file names (taken
relative to the .fit file's folder, as the fit ran there).  Each accepted-MC
event is weighted by the fitted intensity I(x) / N_generated, which is how
AmpTools predicts the yield, so the MC histograms show what the fit expects.

Usage:
  etapi_fit_plotter.py bin_10.fit                      # bin_10_fitplots.root + .pdf
  etapi_fit_plotter.py bin_10.fit -o out.root -t "1.0 < M < 1.1"
  etapi_fit_plotter.py EtaPi_fit/bin_*/bin_*.fit -p all_bins.pdf   # one page per fit
  etapi_fit_plotter.py bin_10.fit --check              # compare with the fit's own
                                                       # normalization integrals

Options kept from the C++ plotter: -o (ROOT file), -p (PDF), -t (title),
-s (total only, no per-sum histograms; the default here), -G (skip generated MC).
"""

import argparse
import math
import os
import re
import sys
from pathlib import Path

import numpy as np

TARGET_MASS = 0.938272       # as in EtaPiDeltaPlotGenerator
# final-state indices (0 = Proton, 1 = EtaPrime, 2 = Pi-, 3 = Pi+)
RES_IDX = (1, 2)             # the eta' pi- system the Zlm waves describe
ZLM_RECOIL_IDX = (0,)        # what Zlm (halld_sim) uses as the recoil: pKin[1] = Proton
# EtaPiDeltaPlotGenerator's recoil for t and the GJ frame is every other
# particle, i.e. p pi+ (indices 0 and 3)
# name: (x-axis label, nbins, low, high).  EtaPiMass, DeltaPPMass, cosTheta,
# Phi and t match the names and binning of EtaPiDeltaPlotGenerator
# (cosTheta/Phi in the GJ frame); cosThetaHel/PhiHel/BigPhi are the
# helicity-frame angles the Zlm amplitudes are built from.
VARIABLES = {
    "EtaPiMass":   (r"$M(\eta'\pi^-)$ [GeV]", 80, 0.4, 3.0),
    "DeltaPPMass": (r"$M(p\pi^+)$ [GeV]", 60, 0.8, 2.0),
    "cosTheta":    (r"$\cos\theta_{GJ}$", 50, -1.0, 1.0),
    "Phi":         (r"$\phi_{GJ}$ [deg]", 50, -180.0, 180.0),
    "t":           (r"$-t$ (to $p\pi^+$) [GeV$^2$]", 100, 0.0, 2.0),
    "cosThetaHel": (r"$\cos\theta_{hel}$", 50, -1.0, 1.0),
    "PhiHel":      (r"$\phi_{hel}$ [deg]", 50, -180.0, 180.0),
    "BigPhi":      (r"$\Phi$ (polarization) [deg]", 50, -180.0, 180.0),
}


# --------------------------------------------------------------------------
# .fit file
# --------------------------------------------------------------------------

def parse_fit(path):
    """Read the parts of an AmpTools .fit file needed here."""
    lines = Path(path).read_text().splitlines()
    sec = {}
    for i, l in enumerate(lines):
        if l.startswith("+++ "):
            sec.setdefault(l.strip("+ ").strip(), i)

    fit = {"path": Path(path)}

    # amplitudes and their scale values
    i = sec["Reactions, Amplitudes, and Scale Parameters"] + 1
    n_reac = int(lines[i]); i += 1
    fit["scales"] = {}
    for _ in range(n_reac):
        reac, n_amp = lines[i].split(); i += 1
        for _ in range(int(n_amp)):
            name, _par, val = lines[i].split(); i += 1
            fit["scales"][name] = float(val)
        fit.setdefault("reaction", reac)

    i = sec["Likelihood Total and Partial Sums"] + 1
    fit["likelihood"] = float(lines[i])

    i = sec["Parameter Values and Errors"] + 1
    n_par = int(lines[i]); i += 1
    fit["params"] = {}
    for _ in range(n_par):
        name, val = lines[i].split(); i += 1
        fit["params"][name] = float(val)

    # normalization integrals (first reaction)
    i = sec["Normalization Integrals"] + 1
    i += 1                                             # reaction name
    n_gen, n_acc = (float(x) for x in lines[i].split()); i += 1
    n_amp = int(lines[i]); i += 1
    ni_amps = [lines[i + k].strip() for k in range(n_amp)]; i += n_amp

    def matrix(i0):
        m = np.zeros((n_amp, n_amp), complex)
        for r in range(n_amp):
            vals = re.findall(r"\(([^,]+),([^)]+)\)", lines[i0 + r])
            m[r] = [complex(float(a), float(b)) for a, b in vals]
        return m
    fit["ni"] = {"n_gen": n_gen, "n_acc": n_acc, "amps": ni_amps,
                 "gen": matrix(i), "acc": matrix(i + n_amp)}

    # the equivalent config file at the end
    fit["amplitudes"], fit["files"] = {}, {}
    start = next(k for k, l in enumerate(lines) if l.startswith("### FIT CONFIGURATION"))
    for l in lines[start + 1:]:
        tok = l.split()
        if not tok or tok[0].startswith("#"):
            continue
        if tok[0] == "amplitude":
            # amplitude <reac::sum::amp> <class> args...  (several lines = product)
            fit["amplitudes"].setdefault(tok[1], []).append((tok[2], tok[3:]))
        elif tok[0] in ("data", "genmc", "accmc", "bkgnd"):
            fit["files"][tok[0]] = (tok[2], tok[3:])
        elif tok[0] == "fit":
            fit["name"] = tok[1]
    return fit


def production(fit, amp):
    p = fit["params"]
    return complex(p.get(amp + "_re", 0.0), p.get(amp + "_im", 0.0)) * fit["scales"].get(amp, 1.0)


# --------------------------------------------------------------------------
# kinematics (ROOTDataReader trees) and the Zlm amplitude
# --------------------------------------------------------------------------

def read_tree(path, tree="kin"):
    import uproot
    import awkward as ak
    with uproot.open(path) as f:
        t = f[tree]
        names = set(t.keys())
        a = t.arrays([k for k in ("E_Beam", "Px_Beam", "Py_Beam", "Pz_Beam",
                                  "E_FinalState", "Px_FinalState", "Py_FinalState",
                                  "Pz_FinalState", "Weight") if k in names])
    beam = np.stack([ak.to_numpy(a[k]) for k in ("Px_Beam", "Py_Beam", "Pz_Beam", "E_Beam")], -1)
    fs = np.stack([ak.to_numpy(a[k]) for k in ("Px_FinalState", "Py_FinalState",
                                             "Pz_FinalState", "E_FinalState")], -1)
    w = ak.to_numpy(a["Weight"]) if "Weight" in a.fields else np.ones(len(beam))
    return beam.astype(float), fs.astype(float), w.astype(float)


def mass(p):
    return np.sqrt(np.maximum(p[..., 3] ** 2 - np.sum(p[..., :3] ** 2, -1), 0.0))


def unit(v):
    return v / np.linalg.norm(v, axis=-1, keepdims=True)


def boost(p, b):
    """Boost four-vectors p (px,py,pz,E) by velocity b (like TLorentzVector::Boost)."""
    b2 = np.sum(b * b, -1, keepdims=True)
    gamma = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * p[..., :3], -1, keepdims=True)
    gamma2 = np.where(b2 > 0, (gamma - 1.0) / np.where(b2 > 0, b2, 1), 0.0)
    xyz = p[..., :3] + gamma2 * bp * b + gamma * b * p[..., 3:4]
    e = gamma * (p[..., 3:4] + bp)
    return np.concatenate([xyz, e], -1)


def kinematics(beam, fs, pol_angle_deg):
    """Same frame and angles as Zlm::calcUserVars (recoil = pKin[1], the proton;
    resonance = pKin[2] + pKin[3], EtaPrime + Pi-)."""
    recoil = fs[:, list(ZLM_RECOIL_IDX)].sum(1)
    p1, p2 = fs[:, RES_IDX[0]], fs[:, RES_IDX[1]]
    res = p1 + p2
    b = -res[:, :3] / res[:, 3:4]
    recoil_res, p1_res = boost(recoil, b), boost(p1, b)

    z = -unit(recoil_res[:, :3])                                   # helicity frame
    y = unit(np.cross(unit(beam[:, :3]), -unit(recoil[:, :3])))
    x = np.cross(y, z)
    v = p1_res[:, :3]
    vx, vy, vz = (np.sum(v * a, -1) for a in (x, y, z))
    cos_theta = vz / np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
    phi = np.arctan2(vy, vx)

    a = math.radians(pol_angle_deg)
    eps = np.array([math.cos(a), math.sin(a), 0.0])
    big_phi = np.arctan2(y @ eps, np.sum(unit(beam[:, :3]) * np.cross(eps, y), -1))

    return {"cosTheta": cos_theta, "Phi": phi, "bigPhi": big_phi}


def plot_variables(beam, fs, pol_angle_deg=0.0):
    """The EtaPiDeltaPlotGenerator variables, with p pi+ as the recoil."""
    others = [k for k in range(fs.shape[1]) if k not in RES_IDX]
    recoil = fs[:, others].sum(1)
    p1, p2 = fs[:, RES_IDX[0]], fs[:, RES_IDX[1]]
    res = p1 + p2

    cm = recoil + res
    b_cm = -cm[:, :3] / cm[:, 3:4]
    beam_cm, res_cm, eta_cm = boost(beam, b_cm), boost(res, b_cm), boost(p1, b_cm)
    b_gj = -res_cm[:, :3] / res_cm[:, 3:4]
    beam_gj, eta_gj = boost(beam_cm, b_gj), boost(eta_cm, b_gj)

    z = unit(beam_gj[:, :3])
    y = unit(np.cross(beam_cm[:, :3], res_cm[:, :3]))
    x = np.cross(y, z)
    v = eta_gj[:, :3]
    vx, vy, vz = (np.sum(v * a, -1) for a in (x, y, z))

    q = recoil - np.array([0.0, 0.0, 0.0, TARGET_MASS])
    hel = kinematics(beam, fs, pol_angle_deg)
    return {
        "EtaPiMass": mass(res),
        "DeltaPPMass": mass(recoil),
        "cosTheta": vz / np.sqrt(vx ** 2 + vy ** 2 + vz ** 2),
        "Phi": np.degrees(np.arctan2(vy, vx)),
        "t": -(q[:, 3] ** 2 - np.sum(q[:, :3] ** 2, -1)),
        "cosThetaHel": hel["cosTheta"],
        "PhiHel": np.degrees(hel["Phi"]),
        "BigPhi": np.degrees(hel["bigPhi"]),
    }


def ylm(l, m, cos_theta, phi):
    """Spherical harmonic Y_l^m with the Condon-Shortley phase."""
    am = abs(m)
    x = np.clip(cos_theta, -1.0, 1.0)
    pmm = np.ones_like(x)
    if am > 0:
        pmm = (-1) ** am * math.prod(range(1, 2 * am, 2)) * (1 - x * x) ** (am / 2)
    if l == am:
        plm = pmm
    else:
        pm1 = x * (2 * am + 1) * pmm
        plm, prev = pm1, pmm
        for ll in range(am + 2, l + 1):
            plm, prev = ((2 * ll - 1) * x * plm - (ll + am - 1) * prev) / (ll - am), plm
    norm = math.sqrt((2 * l + 1) / (4 * math.pi) * math.factorial(l - am) / math.factorial(l + am))
    y = norm * plm * np.exp(1j * am * phi)
    return y if m >= 0 else (-1) ** am * np.conj(y)


def to_float(s):
    """Zlm reads its numeric arguments with atof (a bare word gives 0)."""
    try:
        return float(s)
    except ValueError:
        return 0.0


def zlm_amplitude(args, kin_cache, beam, fs):
    if len(args) != 6:
        sys.exit(f"Zlm with {len(args)} arguments is not supported here (need 6: "
                 "J m r s polAngle polFraction)")
    j, m, r, s = (int(a) for a in args[:4])
    pol_angle, pol_frac = to_float(args[4]), to_float(args[5])
    if pol_angle not in kin_cache:
        kin_cache[pol_angle] = kinematics(beam, fs, pol_angle)
    k = kin_cache[pol_angle]
    rotated = ylm(j, m, k["cosTheta"], k["Phi"]) * np.exp(-1j * k["bigPhi"])
    part = rotated.real if r == 1 else rotated.imag
    return math.sqrt(1 + s * pol_frac) * part


def breakup_momentum(m0, m1, m2):
    """As halld_sim breakupMomentum (sqrt of the absolute value, also below threshold)."""
    q2 = (m0 ** 4 + m1 ** 4 + m2 ** 4 - 2 * m0 ** 2 * m1 ** 2
          - 2 * m0 ** 2 * m2 ** 2 - 2 * m1 ** 2 * m2 ** 2)
    return np.sqrt(np.abs(q2)) / (2.0 * m0)


def barrier_factor(q, spin):
    """Blatt-Weisskopf barrier factor as in halld_sim barrierFactor.cc."""
    z = q * q / (0.1973 * 0.1973)
    if spin == 0:
        return np.ones_like(z)
    if spin == 1:
        return np.sqrt(2.0 * z / (z + 1.0))
    if spin == 2:
        return np.sqrt(13.0 * z * z / ((z - 3.0) ** 2 + 9.0 * z))
    if spin == 3:
        return np.sqrt(277.0 * z ** 3 / (z * (z - 15.0) ** 2 + 9.0 * (2.0 * z - 5.0) ** 2))
    if spin == 4:
        return np.sqrt(12746.0 * z ** 4 / ((z * z - 45.0 * z + 105.0) ** 2
                                           + 25.0 * z * (2.0 * z - 21.0) ** 2))
    sys.exit(f"BreitWigner with spin {spin} is not supported here")


def breit_wigner_amplitude(args, kin_cache, beam, fs, params=None):
    """halld_sim BreitWigner: mass width spin daughters1 daughters2, where each
    daughters string lists pKin indices (0 = beam, so 1 = Proton, 4 = Pi+)."""
    if len(args) != 5:
        sys.exit(f"BreitWigner with {len(args)} arguments is not supported here "
                 "(need 5: mass width spin daughters1 daughters2)")

    def value(a):
        if a.startswith("[") and a.endswith("]"):
            return (params or {})[a[1:-1]]
        return float(a)
    m0, g0, spin = value(args[0]), value(args[1]), int(args[2])
    d1 = sum(fs[:, int(c) - 1] for c in args[3])
    d2 = sum(fs[:, int(c) - 1] for c in args[4])
    m, m1, m2 = mass(d1 + d2), mass(d1), mass(d2)
    q0 = np.abs(breakup_momentum(m0, m1, m2))
    q = np.abs(breakup_momentum(m, m1, m2))
    f0, f = barrier_factor(q0, spin), barrier_factor(q, spin)
    width = g0 * (m0 / m) * (q / q0) * (f * f) / (f0 * f0)
    return f * math.sqrt(m0 * g0 / 3.1416) / ((m0 * m0 - m * m) - 1j * (m0 * width))


AMPLITUDE_CLASSES = {"Zlm": zlm_amplitude, "BreitWigner": breit_wigner_amplitude}


def amplitudes(fit, beam, fs):
    """Complex amplitude per event for every amplitude in the fit."""
    kin_cache, out = {}, {}
    for name, factors in fit["amplitudes"].items():
        a = np.ones(len(beam), complex)
        for cls, args in factors:
            if cls not in AMPLITUDE_CLASSES:
                sys.exit(f"amplitude class {cls} is not implemented in this script")
            if cls == "BreitWigner":
                a = a * breit_wigner_amplitude(args, kin_cache, beam, fs, fit.get("params"))
            else:
                a = a * AMPLITUDE_CLASSES[cls](args, kin_cache, beam, fs)
        out[name] = a
    pol = next(iter(kin_cache)) if kin_cache else 0.0
    return out, plot_variables(beam, fs, pol)


def sums_of(fit):
    sums = {}
    for name in fit["amplitudes"]:
        reac, s, _ = name.split("::")
        sums.setdefault(s, []).append(name)
    return sums


def intensity(fit, amps, only_sum=None):
    total = 0.0
    for s, names in sums_of(fit).items():
        if only_sum and s != only_sum:
            continue
        coherent = sum(production(fit, n) * amps[n] for n in names)
        total = total + np.abs(coherent) ** 2
    return total


# --------------------------------------------------------------------------
# histograms, ROOT and PDF output
# --------------------------------------------------------------------------

def hist(values, weights, edges):
    h, _ = np.histogram(values, edges, weights=weights)
    e2, _ = np.histogram(values, edges, weights=weights ** 2)
    return h, np.sqrt(e2)


def make_th1(name, counts, errors, edges):
    import uproot
    under_over = lambda a: np.concatenate([[0.0], a, [0.0]])
    centers = 0.5 * (edges[1:] + edges[:-1])
    return uproot.writing.identify.to_TH1x(
        fName=name, fTitle=name, data=under_over(counts),
        fEntries=float(np.sum(counts)), fTsumw=float(np.sum(counts)),
        fTsumw2=float(np.sum(errors ** 2)), fTsumwx=float(np.sum(counts * centers)),
        fTsumwx2=float(np.sum(counts * centers ** 2)), fSumw2=under_over(errors ** 2),
        fXaxis=uproot.writing.identify.to_TAxis(
            fName="xaxis", fTitle="", fNbins=len(counts),
            fXmin=float(edges[0]), fXmax=float(edges[-1]), fXbins=edges.astype(float)))


def process(fit_path, args):
    fit = parse_fit(fit_path)
    base = Path(fit_path).resolve().parent
    tree = lambda key: fit["files"][key][1][0] if len(fit["files"][key][1]) > 1 else "kin"
    path = lambda key: str(base / fit["files"][key][1][0]) if fit["files"][key][1] else None

    out_hists, edges = {}, {}
    beam, fs, w = read_tree(path("data"), "kin")
    _, kin_data = amplitudes({"amplitudes": {}}, beam, fs)
    for var, (_, nb, lo, hi) in VARIABLES.items():
        x = kin_data[var]
        if args.zoom and var in ("EtaPiMass", "t"):
            # zoom to the data in this (narrow) bin
            pad = 0.02 * (x.max() - x.min() or 1.0)
            lo, hi, nb = max(lo, x.min() - pad), min(hi, x.max() + pad), 40
        edges[var] = np.linspace(lo, hi, nb + 1)
        out_hists[var + "dat"] = hist(x, w, edges[var])

    n_gen = fit["ni"]["n_gen"]
    kinds = [("accmc", "acc")] + ([("genmc", "gen")] if args.with_gen else [])
    check = None
    for key, suffix in kinds:
        beam, fs, w = read_tree(path(key), "kin")
        amps, kin = amplitudes(fit, beam, fs)
        weight = w * intensity(fit, amps) / n_gen
        for var in VARIABLES:
            out_hists[var + suffix] = hist(kin[var], weight, edges[var])
        if args.per_sum:
            for s in sums_of(fit):
                ws = w * intensity(fit, amps, s) / n_gen
                for var in VARIABLES:
                    out_hists[f"{var}{suffix}_{s}"] = hist(kin[var], ws, edges[var])
        if key == "accmc":
            check = (amps, w)

    if args.check and check:
        report_check(fit, *check)

    # ROOT file
    import uproot
    root_out = args.output or str(Path(fit_path).with_name(Path(fit_path).stem + "_fitplots.root"))
    with uproot.recreate(root_out) as f:
        for name, (c, e) in out_hists.items():
            var = next(v for v in VARIABLES if name.startswith(v))
            f[name] = make_th1(name, c, e, edges[var])

    title = args.title or Path(fit_path).stem
    return fit, out_hists, edges, title, root_out


def report_check(fit, amps, w):
    """Recompute the accepted-MC normalization integrals and compare with the fit's."""
    ni = fit["ni"]
    names = [n for n in ni["amps"] if n in amps]
    idx = [ni["amps"].index(n) for n in names]
    A = np.array([amps[n] for n in names])
    mine = (A * w) @ np.conj(A).T / ni["n_gen"]
    ref = ni["acc"][np.ix_(idx, idx)]
    # AmpTools only stores integrals between amplitudes in the same coherent sum
    sum_of = [n.split("::")[1] for n in names]
    same = np.array([[a == b for b in sum_of] for a in sum_of])
    scale = np.max(np.abs(np.diag(ref)))
    worst = np.max(np.abs(np.where(same, mine - ref, 0))) / scale
    n_pred = sum(production(fit, a) * np.conj(production(fit, b)) * ref[i, k]
                 for i, a in enumerate(names) for k, b in enumerate(names)).real
    print(f"  check: largest difference from the fit's normalization integrals = "
          f"{worst:.2e} of the largest diagonal term "
          f"({'OK' if worst < 1e-4 else 'MISMATCH - amplitudes or files differ from the fit'})")
    print(f"  check: predicted accepted yield from the fit = {n_pred:.1f}")


def draw(page_pdf, out_hists, edges, title, likelihood):
    import matplotlib.pyplot as plt
    ncols = 3
    nrows = math.ceil(len(VARIABLES) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4 * nrows))
    for ax, var in zip(axes.flat, VARIABLES):
        e = edges[var]
        c = 0.5 * (e[1:] + e[:-1])
        acc, _ = out_hists[var + "acc"]
        ax.stairs(acc, e, fill=True, color="#a8c5e8", label="fit (acc. MC)")
        d, de = out_hists[var + "dat"]
        ax.errorbar(c, d, yerr=de, fmt="o", ms=2.5, color="k", label="data")
        ax.set_xlabel(VARIABLES[var][0])
        ax.set_ylabel("events")
        ax.set_ylim(bottom=0)
    for ax in list(axes.flat)[len(VARIABLES):]:
        ax.axis("off")
    axes.flat[0].legend(loc="lower center", frameon=False, fontsize=9)
    fig.suptitle(f"{title}    (likelihood = {likelihood:.1f})")
    fig.tight_layout()
    page_pdf.savefig(fig)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fits", nargs="+", help=".fit file(s)")
    p.add_argument("-o", "--output", help="ROOT output (one .fit only; default <fit>_fitplots.root)")
    p.add_argument("-p", "--pdf", help="PDF output; with several .fit files, one page each")
    p.add_argument("-t", "--title", help="title on the plot (default: the fit name)")
    p.add_argument("-s", "--total-only", dest="per_sum", action="store_false",
                   help="total fit only (the default; kept for compatibility)")
    p.add_argument("--per-sum", dest="per_sum", action="store_true",
                   help="also write histograms for each coherent sum")
    p.add_argument("-G", "--no-gen", dest="with_gen", action="store_false",
                   help="skip generated MC (the default; kept for compatibility)")
    p.add_argument("--with-gen", dest="with_gen", action="store_true",
                   help="also histogram generated MC (acceptance corrected)")
    p.add_argument("--full-range", dest="zoom", action="store_false",
                   help="use the plot generator's full mass and t ranges and binning "
                        "instead of 40 bins over the range of the data in this bin")
    p.add_argument("--check", action="store_true",
                   help="compare with the fit's normalization integrals")
    p.set_defaults(per_sum=False, with_gen=False)
    args = p.parse_args()
    if args.output and len(args.fits) > 1:
        sys.exit("-o only works with a single .fit file")

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages

    pdf_out = args.pdf or (Path(args.output).with_suffix(".pdf") if args.output
                           else Path(args.fits[0]).with_name(Path(args.fits[0]).stem + "_fitplots.pdf"))
    with PdfPages(pdf_out) as pdf:
        for fp in args.fits:
            print(f"{fp}:")
            fit, out_hists, edges, title, root_out = process(fp, args)
            draw(pdf, out_hists, edges, title, fit["likelihood"])
            print(f"  wrote {root_out}")
    print(f"Wrote {pdf_out}")


if __name__ == "__main__":
    main()
