#!/usr/bin/env python3
"""
Fit-quality plots for a TwoPiAngles SDME fit (amptools_sdme.py, or any
AmpTools .fit file whose amplitudes are TwoPiAngles with a fixed polarization):
data compared with the accepted MC weighted by the fitted intensity, as in
Fig. 5 of arXiv:2305.09047.

For every reaction (polarization orientation) in the fit, and for their sum
when there are several, one page with
    cos(theta_H), phi_H, Phi, Phi - phi_H, M(pi- pi0), -t
where theta_H, phi_H are the pi- angles in the helicity frame and Phi is the
angle of the beam polarization vector (polAngle included) with respect to the
production plane, all computed exactly as TwoPiAngles does.  Data events carry
their Weight branch (accidental subtraction); the fit histogram is
    |V|^2 s^2 * 3/(4 pi) W(angles; SDMEs) * weight / N_gen
summed over the accepted MC, i.e. the number of events the fit predicts.

Usage (inside the folder of the fit, since the file paths in the .fit file are
relative to it):
    python3 sdme_fit_plotter.py all.fit [-o fitplots.pdf] [-t "title"] [--check]
--check compares the predicted yield with the fit's own normalization
integrals and prints OK or MISMATCH.
"""

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

SDME_ORDER = ["rho000", "rho100", "rho1m10", "rho111", "rho001",
              "rho101", "rho1m11", "rho102", "rho1m12"]


# --------------------------------------------------------------------------
# .fit file
# --------------------------------------------------------------------------

def parse_fit(path):
    lines = Path(path).read_text().splitlines()
    sec = {}
    for i, l in enumerate(lines):
        if l.startswith("+++ "):
            sec.setdefault(l.strip("+ ").strip(), i)
    fit = {"scale": {}, "params": {}, "ni": {}, "reactions": {}}

    i = sec["Reactions, Amplitudes, and Scale Parameters"] + 1
    n_reac = int(lines[i]); i += 1
    for _ in range(n_reac):
        reac, n_amp = lines[i].split(); i += 1
        for _ in range(int(n_amp)):
            amp, _par, val = lines[i].split(); i += 1
            fit["scale"][amp] = float(val)

    i = sec["Parameter Values and Errors"] + 1
    n_par = int(lines[i]); i += 1
    for _ in range(n_par):
        name, val = lines[i].split()[:2]; i += 1
        fit["params"][name] = float(val)

    # normalization integrals: per reaction  N_gen N_acc, n_amp, names, gen, acc
    i = sec["Normalization Integrals"] + 1
    while i < len(lines) and not lines[i].startswith("+++"):
        reac = lines[i].strip(); i += 1
        n_gen, n_acc = (float(x) for x in lines[i].split()); i += 1
        n_amp = int(lines[i]); i += 1
        amps = [lines[i + k].strip() for k in range(n_amp)]; i += n_amp
        acc = []
        for r in range(n_amp):
            vals = re.findall(r"\(([^,]+),([^)]+)\)", lines[i + n_amp + r])
            acc.append([complex(float(a), float(b)) for a, b in vals])
        i += 2 * n_amp
        fit["ni"][reac] = {"n_gen": n_gen, "n_acc": n_acc, "amps": amps, "acc": acc}

    start = next(k for k, l in enumerate(lines) if l.startswith("### FIT CONFIGURATION"))
    for l in lines[start + 1:]:
        tok = l.split()
        if not tok or tok[0].startswith("#"):
            continue
        if tok[0] == "reaction":
            fit["reactions"].setdefault(tok[1], {})
        elif tok[0] in ("data", "accmc", "genmc", "bkgnd"):
            # data <reaction> ROOTDataReader <file> [tree ...]
            fit["reactions"].setdefault(tok[1], {})[tok[0]] = tok[3]
        elif tok[0] == "amplitude" and tok[2] == "TwoPiAngles":
            reac = tok[1].split("::")[0]
            fit["reactions"].setdefault(reac, {}).setdefault("amps", []).append(
                (tok[1], tok[3:]))
    return fit


def value(arg, params):
    """An amplitude argument: [parameter] or a number."""
    if arg.startswith("[") and arg.endswith("]"):
        return params[arg[1:-1]]
    return float(arg)


def production(fit, amp):
    p = fit["params"]
    return complex(p.get(amp + "_re", 0.0), p.get(amp + "_im", 0.0))


# --------------------------------------------------------------------------
# kinematics (as TwoPiAngles::calcUserVars)
# --------------------------------------------------------------------------

def read_tree(path, tree="kin"):
    import uproot
    with uproot.open(path) as f:
        t = f[tree]
        keys = t.keys()
        a = t.arrays([k for k in ("E_Beam", "Px_Beam", "Py_Beam", "Pz_Beam",
                                  "E_FinalState", "Px_FinalState", "Py_FinalState",
                                  "Pz_FinalState", "Weight") if k in keys], library="np")
    n = len(a["E_Beam"])
    beam = np.stack([a["Px_Beam"], a["Py_Beam"], a["Pz_Beam"], a["E_Beam"]], axis=1).astype(float)

    def fs(k):
        if n == 0:
            return np.zeros((0, 4))
        return np.stack([np.stack([a["Px_FinalState"][e][k], a["Py_FinalState"][e][k],
                                   a["Pz_FinalState"][e][k], a["E_FinalState"][e][k]])
                         for e in range(n)]).astype(float)
    if n and isinstance(a["E_FinalState"][0], np.ndarray):
        # fixed-size arrays: vectorize
        E, X, Y, Z = (np.stack(a[k]) if a[k].dtype == object else a[k]
                      for k in ("E_FinalState", "Px_FinalState", "Py_FinalState", "Pz_FinalState"))
        parts = [np.stack([X[:, k], Y[:, k], Z[:, k], E[:, k]], axis=1).astype(float)
                 for k in range(3)]
    else:
        parts = [fs(k) for k in range(3)]
    w = a["Weight"].astype(float) if "Weight" in a else np.ones(n)
    return beam, parts, w


def boost(p, b):
    """Boost four-vectors p (N,4: px,py,pz,E) by velocity b (N,3)."""
    b2 = np.sum(b * b, axis=1)
    g = 1.0 / np.sqrt(1.0 - b2)
    bp = np.sum(b * p[:, :3], axis=1)
    g2 = np.where(b2 > 0, (g - 1.0) / np.where(b2 > 0, b2, 1), 0.0)
    v = p[:, :3] + (g2 * bp)[:, None] * b + (g * p[:, 3])[:, None] * b
    return np.column_stack([v, g * (p[:, 3] + bp)])


def unit(v):
    return v / np.linalg.norm(v, axis=1)[:, None]


def kinematics(beam, parts):
    recoil, p1, p2 = parts
    res = p1 + p2
    bvec = -res[:, :3] / res[:, 3:4]
    recoil_res = boost(recoil, bvec)
    p1_res = boost(p1, bvec)
    y = unit(np.cross(unit(beam[:, :3]), -unit(recoil[:, :3])))
    z = -unit(recoil_res[:, :3])
    x = unit(np.cross(y, z))
    px, py, pz = (np.sum(p1_res[:, :3] * a, axis=1) for a in (x, y, z))
    theta = np.arctan2(np.hypot(px, py), pz)
    phi = np.arctan2(py, px)
    eps = np.array([1.0, 0.0, 0.0])
    bigphi0 = np.arctan2(y @ eps, np.sum(unit(beam[:, :3]) * np.cross(eps, y), axis=1))
    m = np.sqrt(np.maximum(res[:, 3] ** 2 - np.sum(res[:, :3] ** 2, axis=1), 0))
    q = beam - res
    mt = -(q[:, 3] ** 2 - np.sum(q[:, :3] ** 2, axis=1))
    return {"theta": theta, "phi": phi, "bigphi0": bigphi0, "mass": m, "t": mt}


def twopi_W(k, rho, pol_angle_deg, pgamma):
    """3/(4 pi) W as in TwoPiAngles::calcAmplitude (|amplitude|^2)."""
    r = dict(zip(SDME_ORDER, rho))
    th, phi = k["theta"], k["phi"]
    big = np.radians(pol_angle_deg) + k["bigphi0"]
    ct, s2, s2t = np.cos(th), np.sin(th) ** 2, np.sin(2 * th)
    W = (0.5 * (1 - r["rho000"]) + 0.5 * (3 * r["rho000"] - 1) * ct ** 2
         - math.sqrt(2) * r["rho100"] * s2t * np.cos(phi) - r["rho1m10"] * s2 * np.cos(2 * phi))
    W -= pgamma * np.cos(2 * big) * (r["rho111"] * s2 + r["rho001"] * ct ** 2
                                     - math.sqrt(2) * r["rho101"] * s2t * np.cos(phi)
                                     - r["rho1m11"] * s2 * np.cos(2 * phi))
    W -= pgamma * np.sin(2 * big) * (math.sqrt(2) * r["rho102"] * s2t * np.sin(phi)
                                     + r["rho1m12"] * s2 * np.sin(2 * phi))
    return np.abs(W * 3.0 / (4.0 * math.pi))


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def variables(k, pol_angle_deg):
    big = wrap(np.radians(pol_angle_deg) + k["bigphi0"])
    return {
        "costheta": np.cos(k["theta"]),
        "phi": np.degrees(k["phi"]),
        "Phi": np.degrees(big),
        "Phi_minus_phi": np.degrees(wrap(big - k["phi"])),
        "mass": k["mass"],
        "t": k["t"],
    }


PANELS = [
    ("costheta", r"$\cos\vartheta_H$", 40, None),
    ("phi", r"$\varphi_H$ [deg]", 36, (-180, 180)),
    ("Phi", r"$\Phi$ [deg]", 36, (-180, 180)),
    ("Phi_minus_phi", r"$\Phi-\varphi_H$ [deg]", 36, (-180, 180)),
    ("mass", r"$M(\pi^-\pi^0)$ [GeV]", 35, None),
    ("t", r"$-t$ [GeV$^2$]", 25, None),
]


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def process(fit, tree):
    """Per reaction: data and fit variables with weights."""
    out = {}
    cache = {}

    def load(path):
        if path not in cache:
            beam, parts, w = read_tree(path, tree)
            cache[path] = (kinematics(beam, parts), w)
        return cache[path]

    for reac, info in fit["reactions"].items():
        if "data" not in info or "accmc" not in info or not info.get("amps"):
            continue
        kd, wd = load(info["data"])
        ka, wa = load(info["accmc"])
        ni = fit["ni"].get(reac)
        n_gen = ni["n_gen"] if ni else len(wa)
        pred = np.zeros(len(wa))
        pol_angle = None
        for amp, args in info["amps"]:
            rho = [value(a, fit["params"]) for a in args[:9]]
            pol_angle = value(args[9], fit["params"])
            if len(args) != 11:
                sys.exit(f"{amp}: only TwoPiAngles with a fixed polarization fraction is supported")
            pgamma = float(args[10])
            vs = abs(production(fit, amp) * fit["scale"].get(amp, 1.0)) ** 2
            pred += vs * twopi_W(ka, rho, pol_angle, pgamma)
        out[reac] = {
            "data": (variables(kd, pol_angle), wd),
            "fit": (variables(ka, pol_angle), pred * wa / n_gen),
            "pol_angle": pol_angle,
            "check": (pred * wa / n_gen).sum(),
            "check_ni": sum(abs(production(fit, amp) * fit["scale"].get(amp, 1.0)) ** 2
                            * ni["acc"][ni["amps"].index(amp)][ni["amps"].index(amp)].real
                            for amp, _ in info["amps"]) if ni else None,
        }
    return out


def hist(x, w, edges):
    h, _ = np.histogram(x, bins=edges, weights=w)
    e2, _ = np.histogram(x, bins=edges, weights=w * w)
    return h, np.sqrt(e2)


def draw_page(pdf, pages, title):
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for ax, (key, label, nb, rng) in zip(axes.flat, PANELS):
        allx = np.concatenate([p["data"][0][key] for p in pages] +
                              [p["fit"][0][key] for p in pages])
        lo, hi = rng if rng else ((allx.min(), allx.max()) if len(allx) else (0, 1))
        if key == "costheta":
            lo, hi = -1, 1
        edges = np.linspace(lo, hi, nb + 1)
        hd = np.zeros(nb); ed2 = np.zeros(nb); hf = np.zeros(nb)
        for p in pages:
            h, e = hist(p["data"][0][key], p["data"][1], edges)
            hd += h; ed2 += e * e
            hf += np.histogram(p["fit"][0][key], bins=edges, weights=p["fit"][1])[0]
        ctr = 0.5 * (edges[1:] + edges[:-1])
        ax.fill_between(edges, np.append(hf, hf[-1]), step="post", color="#8fd18f",
                        alpha=0.8, label="fit (weighted acc. MC)", lw=0)
        ax.errorbar(ctr, hd, yerr=np.sqrt(ed2), fmt="o", ms=2.5, color="black",
                    lw=0.8, label="data")
        ax.set_xlim(lo, hi)
        ax.set_ylim(bottom=min(0, (hd - np.sqrt(ed2)).min() * 1.1 if nb else 0))
        ax.set_xlabel(label)
        ax.set_ylabel(f"events / {(hi - lo) / nb:.3g}")
        ax.tick_params(direction="in", top=True, right=True)
    axes.flat[0].legend(frameon=False, fontsize=8, loc="lower center")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fit", help="AmpTools .fit file")
    ap.add_argument("-o", "--output", default="fitplots.pdf", help="output PDF")
    ap.add_argument("-t", "--title", default="", help="title on every page")
    ap.add_argument("--tree", default="kin")
    ap.add_argument("--check", action="store_true",
                    help="compare the predicted yield with the normalization integrals")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    from matplotlib.backends.backend_pdf import PdfPages

    fit = parse_fit(args.fit)
    res = process(fit, args.tree)
    if not res:
        sys.exit(f"{args.fit}: no reaction with data, accmc and TwoPiAngles amplitudes")

    if args.check:
        bad = False
        for reac, r in res.items():
            ok = r["check_ni"] is None or abs(r["check"] - r["check_ni"]) <= 1e-6 * max(1, abs(r["check_ni"]))
            bad |= not ok
            print(f"{reac}: predicted {r['check']:.6g}, from normalization integrals "
                  f"{r['check_ni']:.6g} {'OK' if ok else 'MISMATCH'}")
        print("OK" if not bad else "MISMATCH")

    with PdfPages(args.output) as pdf:
        if len(res) > 1:
            draw_page(pdf, list(res.values()), f"{args.title}  all orientations".strip())
        for reac, r in res.items():
            draw_page(pdf, [r], f"{args.title}  {reac} (pol = {r['pol_angle']:g} deg)".strip())
    for reac, r in res.items():
        nd = r["data"][1].sum()
        print(f"{reac}: data {nd:.1f} (weighted), fit {r['check']:.1f}")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
