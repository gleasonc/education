#!/usr/bin/env python3
"""
Automate binned AmpTools SDME fits of gamma p(n) -> rho- p (p missing),
rho- -> pi- pi0 (reaction Pi0Pim Beam Proton Pi- Pi0, TwoPiAngles), with
the four beam-polarization orientations 0/45/90/135:

  1. bin     split each sample (data, acc MC, gen MC) into M(pi- pi0) x -t bins
             with split_bins.C, and the samples listed in `split` (normally
             only DATAFILE) also by orientation, using the `pol` branch the v2
             DSelector writes.  Samples that are not split are used whole for
             every orientation.
  2. fit     in each bin, run
               all/      one fit of all four orientations with shared SDMEs
                         (fit_pi0pim_SDME_4pol_TEMPLATE.cfg, `cfg`)
               pol00/ pol45/ pol90/ pol135/
                         one fit per orientation
                         (fit_pi0pim_SDME_1pol_TEMPLATE.cfg, `cfg_single`)
             Leave `cfg` or `cfg_single` empty to skip that kind of fit.
  3. gather  read every .fit file: the nine SDMEs with errors, the
             natural / unnatural-parity components, the parity asymmetry
             P_sigma and the beam asymmetry Sigma (errors from the fit's
             covariance matrix), into results.csv (one row per bin and fit)
  4. plot    plots/sdmes.pdf (3x3 SDMEs vs -t, one colour per orientation and
             one for the combined fit), plots/parity_components.pdf
             (natural / unnatural, as Fig. 7 of arXiv:2305.09047) and
             plots/asymmetries.pdf (P_sigma, Sigma); the same with only the
             combined fit in plots/*_total.pdf
  5. plotfit data vs. the fit (accepted MC weighted by the fitted intensity) in
             cos(theta_H), phi_H, Phi, Phi - phi_H, M(pi- pi0) and -t, as Fig. 5
             of arXiv:2305.09047 (sdme_fit_plotter.py): <bin>/<fit>/fitplots.pdf,
             merged per fit into plots/fit_quality_<fit>.pdf.  `--no-plotfit`
             (or plotfit = false in the INI) leaves it out of `all`.

Natural / unnatural components (arXiv:2305.09047 Eq. 22,
rho^{N,U}_{ik} = 1/2 [rho^0_{ik} -+ (-1)^i rho^1_{-ik}]):
  rho^N,U_00   = 1/2 (rho^0_00 -+ rho^1_00)
  rho^N,U_11   = 1/2 (rho^0_11 +- rho^1_1-1),  rho^0_11 = (1 - rho^0_00)/2
  rho^N,U_1-1  = 1/2 (rho^0_1-1 +- rho^1_11)
  Re rho^N,U_10 = 1/2 (Re rho^0_10 -+ Re rho^1_10)
Parity asymmetry (Eq. 23):   P_sigma = 2 rho^1_1-1 - rho^1_00
Beam asymmetry (decay angles integrated, Schilling et al.):
                             Sigma   = 2 rho^1_11 + rho^1_00

Each step is a subcommand; `all` runs them in order.  Settings go in an INI
file (see example_sdme.ini) passed with -c; command-line options override it:

  amptools_sdme.py all -c my_sdme.ini
  amptools_sdme.py fit -c my_sdme.ini -r 5
"""

import argparse
import configparser
import csv
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

BINS_JSON = "bins.json"
RESULTS_CSV = "results.csv"
COMBINED = "all"

# TwoPiAngles parameters, in the order of the plots (as the reference plot)
SDMES = [
    ("rho000",  r"$\rho^0_{00}$"),
    ("rho100",  r"Re$\rho^0_{10}$"),
    ("rho1m10", r"Re$\rho^0_{1-1}$"),
    ("rho111",  r"Re$\rho^1_{11}$"),
    ("rho001",  r"Re$\rho^1_{00}$"),
    ("rho101",  r"Re$\rho^1_{10}$"),
    ("rho1m11", r"$\rho^1_{1-1}$"),
    ("rho102",  r"Im$\rho^2_{10}$"),
    ("rho1m12", r"Im$\rho^2_{1-1}$"),
]

# derived quantities: name -> (constant, {sdme: coefficient}, label)
DERIVED = {
    "rhoN_11":   (0.25, {"rho000": -0.25, "rho1m11": 0.5}, r"$\rho^N_{11}$"),
    "rhoN_00":   (0.0, {"rho000": 0.5, "rho001": -0.5}, r"$\rho^N_{00}$"),
    "rhoN_1m1":  (0.0, {"rho1m10": 0.5, "rho111": 0.5}, r"$\rho^N_{1-1}$"),
    "ReRhoN_10": (0.0, {"rho100": 0.5, "rho101": -0.5}, r"Re$\rho^N_{10}$"),
    "rhoU_11":   (0.25, {"rho000": -0.25, "rho1m11": -0.5}, r"$\rho^U_{11}$"),
    "rhoU_00":   (0.0, {"rho000": 0.5, "rho001": 0.5}, r"$\rho^U_{00}$"),
    "rhoU_1m1":  (0.0, {"rho1m10": 0.5, "rho111": -0.5}, r"$\rho^U_{1-1}$"),
    "ReRhoU_10": (0.0, {"rho100": 0.5, "rho101": 0.5}, r"Re$\rho^U_{10}$"),
    "P_sigma":   (0.0, {"rho1m11": 2.0, "rho001": -1.0}, r"$P_\sigma$"),
    "Sigma":     (0.0, {"rho111": 2.0, "rho001": 1.0}, r"$\Sigma$"),
}
# value for s-channel helicity conservation + natural-parity exchange
SCHC_NPE = {"rhoN_11": 0.5, "P_sigma": 1.0, "Sigma": 1.0}


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def run(cmd, cwd=None, log=None, dry=False):
    """Run a command (list or string). Returns the exit code."""
    shown = cmd if isinstance(cmd, str) else " ".join(shlex.quote(c) for c in cmd)
    print(f"[{cwd or '.'}] $ {shown}")
    if dry:
        return 0
    with open(log, "w") if log else open(os.devnull, "w") as out:
        proc = subprocess.run(cmd, cwd=cwd, shell=isinstance(cmd, str),
                              stdout=out if log else None,
                              stderr=subprocess.STDOUT if log else None)
    return proc.returncode


def uniform_edges(lo, hi, n):
    step = (hi - lo) / n
    return [lo + i * step for i in range(n + 1)]


def bin_name(i, j):
    return f"bin_m{i:03d}_t{j:02d}"


def main_fit_file(d, name):
    """With -r the per-attempt results are <name>_<k>.fit; prefer <name>.fit."""
    f = Path(d) / f"{name}.fit"
    if f.exists():
        return f
    fitfiles = sorted(Path(d).glob("*.fit"),
                      key=lambda p: (bool(re.search(r"_\d+$", p.stem)), p.name))
    return fitfiles[0] if fitfiles else None


def n_events(path, tree="kin"):
    """Number of entries in a ROOT file's tree, or None if it can't be read
    (no uproot, missing file)."""
    try:
        import uproot
        with uproot.open(path) as f:
            try:
                return f[tree].num_entries
            except KeyError:
                return 0
    except Exception:
        return None


def load_bins(workdir):
    path = Path(workdir) / BINS_JSON
    if not path.exists():
        sys.exit(f"{path} not found; run the 'bin' step first.")
    return json.loads(path.read_text())


def select_bins(bins, names):
    """Only the bins named in --bin (all bins without it)."""
    if not names:
        return bins
    unknown = [n for n in names if n not in {b["name"] for b in bins}]
    if unknown:
        sys.exit(f"--bin: no bin {' '.join(unknown)} in {BINS_JSON}")
    return [b for b in bins if b["name"] in names]


def parse_orientations(items):
    """['00:0', '45:45', ...] -> [('00', 0), ('45', 45), ...]"""
    out = []
    for it in items:
        if ":" not in it:
            sys.exit(f"orientations expects SUFFIX:POLVALUE, got '{it}'")
        suf, val = it.split(":", 1)
        out.append((suf, int(val)))
    return out


# --------------------------------------------------------------------------
# step 1: binning
# --------------------------------------------------------------------------

def step_bin(args):
    work = Path(args.workdir)
    tmp = work / "_split"
    tmp.mkdir(parents=True, exist_ok=True)

    mlo, mhi, nm = args.mass
    nm = int(nm)
    m_edges = uniform_edges(mlo, mhi, nm)
    if args.t:
        tlo, thi, nt = args.t
        nt = int(nt)
        if args.t_log:
            r = (thi / tlo) ** (1.0 / nt)
            t_edges = [tlo * r ** k for k in range(nt + 1)]
        else:
            t_edges = uniform_edges(tlo, thi, nt)
    else:
        nt, t_edges = 1, None

    samples = {}
    for item in (x for group in args.sample for x in group):
        if "=" not in item:
            sys.exit(f"--sample expects TOKEN=FILE, got '{item}'")
        token, path = item.split("=", 1)
        samples[token] = path
    orients = parse_orientations(args.orientations)
    unknown = [s for s in args.split if s not in samples]
    if unknown:
        sys.exit(f"split: {' '.join(unknown)} not among the samples ({' '.join(samples)})")

    # (token of the output file, input file, pol value or None)
    jobs = []
    for token, path in samples.items():
        if token in args.split:
            jobs += [(f"{token}_{suf}", path, val) for suf, val in orients]
        else:
            jobs.append((token, path, None))

    tree = args.tree or "kin"
    for out_tok, path, pol in jobs:
        src = Path(path).resolve()
        call = (f'{Path(args.splitter).resolve()}("{src}", "{tmp / out_tok}", '
                f'{mlo}, {mhi}, {nm}, '
                + (f'{tlo}, {thi}, {nt}, {"true" if args.t_log else "false"}, '
                   if t_edges else '0, 0, 0, false, ')
                + f'"{tree}", "{args.mass_indices}", '
                + (f'"{args.pol_branch}", {pol})' if pol is not None else '"", 0)'))
        rc = run([args.root_exe, "-l", "-b", "-q", call], dry=args.dry_run)
        if rc:
            sys.exit(f"split_bins.C failed on {src} (exit {rc})")
        for i in range(nm):
            for j in range(nt):
                part = tmp / f"{out_tok}_m{i}_t{j}.root"
                dest = work / bin_name(i, j) / f"{out_tok.lower()}.root"
                if not args.dry_run:
                    if not part.exists():
                        sys.exit(f"split_bins.C did not write {part}")
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(part), dest)

    # token -> file in the bin folder; every sample gets a TOKEN_<suf> entry
    files = {}
    for token in samples:
        for suf, _ in orients:
            name = f"{token}_{suf}" if token in args.split else token
            files[f"{token}_{suf}"] = f"{name.lower()}.root"
        if token not in args.split:
            files[token] = f"{token.lower()}.root"

    bins = []
    for i in range(nm):
        for j in range(nt):
            bins.append({
                "name": bin_name(i, j), "i_mass": i, "i_t": j,
                "mass_lo": m_edges[i], "mass_hi": m_edges[i + 1],
                "t_lo": t_edges[j] if t_edges else None,
                "t_hi": t_edges[j + 1] if t_edges else None,
                "samples": list(samples), "orientations": [s for s, _ in orients],
                "files": files,
            })
    if not args.dry_run:
        (work / BINS_JSON).write_text(json.dumps(bins, indent=1))
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"Made {len(bins)} bins ({nm} mass x {nt} t) x {len(orients)} orientations in {work}/")
    if not args.dry_run:
        for b in bins:
            empty = sorted({f for f in b["files"].values()
                            if n_events(work / b["name"] / f, tree) == 0})
            if empty:
                print(f"  ! {b['name']}: no events in {', '.join(empty)}; "
                      f"fits that need them will be skipped")


# --------------------------------------------------------------------------
# step 2: fitting
# --------------------------------------------------------------------------

def fill_template(template, subs):
    pattern = re.compile(r"\b(" + "|".join(map(re.escape, subs)) + r")\b")
    return pattern.sub(lambda m: subs[m.group(1)], template)


def fit_jobs(b, args):
    """The fits of one bin: (folder name, config text, files it reads)."""
    jobs = []
    if args.cfg:
        subs = {tok: f"../{f}" for tok, f in b["files"].items()}
        subs["FITNAME"] = COMBINED
        used = [f for tok, f in b["files"].items()
                if re.search(rf"\b{re.escape(tok)}\b", args.cfg_text)]
        jobs.append((COMBINED, fill_template(args.cfg_text, subs), used))
    if args.cfg_single:
        for suf in b["orientations"]:
            name = f"pol{suf}"
            subs = {s: f"../{b['files'][f'{s}_{suf}']}" for s in b["samples"]}
            subs.update({"FITNAME": name, "POLANGLE": f"polAngle_{suf}",
                         "POLVAL": f"polVal_{suf}"})
            used = [b["files"][f"{s}_{suf}"] for s in b["samples"]]
            jobs.append((name, fill_template(args.cfg_single_text, subs), used))
    return jobs


def set_start_values(cfg, params):
    """Start values of the SDME parameters in a fit config, from a previous fit."""
    names = {n for n, _ in SDMES}

    def repl(m):
        name = m.group(2)
        if name not in names or name not in params:
            return m.group(0)
        return f"{m.group(1)}{params[name]:.6g}"

    return re.sub(r"^(\s*parameter\s+(\w+)\s+)\S+", repl, cfg, flags=re.M)


def step_fit(args):
    work = Path(args.workdir)
    bins = select_bins(load_bins(work), args.bin)
    args.cfg_text = Path(args.cfg).read_text() if args.cfg else ""
    args.cfg_single_text = Path(args.cfg_single).read_text() if args.cfg_single else ""
    extra = []
    if args.rand is not None:
        extra += ["-r", str(args.rand)]
    if args.rand_seed is not None:
        extra += ["-rs", str(args.rand_seed)]
    if args.max_iter is not None:
        extra += ["-m", str(args.max_iter)]
    if args.min_ematrix is not None:
        extra += ["-e", str(args.min_ematrix)]
    if args.hesse:
        extra.append("-H")
    if args.minos:
        extra.append("-n")
    extra += shlex.split(args.fit_args or "")
    # single-orientation fits that start from the combined result run once:
    # randomizing them would throw the start values away
    start_combined = args.start_from_combined and args.cfg and args.cfg_single
    extra_single = extra
    if start_combined:
        extra_single, skip = [], False
        for a in extra:
            if skip:
                skip = False
            elif a in ("-r", "-rs"):
                skip = True
            else:
                extra_single.append(a)
    tree = args.tree or "kin"

    def fit_one(job):
        b, name, cfg, used = job
        d = work / b["name"] / name
        label = f"{b['name']}/{name}"
        empty = [] if args.dry_run else [f for f in used
                                         if n_events(work / b["name"] / f, tree) == 0]
        if empty:
            print(f"  - {label}: skipped, no events in {', '.join(empty)}")
            return label, "skipped"
        opts = extra
        if start_combined and name != COMBINED:
            opts = extra_single
            comb = main_fit_file(work / b["name"] / COMBINED, COMBINED)
            if comb:
                cfg = set_start_values(cfg, parse_fit(comb)["params"])
        d.mkdir(parents=True, exist_ok=True)
        for old in d.glob("*.fit"):
            old.unlink()                # a failed fit must not leave an old result behind
        (d / "fit.cfg").write_text(cfg)
        rc = run([args.fit_exe, "-c", "fit.cfg"] + opts, cwd=d,
                 log=d / "fit.log", dry=args.dry_run)
        if rc:
            print(f"  ! {label}: fit exited {rc}, see {d/'fit.log'}")
        elif not args.dry_run and not main_fit_file(d, name):
            print(f"  ! {label}: no converged fit (no .fit written), see {d/'fit.log'}")
            return label, "noconv"
        return label, rc

    jobs = [(b, *j) for b in bins for j in fit_jobs(b, args)]
    if not jobs:
        sys.exit("nothing to fit: set cfg and/or cfg_single")
    # the combined fits first, so the single-orientation fits can start from them
    phases = ([[j for j in jobs if j[1] == COMBINED], [j for j in jobs if j[1] != COMBINED]]
              if start_combined else [jobs])
    results = []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for phase in phases:
            results += list(pool.map(fit_one, phase))
    skipped = [n for n, rc in results if rc == "skipped"]
    noconv = [n for n, rc in results if rc == "noconv"]
    bad = [n for n, rc in results if rc and rc not in ("skipped", "noconv")]
    print(f"Fits done: {len(results) - len(bad) - len(skipped) - len(noconv)} ok, "
          f"{len(bad)} failed" + (f": {', '.join(bad)}" if bad else "")
          + (f", {len(noconv)} not converged: {', '.join(noconv)}" if noconv else "")
          + (f", {len(skipped)} skipped (empty bins)" if skipped else ""))


# --------------------------------------------------------------------------
# step 3: gathering
# --------------------------------------------------------------------------

def parse_fit(path):
    """Likelihood, parameter values and covariance matrix of an AmpTools .fit file."""
    lines = Path(path).read_text().splitlines()
    sec = {}
    for i, l in enumerate(lines):
        if l.startswith("+++ "):
            sec.setdefault(l.strip("+ ").strip(), i)
    out = {"likelihood": None, "params": {}, "cov": None}
    if "Likelihood Total and Partial Sums" in sec:
        out["likelihood"] = float(lines[sec["Likelihood Total and Partial Sums"] + 1].split()[0])
    i = sec["Parameter Values and Errors"] + 1
    n = int(lines[i]); i += 1
    names = []
    for _ in range(n):
        name, val = lines[i].split()[:2]; i += 1
        names.append(name)
        out["params"][name] = float(val)
    try:
        cov = [[float(x) for x in lines[i + r].split()] for r in range(n)]
        if all(len(row) == n for row in cov):
            out["cov"] = {(a, c): cov[r][k] for r, a in enumerate(names)
                          for k, c in enumerate(names)}
    except (ValueError, IndexError):
        pass
    return out


def sdme_values(fit):
    """SDMEs, derived quantities and their errors for one parsed fit."""
    p, cov = fit["params"], fit["cov"]
    row = {}
    for name, _ in SDMES:
        if name not in p:
            return None
        row[name] = p[name]
        row[name + "_err"] = math.sqrt(max(cov[(name, name)], 0.0)) if cov else float("nan")
    for name, (const, coef, _) in DERIVED.items():
        row[name] = const + sum(c * p[k] for k, c in coef.items())
        if cov:
            var = sum(ca * cb * cov[(a, b)] for a, ca in coef.items() for b, cb in coef.items())
            row[name + "_err"] = math.sqrt(max(var, 0.0))
        else:
            row[name + "_err"] = float("nan")
    for k, v in p.items():
        if k.startswith("parScale"):
            row[k] = v
            if cov:
                row[k + "_err"] = math.sqrt(max(cov[(k, k)], 0.0))
    return row


def step_gather(args):
    work = Path(args.workdir)
    bins = load_bins(work)
    rows = []
    for b in bins:
        for name in [COMBINED] + [f"pol{s}" for s in b["orientations"]]:
            d = work / b["name"] / name
            if not d.is_dir():
                continue
            row = {"bin": b["name"], "fit": name, "i_mass": b["i_mass"], "i_t": b["i_t"],
                   "mass_lo": b["mass_lo"], "mass_hi": b["mass_hi"],
                   "mass": 0.5 * (b["mass_lo"] + b["mass_hi"]),
                   "t_lo": b["t_lo"], "t_hi": b["t_hi"],
                   "t": None if b["t_lo"] is None else 0.5 * (b["t_lo"] + b["t_hi"])}
            fitfile = main_fit_file(d, name)
            vals = None
            if fitfile:
                fit = parse_fit(fitfile)
                vals = sdme_values(fit)
                row["likelihood"] = fit["likelihood"]
            else:
                print(f"  ! no .fit file in {d}")
            row["fit_ok"] = int(vals is not None)
            row.update(vals or {})
            rows.append(row)
    if not rows:
        sys.exit("no fit folders found; run the 'fit' step first")
    keys = list(dict.fromkeys(k for r in rows for k in r))
    out = work / RESULTS_CSV
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out} ({len(rows)} fits, {len(keys)} columns)")


# --------------------------------------------------------------------------
# step 4: plotting
# --------------------------------------------------------------------------

# colours of the reference plot: pol=0 black, 45 green, 90 red, 135 blue
FIT_STYLE = {
    "pol00": ("black", "pol=0"), "pol0": ("black", "pol=0"),
    "pol45": ("#009900", "pol=45"), "pol90": ("red", "pol=90"),
    "pol135": ("blue", "pol=135"), COMBINED: ("#9900cc", "all pol."),
}


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def read_reference(path):
    """Reference SDMEs vs -t (CSV: t or t_lo,t_hi, plus SDME columns; '#' lines are
    comments), with the derived quantities added."""
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith("#")]
    ref = []
    for r in csv.DictReader(lines):
        r = {k.strip(): num(v) for k, v in r.items()}
        if all(r.get(n) is not None for n, _ in SDMES):
            for name, (const, coef, _) in DERIVED.items():
                r[name] = const + sum(c * r[k] for k, c in coef.items())
        ref.append(r)
    return sorted(ref, key=lambda r: r.get("t") or 0.5 * (r["t_lo"] + r["t_hi"]))


def plot_grid(rows, quantities, ncols, out_pdf, args, ylims, ref_lines, reference=None):
    """One page per mass bin (x = -t) or a single page (x = mass, no t bins):
    a panel per quantity, one series per fit. `reference` (rows of
    read_reference) is drawn in gray: a box per bin if it has t_lo/t_hi,
    otherwise a line."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    use_t = any(r["t"] not in ("", None) for r in rows)
    xkey, lo, hi = ("t", "t_lo", "t_hi") if use_t else ("mass", "mass_lo", "mass_hi")
    pages = sorted({int(r["i_mass"]) for r in rows}) if use_t else [None]
    fits = [f for f in [COMBINED] + sorted({r["fit"] for r in rows} - {COMBINED},
                                          key=lambda s: int(re.sub(r"\D", "", s) or 0))
            if any(r["fit"] == f for r in rows)]
    if args.fits:
        fits = [f for f in fits if f in args.fits]
    nrows = math.ceil(len(quantities) / ncols)
    with PdfPages(out_pdf) as pdf:
        for page in pages:
            prow = [r for r in rows if page is None or int(r["i_mass"]) == page]
            fig, axes = plt.subplots(nrows, ncols, figsize=(3.6 * ncols, 3.4 * nrows),
                                     squeeze=False)
            for ax, (q, label) in zip(axes.flat, quantities):
                xs = sorted({float(r[xkey]) for r in prow})
                for f in fits:
                    pts = [r for r in prow if r["fit"] == f and num(r.get(q)) is not None]
                    if not pts:
                        continue
                    color, name = FIT_STYLE.get(f, (None, f))
                    x = [float(r[xkey]) for r in pts]
                    xerr = [0.5 * (float(r[hi]) - float(r[lo])) for r in pts]
                    ax.errorbar(x, [num(r[q]) for r in pts],
                                yerr=[num(r.get(q + "_err")) or 0.0 for r in pts],
                                xerr=xerr, fmt="o", ms=3, capsize=1.5, lw=1,
                                color=color, label=name)
                ax.axhline(0.0, color="#1f3fbf", lw=0.8)
                rpts = [r for r in (reference or []) if r.get(q) is not None] if use_t else []
                if rpts and rpts[0].get("t_lo") is not None:
                    for i, r in enumerate(rpts):
                        ax.fill_between([r["t_lo"], r["t_hi"]], r[q] - 0.01, r[q] + 0.01,
                                        color="gray", alpha=0.5, lw=0, zorder=0,
                                        label=args.reference_label if i == 0 else None)
                elif rpts:
                    ax.plot([r["t"] for r in rpts], [r[q] for r in rpts], color="gray",
                            lw=1.5, marker="s", ms=2.5, zorder=0, label=args.reference_label)
                if q in ref_lines:
                    ax.axhline(ref_lines[q], color="gray", lw=0.8, ls="--")
                # the given range, widened when points (with errors) fall outside it
                ylo, yhi = ylims.get(q, ylims["default"])
                vals = [(num(r[q]), num(r.get(q + "_err")) or 0.0) for r in prow
                        if (not args.fits or r["fit"] in args.fits) and num(r.get(q)) is not None]
                if vals:
                    lo_d, hi_d = min(v - e for v, e in vals), max(v + e for v, e in vals)
                    pad = 0.1 * (yhi - ylo)
                    ylo, yhi = min(ylo, lo_d - pad), max(yhi, hi_d + pad)
                ax.set_ylim(ylo, yhi)
                if xs:
                    ax.set_xlim(min(float(r[lo]) for r in prow), max(float(r[hi]) for r in prow))
                ax.text(0.06, 0.92, label, transform=ax.transAxes, fontsize=15,
                        fontweight="bold", va="top", zorder=5,
                        bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
                ax.set_xlabel(r"$-t$ [GeV$^2$]" if use_t else r"$M(\pi^-\pi^0)$ [GeV]")
                ax.tick_params(direction="in", top=True, right=True)
            for ax in list(axes.flat)[len(quantities):]:
                ax.axis("off")
            handles, labels = axes.flat[0].get_legend_handles_labels()
            for ax in axes.flat:
                h, l = ax.get_legend_handles_labels()
                if len(h) > len(handles):
                    handles, labels = h, l
            fig.legend(handles, labels, loc="upper center", ncol=len(labels),
                       frameon=False, fontsize=11)
            if page is not None:
                r0 = prow[0]
                fig.suptitle(f"{float(r0['mass_lo']):.3f} < M($\\pi^-\\pi^0$) < "
                             f"{float(r0['mass_hi']):.3f} GeV", y=0.005, va="bottom",
                             fontsize=10)
            fig.tight_layout(rect=(0, 0.02, 1, 0.96))
            pdf.savefig(fig)
            if args.format != "pdf":
                suffix = f"_m{page:03d}" if page is not None else ""
                fig.savefig(out_pdf.with_name(f"{out_pdf.stem}{suffix}.{args.format}"), dpi=150)
            plt.close(fig)
    print(f"Wrote {out_pdf}")


def step_plot(args):
    work = Path(args.workdir)
    with open(work / RESULTS_CSV) as f:
        rows = [r for r in csv.DictReader(f) if r.get("fit_ok") == "1"]
    if not rows:
        sys.exit("results.csv has no successful fits")
    outdir = work / "plots"
    outdir.mkdir(exist_ok=True)
    y = (-args.ylim, args.ylim)
    comps = [(k, DERIVED[k][2]) for k in ("rhoN_11", "rhoN_00", "rhoN_1m1", "ReRhoN_10",
                                          "rhoU_11", "rhoU_00", "rhoU_1m1", "ReRhoU_10")]
    # every fit, then (<name>_total.pdf) only the combined fit of all orientations
    total = [r for r in rows if r["fit"] == COMBINED]
    ref = read_reference(args.reference) if args.reference else None
    for sel, suffix in ((rows, ""), (total, "_total")):
        if not sel or (suffix and len(total) == len(rows)):
            continue
        plot_grid(sel, SDMES, 3, outdir / f"sdmes{suffix}.pdf", args, {"default": y}, {}, ref)
        plot_grid(sel, comps, 4, outdir / f"parity_components{suffix}.pdf", args,
                  {"default": (-0.2, 0.2), "rhoN_11": (0.3, 0.7)}, SCHC_NPE, ref)
        plot_grid(sel, [(k, DERIVED[k][2]) for k in ("P_sigma", "Sigma")], 2,
                  outdir / f"asymmetries{suffix}.pdf", args, {"default": (-0.5, 1.5)},
                  SCHC_NPE, ref)


# --------------------------------------------------------------------------
# step 5: fit quality (sdme_fit_plotter.py)
# --------------------------------------------------------------------------

def merge_pdfs(pdfs, out):
    """Merge per-bin PDFs into one file with pypdf, pdfunite or ghostscript."""
    try:
        from pypdf import PdfWriter
        w = PdfWriter()
        for p in pdfs:
            w.append(str(p))
        w.write(str(out))
        return True
    except ImportError:
        pass
    if shutil.which("pdfunite"):
        return run(["pdfunite"] + [str(p) for p in pdfs] + [str(out)]) == 0
    if shutil.which("gs"):
        return run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-sDEVICE=pdfwrite",
                    f"-sOutputFile={out}"] + [str(p) for p in pdfs]) == 0
    return False


def step_plotfit(args):
    work = Path(args.workdir)
    bins = select_bins(load_bins(work), args.bin)
    jobs = []
    for b in bins:
        for name in [COMBINED] + [f"pol{s}" for s in b["orientations"]]:
            d = work / b["name"] / name
            fitfile = main_fit_file(d, name) if d.is_dir() else None
            if fitfile:
                jobs.append((b, name, d, fitfile))
    if not jobs:
        sys.exit("no .fit files found; run the 'fit' step first")

    def one(job):
        b, name, d, fitfile = job
        out = d / "fitplots.pdf"
        if not args.redo_plots and out.exists() and out.stat().st_mtime >= fitfile.stat().st_mtime:
            return name, out                     # already made (e.g. by a farm job)
        title = f"{b['name']}/{name}: {b['mass_lo']:.3f} < M < {b['mass_hi']:.3f} GeV"
        if b["t_lo"] is not None:
            title += f", {b['t_lo']:.2f} < -t < {b['t_hi']:.2f}"
        rc = run([sys.executable, args.plotter, fitfile.name, "-o", "fitplots.pdf",
                  "-t", title, "--tree", args.tree or "kin"],
                 cwd=d, log=d / "plotfit.log", dry=args.dry_run)
        if rc:
            print(f"  ! {b['name']}/{name}: plotter exited {rc}, see {d/'plotfit.log'}")
            return name, None
        return name, d / "fitplots.pdf"

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        done = [r for r in pool.map(one, jobs) if r[1] and r[1].exists()]
    if args.dry_run or not done or args.bin:
        return
    outdir = work / "plots"
    outdir.mkdir(exist_ok=True)
    for name in dict.fromkeys(n for n, _ in done):
        pdfs = [p for n, p in done if n == name]
        out = outdir / f"fit_quality_{name}.pdf"
        if merge_pdfs(pdfs, out):
            print(f"Wrote {out} ({len(pdfs)} bins)")
        else:
            print(f"Per-bin fit-quality plots in <bin>/{name}/fitplots.pdf; install pypdf "
                  f"to also get {out.name}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description=__doc__, fromfile_prefix_chars="@",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step", required=True)
    here = Path(__file__).resolve().parent

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-c", "--config", help="INI file with settings (see example_sdme.ini)")
    common.add_argument("--workdir", default="SDME_fit", help="where bin directories go")
    common.add_argument("--dry-run", action="store_true", help="print commands only")
    common.add_argument("--jobs", type=int, default=1, help="fits run in parallel")
    common.add_argument("--tree", help="tree name (default kin)")

    b = argparse.ArgumentParser(add_help=False)
    b.add_argument("--sample", action="append", nargs="+", default=[], metavar="TOKEN=FILE",
                   help="input sample(s), e.g. --sample DATAFILE=data.root ACCMCFILE=acc.root "
                        "GENMCFILE=gen.root")
    b.add_argument("--orientations", nargs="+", default=["00:0", "45:45", "90:90", "135:135"],
                   metavar="SUFFIX:POL",
                   help="orientations: token suffix and value of the pol branch "
                        "(default 00:0 45:45 90:90 135:135)")
    b.add_argument("--split", nargs="+", default=["DATAFILE"], metavar="TOKEN",
                   help="samples split by orientation (default DATAFILE); the others "
                        "are used whole for every orientation")
    b.add_argument("--pol-branch", default="pol", help="Int_t branch with the orientation")
    b.add_argument("--mass", nargs=3, type=float, metavar=("LO", "HI", "N"),
                   help="uniform M(pi- pi0) binning")
    b.add_argument("--t", nargs=3, type=float, metavar=("LO", "HI", "N"),
                   help="uniform -t binning (omit for mass bins only)")
    b.add_argument("--t-log", action="store_true", help="logarithmic t bins; LO must be > 0")
    b.add_argument("--splitter", default=str(here / "split_bins.C"), help="split_bins.C")
    b.add_argument("--mass-indices", default="1,2",
                   help="final-state indices whose mass is binned (default 1,2 = Pi- Pi0)")
    b.add_argument("--root-exe", default="root", help="ROOT executable for split_bins.C")

    f = argparse.ArgumentParser(add_help=False)
    f.add_argument("--cfg", help="template for the combined fit of all orientations "
                   "(fit_pi0pim_SDME_4pol_TEMPLATE.cfg); empty = skip")
    f.add_argument("--cfg-single", help="template for one fit per orientation "
                   "(fit_pi0pim_SDME_1pol_TEMPLATE.cfg); empty = skip")
    f.add_argument("--fit-exe", default="fit", help="fit executable (fit, fitMPI, ...)")
    f.add_argument("-r", "--rand", type=int, help="fit -r: number of randomized fits")
    f.add_argument("--rand-seed", type=int, help="fit -rs: random seed")
    f.add_argument("--max-iter", type=int, help="fit -m: maximum number of iterations")
    f.add_argument("--min-ematrix", type=int, choices=range(4), help="fit -e")
    f.add_argument("--hesse", action="store_true", help="fit -H: run HESSE")
    f.add_argument("--minos", action="store_true", help="fit -n: use MINOS")
    f.add_argument("--fit-args", default="", help='any other fit options')
    f.add_argument("--no-start-from-combined", dest="start_from_combined",
                   action="store_false",
                   help="do not start the single-orientation fits from the combined fit "
                        "result (and do not run them before it)")

    sel = argparse.ArgumentParser(add_help=False)
    sel.add_argument("--bin", nargs="+", metavar="NAME",
                     help="only these bins, e.g. bin_m000_t03 (used by the farm jobs)")

    g = argparse.ArgumentParser(add_help=False)

    pl = argparse.ArgumentParser(add_help=False)
    pl.add_argument("--format", default="pdf", help="pdf, or png etc. for extra image files")
    pl.add_argument("--ylim", type=float, default=0.5, help="SDME plots: y from -YLIM to YLIM")
    pl.add_argument("--fits", nargs="+", help="plot only these fits, e.g. all pol00 pol90")
    pl.add_argument("--reference", help="CSV of SDMEs vs -t drawn in gray, e.g. published "
                    "values or the truth.csv of a closure test")
    pl.add_argument("--reference-label", default="reference", help="legend label of --reference")

    pf = argparse.ArgumentParser(add_help=False)
    pf.add_argument("--plotter", default=str(here / "sdme_fit_plotter.py"),
                    help="fit-quality plotter (default: sdme_fit_plotter.py next to this script)")
    pf.add_argument("--no-plotfit", dest="plotfit", action="store_false",
                    help="with 'all': skip the fit-quality plots")
    pf.add_argument("--redo-plots", action="store_true",
                    help="remake fitplots.pdf even when it is newer than the .fit file")

    subs = {
        "bin": sub.add_parser("bin", parents=[common, b], help="split samples into bins"),
        "fit": sub.add_parser("fit", parents=[common, f, sel], help="run the fits in each bin"),
        "gather": sub.add_parser("gather", parents=[common, g], help="collect results to CSV"),
        "plot": sub.add_parser("plot", parents=[common, pl], help="plot the results"),
        "plotfit": sub.add_parser("plotfit", parents=[common, pf, sel],
                                  help="data vs. fit plots for every fit"),
        "all": sub.add_parser("all", parents=[common, b, f, g, pl, pf, sel],
                              help="run bin, fit, gather, plot, plotfit"),
    }
    return p, subs


# settings whose relative paths are taken relative to the config file
CONFIG_PATH_KEYS = {"cfg", "cfg_single", "splitter", "plotter", "reference"}


def config_defaults(path, subparser, all_subparsers):
    """
    Read an INI file and turn it into defaults for `subparser`.  Section names
    are only for readability, except [samples], whose entries are TOKEN = file.
    Keys are option names with - or _ (e.g. fit_args, max-iter).
    """
    path = Path(path)
    if not path.exists():
        sys.exit(f"config file {path} not found")
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"), interpolation=None)
    cp.optionxform = str            # keep DATAFILE etc. upper case
    cp.read(path)
    base = path.resolve().parent
    here = Path(__file__).resolve().parent

    def localize(v):
        # relative to the INI file, else (helper scripts) next to this script
        p = Path(os.path.expanduser(v))
        if not p.is_absolute():
            for d in (base, here):
                if (d / p).exists():
                    return str(d / p)
        return str(p)

    actions = {a.dest: a for a in subparser._actions}
    known = {a.dest for sp in all_subparsers.values() for a in sp._actions}
    out = {}
    for section in cp.sections():
        if section.lower() == "farm":
            continue                    # job settings, read by submit_sdme_jobs.py
        if section.lower() == "samples":
            out["sample"] = [[f"{tok}={localize(v)}" for tok, v in cp.items(section)]]
            continue
        for key, raw in cp.items(section):
            dest = key.strip().lower().replace("-", "_")
            if dest == "r":
                dest = "rand"
            if dest not in known:
                sys.exit(f"{path}: unknown setting '{key}' in [{section}]")
            if dest not in actions:
                continue
            a = actions[dest]
            if raw.strip() == "":
                if dest in ("cfg", "cfg_single"):
                    out[dest] = None            # empty = skip that kind of fit
                continue
            if isinstance(a, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
                try:
                    out[dest] = cp.getboolean(section, key)
                except ValueError:
                    sys.exit(f"{path}: '{key}' should be true or false")
                continue
            vals = raw.split() if a.nargs not in (None, "?") else [raw.strip()]
            try:
                vals = [a.type(v) if a.type else v for v in vals]
            except ValueError:
                sys.exit(f"{path}: bad value for '{key}': {raw}")
            if dest in CONFIG_PATH_KEYS:
                vals = [localize(v) for v in vals]
            out[dest] = vals if a.nargs not in (None, "?") else vals[0]
    return out


def main():
    parser, subs = build_parser()
    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("step", nargs="?")
    pre.add_argument("-c", "--config")
    known, _ = pre.parse_known_args()
    if known.config and known.step in subs:
        subs[known.step].set_defaults(**config_defaults(known.config, subs[known.step], subs))
    args = parser.parse_args()
    steps = {"bin": [step_bin], "fit": [step_fit], "gather": [step_gather],
             "plot": [step_plot], "plotfit": [step_plotfit],
             "all": [step_bin, step_fit, step_gather, step_plot]
                    + ([step_plotfit] if getattr(args, "plotfit", True) else [])}[args.step]
    if step_bin in steps:
        if not args.sample or not args.mass:
            sys.exit("bin needs --sample and --mass")
        for name, rng in (("mass", args.mass), ("t", args.t)):
            if rng and (rng[1] <= rng[0] or rng[2] < 1):
                sys.exit(f"{name} = {rng[0]:g} {rng[1]:g} {rng[2]:g}: need LOW HIGH NBINS "
                         f"with HIGH > LOW and NBINS >= 1")
        if args.t_log and (not args.t or args.t[0] <= 0):
            sys.exit("--t-log needs --t with LO > 0")
    if step_fit in steps and not (args.cfg or args.cfg_single):
        sys.exit("fit needs --cfg and/or --cfg-single")
    for step in steps:
        step(args)


if __name__ == "__main__":
    main()
