#!/usr/bin/env python3
"""
Automate a binned AmpTools fit of gamma p -> eta' pi- pi+ p
(reaction EtaPi Beam Proton EtaPrime Pi- Pi+, i.e. eta' pi- recoiling
against a p pi+ / Delta++ system):

  1. bin     split each sample (data, acc MC, gen MC, bkg) into M(eta' pi-) x -t
             bins with split_bins.C (ROOT macro, one pass per sample).
             -t = -(beam - eta' pi-)^2, the momentum transfer to the p pi+ system.
             (split_mass takes every particle after the recoil, which here includes
             the pi+, so it is not the default.)
  2. fit     write a per-bin config from a template and run `fit -c` in each bin
  3. gather  run convert_to_csv.py on every bin's .fit file and collect one
             results.csv with the bin edges added
  4. plot    plot every gathered quantity vs. mass, one series per t bin
  5. plotfit (optional) run etapi_fit_plotter.py on each bin's fit; each bin
             gets fitplots.root + fitplots.pdf (data vs. fit, including M(p pi+)),
             and all bins
             are merged into plots/fit_quality.pdf

Wave sets: list them in a [wavesets] section of the INI file, e.g.
      [wavesets]
      P1D1 = D1+ P1+
      SPD  = D1+ S0+ P1+ D0+ D2+ S0- D1-
      ALLD = S0 P1+ D*                # shorthand, see below
Shorthand: a wave without + or - means both reflectivities (D1 = D1+ D1-),
* and ? are wildcards (D*+ = every positive D wave, incl. Dn1+, Dn2+), and
`all` is every wave in the template.  The listed waves float and every other wave in the template is fixed to 0;
the first listed wave of each reflectivity (+/-) is the real reference wave.
The data is binned once in <workdir>/, and each wave set is fit in its own
<workdir>/<name>/ (bins link to the shared ROOT files) with its own
results.csv and plots; gather also writes <workdir>/results_wavesets.csv
and plot draws every quantity with one curve per wave set in
<workdir>/plots_wavesets/.
--waveset NAME ... runs only some of the sets.

Each step is a subcommand; `all` runs them in order.  Settings can go in an
INI config file (see example_run.ini) passed with -c; options given on the
command line override the file:

  amptools_auto.py all -c myfit.ini
  amptools_auto.py fit -c myfit.ini -r 20        # same, but 20 random starts

Everything can also be given on the command line, e.g.:
  amptools_auto.py all --workdir fits \\
      --sample DATAFILE=data.root --sample ACCMCFILE=acc.root \\
      --sample GENMCFILE=gen.root --sample BKGNDFILE=bkg.root \\
      --mass 1.1 2.6 30 --t 0.1 0.3 1 \\
      --cfg template.cfg -r 10 --jobs 8 \\
      --data-info DATAFILE
"""

import argparse
import configparser
import csv
import json
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


def main_fit_file(d):
    """With -r the per-attempt results are <name>_<k>.fit; prefer <name>.fit."""
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


def empty_samples(d, b, tree="kin"):
    """Samples of bin b (in folder d) whose tree has no events."""
    return [tok for tok, f in b["files"].items() if n_events(d / f, tree) == 0]


def load_bins(workdir):
    path = Path(workdir) / BINS_JSON
    if not path.exists():
        sys.exit(f"{path} not found; run the 'bin' step first.")
    return json.loads(path.read_text())


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

    tree_opt = ["-T", args.tree] if args.tree else []
    use_macro = args.splitter.endswith(".C")
    samples = {}
    for item in (x for group in args.sample for x in group):
        if "=" not in item:
            sys.exit(f"--sample expects TOKEN=FILE, got '{item}'")
        token, path = item.split("=", 1)
        samples[token] = path

    for token, path in samples.items():
        src = Path(path).resolve()
        if use_macro:
            # one pass: <tmp>/<token>_m<i>_t<j>.root
            call = (f'{Path(args.splitter).resolve()}("{src}", "{tmp / token}", '
                    f'{mlo}, {mhi}, {nm}, '
                    + (f'{tlo}, {thi}, {nt}, {"true" if args.t_log else "false"}, '
                       if t_edges else '0, 0, 0, false, ')
                    + f'"{args.tree or "kin"}", "{args.mass_indices}")')
            rc = run([args.root_exe, "-l", "-b", "-q", call], dry=args.dry_run)
            if rc:
                sys.exit(f"split_bins.C failed on {src} (exit {rc})")
            for i in range(nm):
                for j in range(nt):
                    part = tmp / f"{token}_m{i}_t{j}.root"
                    dest = work / bin_name(i, j) / f"{token.lower()}.root"
                    if not args.dry_run:
                        if not part.exists():
                            sys.exit(f"split_bins.C did not write {part}")
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(part), dest)
            continue

        # split_mass / split_t (mass of every final-state particle after the recoil)
        mbase = tmp / f"{token}_m"
        # split_mass writes <base>_<i>.root, i = 0..nm-1 (uniform bins)
        rc = run([args.split_mass, str(src), str(mbase),
                  str(mlo), str(mhi), str(nm)] + tree_opt, dry=args.dry_run)
        if rc:
            sys.exit(f"split_mass failed on {src} (exit {rc})")

        for i in range(nm):
            mfile = Path(f"{mbase}_{i}.root")
            if t_edges:
                tbase = tmp / f"{token}_m{i}_t"
                rc = run([args.split_t, str(mfile), str(tbase),
                          str(tlo), str(thi), str(nt)] + tree_opt
                         + (["-e"] if args.t_log else []),
                         dry=args.dry_run)
                if rc:
                    sys.exit(f"split_t failed on {mfile} (exit {rc})")
            for j in range(nt):
                part = Path(f"{tbase}_{j}.root") if t_edges else mfile
                dest = work / bin_name(i, j) / f"{token.lower()}.root"
                if not args.dry_run:
                    dest.parent.mkdir(parents=True, exist_ok=True)
                    shutil.move(str(part), dest)

    bins = []
    for i in range(nm):
        for j in range(nt):
            bins.append({
                "name": bin_name(i, j), "i_mass": i, "i_t": j,
                "mass_lo": m_edges[i], "mass_hi": m_edges[i + 1],
                "t_lo": t_edges[j] if t_edges else None,
                "t_hi": t_edges[j + 1] if t_edges else None,
                "files": {tok: f"{tok.lower()}.root" for tok in samples},
            })
    if not args.dry_run:
        (work / BINS_JSON).write_text(json.dumps(bins, indent=1))
        shutil.rmtree(tmp, ignore_errors=True)
    print(f"Made {len(bins)} bins ({nm} mass x {nt} t) in {work}/")
    if not args.dry_run:
        for b in bins:
            empty = empty_samples(work / b["name"], b, args.tree or "kin")
            if empty:
                print(f"  ! {b['name']} ({b['mass_lo']:.3f}-{b['mass_hi']:.3f} GeV): "
                      f"no events in {', '.join(empty)}; it will be skipped by the fit")


# --------------------------------------------------------------------------
# step 2: fitting
# --------------------------------------------------------------------------

INIT_RE = re.compile(r"^([ \t]*initialize[ \t]+)(\S+)::(\S+)::(\S+)([ \t].*)?$", re.M)


def template_waves(cfg):
    """Wave names (last part of reaction::sum::wave) in initialize-line order."""
    return list(dict.fromkeys(m.group(4) for m in INIT_RE.finditer(cfg)))


def expand_waves(tokens, known):
    """Expand shorthand: `all`, wildcards (D*+), and names without +/- (D1 ->
    D1+ D1-).  Keeps the listed order; wildcards follow the template order."""
    import fnmatch
    out = []
    for tok in tokens:
        if tok.lower() == "all":
            pat = "*"
        elif tok[-1] in "+-":
            pat = tok
        else:
            pat = tok + "[+-]"
        hits = [w for w in known if fnmatch.fnmatchcase(w, pat)]
        if tok[-1] not in "+-" and "*" not in tok and "?" not in tok:
            hits.sort(key=lambda w: w[-1] != "+")      # D1+ before D1-
        if not hits:
            sys.exit(f"wave '{tok}' matches nothing in the template "
                     f"(known: {' '.join(known)})")
        out += [w for w in hits if w not in out]
    return out


def apply_waveset(cfg, waves):
    """Let `waves` float and fix every other wave to 0.  The first wave of
    each reflectivity (last character + or -) in `waves` is the real
    reference; floating waves start at 10+10i, references at 100."""
    in_cfg = set(template_waves(cfg))
    missing = [w for w in waves if w not in in_cfg]
    if missing:
        sys.exit(f"wave(s) {' '.join(missing)} not in the template's initialize lines "
                 f"(known: {' '.join(sorted(in_cfg))})")
    ref = {}
    for w in waves:
        ref.setdefault(w[-1], w)

    def fix(m):
        wave = m.group(4)
        if wave not in waves:
            init = "cartesian 0.0 0.0 fixed"
        elif ref.get(wave[-1]) == wave:
            init = "cartesian 100.0 0.0 real"
        else:
            init = "cartesian 10.0 10.0"
        return f"{m.group(1)}{m.group(2)}::{m.group(3)}::{wave} {init}"
    return INIT_RE.sub(fix, cfg)


def make_config(template, b, template_dir, param_init=True, waves=None):
    """Replace sample tokens (DATAFILE, ...) and bin tokens in the template.
    With param_init False, `include param_init...` lines are commented out.
    With `waves`, only those waves float (see apply_waveset)."""
    subs = dict(b["files"])
    subs.update({
        "FITNAME": b["name"], "BINNAME": b["name"], "NIFILE": f"{b['name']}.ni",
        "MASSLOW": f"{b['mass_lo']:g}", "MASSHIGH": f"{b['mass_hi']:g}",
        "TLOW": "" if b["t_lo"] is None else f"{b['t_lo']:g}",
        "THIGH": "" if b["t_hi"] is None else f"{b['t_hi']:g}",
    })
    pattern = re.compile(r"\b(" + "|".join(map(re.escape, subs)) + r")\b")
    cfg = pattern.sub(lambda m: subs[m.group(1)], template)
    if waves:
        cfg = apply_waveset(cfg, waves)

    # fits run inside the bin directory, so point `include <file>` lines at
    # files that sit next to the template
    def fix_include(m):
        inc = template_dir / m.group(2)
        if not param_init and inc.name.startswith("param_init"):
            return "# param_init.cfg left out by amptools_auto (--param-init)"
        return m.group(1) + str(inc.resolve()) if inc.exists() else m.group(0)
    return re.sub(r"^([ \t]*include[ \t]+)(\S+)", fix_include, cfg, flags=re.M)


def step_fit(args):
    work = Path(args.workdir)
    bins = load_bins(work)
    template = Path(args.cfg).read_text()
    template_dir = Path(args.cfg).resolve().parent
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
    # param_init.cfg sets starting values; not needed with randomized starts
    if args.param_init == "auto":
        use_init = not args.rand
    else:
        use_init = args.param_init == "yes"
    if not use_init:
        print("Not including param_init.cfg in the fit configs")

    def fit_one(b):
        d = work / b["name"]
        empty = [] if args.dry_run else empty_samples(d, b, getattr(args, "tree", None) or "kin")
        if empty:
            print(f"  - {b['name']}: skipped, no events in {', '.join(empty)}")
            return b["name"], "skipped"
        (d / "fit.cfg").write_text(make_config(template, b, template_dir, use_init,
                                                     getattr(args, "waves", None)))
        rc = run([args.fit_exe, "-c", "fit.cfg"] + extra, cwd=d,
                 log=d / "fit.log", dry=args.dry_run)
        if rc:
            print(f"  ! {b['name']}: fit exited {rc}, see {d/'fit.log'}")
        return b["name"], rc

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        results = list(pool.map(fit_one, bins))
    skipped = [n for n, rc in results if rc == "skipped"]
    bad = [n for n, rc in results if rc and rc != "skipped"]
    print(f"Fits done: {len(results) - len(bad) - len(skipped)} ok, {len(bad)} failed"
          + (f": {', '.join(bad)}" if bad else "")
          + (f", {len(skipped)} skipped (empty bins)" if skipped else ""))


# --------------------------------------------------------------------------
# step 3: gathering
# --------------------------------------------------------------------------

def run_convert(args, files, out_csv, extra):
    """Run convert_to_csv.py on a list of .fit or .root files; return {abspath: row}."""
    if not files:
        return {}
    list_file = out_csv.with_suffix(".list.txt")
    list_file.write_text("\n".join(files) + "\n")
    cmd = [sys.executable, args.convert_script, "-i", str(list_file),
           "-o", str(out_csv), "--no-mass", "--no-sort"] + extra
    rc = run(cmd)
    list_file.unlink()
    if rc or not out_csv.exists():
        sys.exit(f"convert_to_csv.py failed (exit {rc})")
    with open(out_csv) as f:
        return {os.path.abspath(r["file"]): r for r in csv.DictReader(f)}


def step_gather(args):
    work = Path(args.workdir)
    bins = load_bins(work)

    fit_of = {}
    for b in bins:
        d = work / b["name"]
        fitfile = main_fit_file(d)
        if fitfile:
            fit_of[b["name"]] = os.path.abspath(fitfile)
        else:
            print(f"  ! no .fit file in {d}")

    extra = shlex.split(args.convert_args or "")
    if args.acceptance_corrected:
        extra.append("-a")
    fits = run_convert(args, list(fit_of.values()), work / "fits_raw.csv", extra)

    data_of, data = {}, {}
    if args.data_info:
        for b in bins:
            if args.data_info in b["files"]:
                f = os.path.abspath(work / b["name"] / b["files"][args.data_info])
                if os.path.exists(f):
                    data_of[b["name"]] = f
        data = run_convert(args, list(data_of.values()), work / "data_raw.csv",
                           shlex.split(args.data_info_args or ""))

    rows = []
    for b in bins:
        row = {"bin": b["name"], "i_mass": b["i_mass"], "i_t": b["i_t"],
               "mass_lo": b["mass_lo"], "mass_hi": b["mass_hi"],
               "mass": 0.5 * (b["mass_lo"] + b["mass_hi"]),
               "t_lo": b["t_lo"], "t_hi": b["t_hi"]}
        fr = fits.get(fit_of.get(b["name"]))
        row["fit_ok"] = int(fr is not None)
        if fr:
            row.update(fr)
        dr = data.get(data_of.get(b["name"]))
        if dr:
            row.update({f"data_{k}": v for k, v in dr.items() if k != "file"})
        rows.append(row)

    keys = list(dict.fromkeys(k for r in rows for k in r))
    out = work / RESULTS_CSV
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out} ({len(rows)} bins, {len(keys)} columns)")


# --------------------------------------------------------------------------
# step 4: plotting
# --------------------------------------------------------------------------

BIN_COLS = {"bin", "i_mass", "i_t", "mass_lo", "mass_hi", "mass",
            "t_lo", "t_hi", "fit_ok", "file"}


def num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if x != x else x


def step_plot(args):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    work = Path(args.workdir)
    with open(work / RESULTS_CSV) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        sys.exit("results.csv is empty")

    quantities = [k for k in rows[0]
                  if k not in BIN_COLS and not k.endswith("_err")
                  and any(num(r.get(k)) is not None for r in rows)]
    if args.only:
        quantities = [q for q in quantities if q in args.only]
    outdir = work / "plots"
    outdir.mkdir(exist_ok=True)

    t_bins = sorted({int(r["i_t"]) for r in rows})
    for q in quantities:
        fig, ax = plt.subplots(figsize=(7, 4.5))
        for j in t_bins:
            pts = [r for r in rows if int(r["i_t"]) == j and num(r.get(q)) is not None]
            if not pts:
                continue
            x = [float(r["mass"]) for r in pts]
            y = [num(r[q]) for r in pts]
            err = [num(r.get(q + "_err")) or 0.0 for r in pts]
            r0 = pts[0]
            label = (f"{float(r0['t_lo']):.2f} < -t < {float(r0['t_hi']):.2f}"
                     if r0["t_lo"] else None)
            ax.errorbar(x, y, yerr=err, fmt="o", ms=4, capsize=2, label=label)
        ax.set_xlabel(args.mass_label)
        ax.set_ylabel(q)
        if len(t_bins) > 1:
            ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / f"{q}.{args.format}")
        plt.close(fig)
    print(f"Wrote {len(quantities)} plots to {outdir}/")


# --------------------------------------------------------------------------
# step 5 (optional): per-bin fit quality (etapi_fit_plotter.py or etapi_delta_plotter)
# --------------------------------------------------------------------------

PLOTFIT_ROOT = "fitplots.root"


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
    bins = load_bins(work)
    if args.plotter.endswith(".py"):
        # standalone Python plotter (etapi_fit_plotter.py)
        exe = [sys.executable, args.plotter]
        opts = (["--per-sum"] if args.per_sum else []) + (["--with-gen"] if args.with_gen else [])
    else:
        # compiled etapi_delta_plotter
        exe = [args.plotter]
        opts = (["-s"] if not args.per_sum else []) + (["-G"] if not args.with_gen else [])

    def one(b):
        d = work / b["name"]
        fitfile = main_fit_file(d)
        if not fitfile:
            print(f"  ! no .fit file in {d}")
            return None
        title = f"{b['mass_lo']:.3f} < M < {b['mass_hi']:.3f} GeV"
        if b["t_lo"] is not None:
            title += f", {b['t_lo']:.2f} < -t < {b['t_hi']:.2f}"
        rc = run(exe + [fitfile.name, "-o", PLOTFIT_ROOT, "-t", title] + opts,
                 cwd=d, log=d / "plotfit.log", dry=args.dry_run)
        if rc:
            print(f"  ! {b['name']}: plotter exited {rc}, see {d/'plotfit.log'}")
            return None
        return d / PLOTFIT_ROOT.replace(".root", ".pdf")

    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        pdfs = [p for p in pool.map(one, bins) if p and p.exists()]
    if args.dry_run or not pdfs:
        return

    out = work / "plots" / "fit_quality.pdf"
    out.parent.mkdir(exist_ok=True)
    if merge_pdfs(pdfs, out):
        print(f"Per-bin histograms in <bin>/{PLOTFIT_ROOT}; all bins in {out}")
    else:
        print(f"Per-bin histograms and PDFs in <bin>/; install pypdf to also get "
              f"one combined {out.name}")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description=__doc__, fromfile_prefix_chars="@",
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="step", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-c", "--config", help="INI file with settings (see example_run.ini)")
    common.add_argument("--workdir", default="fits", help="where bin directories go")
    common.add_argument("--dry-run", action="store_true", help="print commands only")
    common.add_argument("--jobs", type=int, default=1, help="bins to process in parallel")
    common.add_argument("--waveset", nargs="+", metavar="NAME",
                        help="run only these wave sets from the [wavesets] section")
    common.add_argument("--wavesets", default=None, help=argparse.SUPPRESS)

    b = argparse.ArgumentParser(add_help=False)
    b.add_argument("--sample", action="append", nargs="+", default=[], metavar="TOKEN=FILE",
                   help="input sample(s); TOKEN is the word replaced in the config "
                        "template, e.g. --sample DATAFILE=data.root ACCMCFILE=acc.root")
    b.add_argument("--mass", nargs=3, type=float, metavar=("LO", "HI", "N"),
                   help="uniform mass binning")
    b.add_argument("--t", nargs=3, type=float, metavar=("LO", "HI", "N"),
                   help="uniform -t binning (omit for mass bins only)")
    b.add_argument("--t-log", action="store_true",
                   help="logarithmic t bins (split_t -e); LO must be > 0")
    b.add_argument("--tree", help="tree name (default kin)")
    b.add_argument("--splitter",
                   default=str(Path(__file__).resolve().parent / "split_bins.C"),
                   help="split_bins.C (default, next to this script) or 'split_mass' "
                        "to use split_mass / split_t instead")
    b.add_argument("--mass-indices", default="1,2",
                   help="final-state indices whose mass is binned (split_bins.C; "
                        "default 1,2 = EtaPrime Pi-, proton = 0, Pi+ = 3)")
    b.add_argument("--root-exe", default="root", help="ROOT executable for split_bins.C")
    b.add_argument("--split-mass", default="split_mass", help="split_mass executable")
    b.add_argument("--split-t", default="split_t", help="split_t executable")

    f = argparse.ArgumentParser(add_help=False)
    f.add_argument("--cfg", help="config template, e.g. fit_etaprpim_SPD_TEMPLATE.cfg; "
                   "FITNAME, DATAFILE, ACCMCFILE, GENMCFILE, NIFILE (and any other "
                   "--sample TOKEN) are filled in per bin")
    f.add_argument("--fit-exe", default="fit", help="fit executable (fit, fitMPI, ...)")
    f.add_argument("-r", "--rand", type=int, help="fit -r: number of fits with randomized starting parameters")
    f.add_argument("--rand-seed", type=int, help="fit -rs: random seed for the randomized fits")
    f.add_argument("--max-iter", type=int, help="fit -m: maximum number of fit iterations")
    f.add_argument("--min-ematrix", type=int, choices=range(4), help="fit -e: minimum error matrix status")
    f.add_argument("--hesse", action="store_true", help="fit -H: run HESSE after minimization")
    f.add_argument("--minos", action="store_true", help="fit -n: use MINOS instead of MIGRAD")
    f.add_argument("--param-init", choices=["auto", "yes", "no"], default="auto",
                   help="keep the template's `include param_init...` line: "
                        "auto (default) drops it when -r is used")
    f.add_argument("--fit-args", default="", help='any other fit options, e.g. --fit-args="-s seed.txt"')

    g = argparse.ArgumentParser(add_help=False)
    g.add_argument("--convert-script",
                   default=str(Path(__file__).resolve().parent / "convert_to_csv.py"),
                   help="convert_to_csv.py (default: the one next to this script)")
    g.add_argument("-a", "--acceptance-corrected", action="store_true",
                   help="acceptance-corrected intensities (convert_to_csv.py -a)")
    g.add_argument("--convert-args", default="",
                   help='other convert_to_csv.py options for .fit files, '
                        'e.g. --convert-args="--phases P1+:D1+ --covariance"')
    g.add_argument("--data-info", metavar="TOKEN",
                   help="also summarize this sample's per-bin ROOT file "
                        "(e.g. DATAFILE); adds data_* columns")
    g.add_argument("--data-info-args", default="",
                   help='convert_to_csv.py options for ROOT files, e.g. "--indices 1,2"')

    pl = argparse.ArgumentParser(add_help=False)
    pl.add_argument("--only", nargs="+", help="plot only these quantities")
    pl.add_argument("--format", default="pdf", help="pdf, png, ...")
    pl.add_argument("--mass-label", default=r"$M(\eta'\pi^-)$ [GeV]")

    pf = argparse.ArgumentParser(add_help=False)
    pf.add_argument("--plotfit", action="store_true",
                    help="with 'all': also run the plotfit step")
    pf.add_argument("--plotter",
                    default=str(Path(__file__).resolve().parent / "etapi_fit_plotter.py"),
                    help="etapi_fit_plotter.py (default, next to this script) or a "
                         "compiled etapi_delta_plotter")
    pf.add_argument("--per-sum", action="store_true",
                    help="also make per-sum histograms (slow; off by default)")
    pf.add_argument("--with-gen", action="store_true",
                    help="also project generated MC (slow; off by default)")

    subs = {
        "bin": sub.add_parser("bin", parents=[common, b], help="split samples into bins"),
        "fit": sub.add_parser("fit", parents=[common, f], help="run fit in each bin"),
        "gather": sub.add_parser("gather", parents=[common, g], help="collect results to CSV"),
        "plot": sub.add_parser("plot", parents=[common, pl], help="plot gathered results"),
        "plotfit": sub.add_parser("plotfit", parents=[common, pf],
                                  help="data vs. fit plots for each bin (etapi_delta_plotter)"),
        "all": sub.add_parser("all", parents=[common, b, f, g, pl, pf],
                              help="run bin, fit, gather, plot (and plotfit if enabled)"),
    }
    return p, subs


# --------------------------------------------------------------------------
# config file
# --------------------------------------------------------------------------

# settings whose relative paths are taken relative to the config file
CONFIG_PATH_KEYS = {"cfg", "convert_script", "plotter", "splitter"}


def config_defaults(path, subparser, all_subparsers):
    """
    Read an INI file and turn it into defaults for `subparser`.  Section names
    are only for readability, except [samples], whose entries are TOKEN = file.
    Keys are option names with - or _ (e.g. fit_args, max-iter).
    """
    path = Path(path)
    if not path.exists():
        sys.exit(f"config file {path} not found")
    cp = configparser.ConfigParser(inline_comment_prefixes=("#", ";"),
                                   interpolation=None)
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
        if section.lower() == "wavesets":
            out["wavesets"] = {name: v.split() for name, v in cp.items(section)}
            for name, waves in out["wavesets"].items():
                if not waves or not re.fullmatch(r"[\w.-]+", name):
                    sys.exit(f"{path}: wave set '{name}' needs a plain name and some waves")
            continue
        if section.lower() == "samples":
            out["sample"] = [[f"{tok}={localize(v)}" for tok, v in cp.items(section)]]
            continue
        for key, raw in cp.items(section):
            dest = key.strip().lower().replace("-", "_")
            if dest == "r":
                dest = "rand"
            if dest not in known:
                sys.exit(f"{path}: unknown setting '{key}' in [{section}]")
            if dest not in actions or raw.strip() == "":
                continue
            a = actions[dest]
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
        subs[known.step].set_defaults(
            **config_defaults(known.config, subs[known.step], subs))
    args = parser.parse_args()
    steps = {"bin": [step_bin], "fit": [step_fit], "gather": [step_gather],
             "plot": [step_plot], "plotfit": [step_plotfit],
             "all": [step_bin, step_fit, step_gather, step_plot]
                    + ([step_plotfit] if getattr(args, "plotfit", False) else [])}[args.step]
    if step_bin in steps and (not args.sample or not args.mass):
        sys.exit("bin needs --sample and --mass")
    if step_bin in steps:
        for name, rng in (("mass", args.mass), ("t", args.t)):
            if rng and (rng[1] <= rng[0] or rng[2] < 1):
                sys.exit(f"{name} = {rng[0]:g} {rng[1]:g} {rng[2]:g}: need LOW HIGH NBINS "
                         f"with HIGH > LOW and NBINS >= 1")
    if step_bin in steps and args.t_log and (not args.t or args.t[0] <= 0):
        sys.exit("--t-log needs --t with LO > 0")
    if step_fit in steps and not args.cfg:
        sys.exit("fit needs --cfg")
    wavesets = args.wavesets or {}
    if args.waveset:
        unknown = [n for n in args.waveset if n not in wavesets]
        if unknown:
            sys.exit(f"unknown wave set(s) {' '.join(unknown)}; "
                     f"defined: {' '.join(wavesets) or 'none'}")
        wavesets = {n: wavesets[n] for n in args.waveset}
    if not wavesets:
        for step in steps:
            step(args)
        return

    base = Path(args.workdir)
    if args.cfg and Path(args.cfg).exists():
        known = template_waves(Path(args.cfg).read_text())
        wavesets = {n: expand_waves(w, known) for n, w in wavesets.items()}
    if step_bin in steps:
        step_bin(args)
    for name, waves in wavesets.items():
        print(f"\n=== wave set {name}: {' '.join(waves)} ===")
        sub = argparse.Namespace(**vars(args))
        sub.workdir, sub.waves = str(base / name), waves
        if step_fit in steps:
            link_waveset_bins(base, base / name)
        for step in steps:
            if step is not step_bin:
                step(sub)
    if (step_gather in steps or step_plot in steps) and not args.dry_run:
        combined = combine_waveset_results(base, wavesets)
        if combined and step_plot in steps:
            plot_wavesets(combined, base / "plots_wavesets", args)


def link_waveset_bins(base, setdir):
    """Give a wave set its own bin directories, linking the shared ROOT files."""
    bins = load_bins(base)
    setdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(base / BINS_JSON, setdir / BINS_JSON)
    for b in bins:
        d = setdir / b["name"]
        d.mkdir(exist_ok=True)
        for f in b["files"].values():
            link = d / f
            if not link.exists():
                os.symlink(os.path.relpath(base / b["name"] / f, d), link)


def combine_waveset_results(base, wavesets):
    rows = []
    for name in wavesets:
        path = base / name / RESULTS_CSV
        if path.exists():
            with open(path, newline="") as f:
                rows += [{"waveset": name, **r} for r in csv.DictReader(f)]
    if not rows:
        return None
    keys = list(dict.fromkeys(k for r in rows for k in r))
    out = base / "results_wavesets.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {out} ({len(rows)} rows from {len(wavesets)} wave sets)")
    return out


def plot_wavesets(path, outdir, args):
    """One figure per quantity, a panel per t bin, one curve per wave set."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    sets = list(dict.fromkeys(r["waveset"] for r in rows))
    quantities = [k for k in rows[0]
                  if k not in BIN_COLS | {"waveset"} and not k.endswith("_err")
                  and any(num(r.get(k)) is not None for r in rows)]
    if args.only:
        quantities = [q for q in quantities if q in args.only]
    outdir.mkdir(exist_ok=True)
    t_bins = sorted({int(r["i_t"]) for r in rows})
    # small horizontal offsets so error bars of different sets don't overlap
    widths = [float(r["mass_hi"]) - float(r["mass_lo"]) for r in rows]
    step = 0.6 * min(widths) / max(len(sets), 1) if widths else 0.0
    for q in quantities:
        fig, axes = plt.subplots(1, len(t_bins), figsize=(6 * len(t_bins), 4.5),
                                 squeeze=False, sharey=True)
        for ax, j in zip(axes[0], t_bins):
            for k, name in enumerate(sets):
                pts = [r for r in rows if r["waveset"] == name and int(r["i_t"]) == j
                       and num(r.get(q)) is not None]
                if not pts:
                    continue
                shift = (k - (len(sets) - 1) / 2) * step
                ax.errorbar([float(r["mass"]) + shift for r in pts],
                            [num(r[q]) for r in pts],
                            yerr=[num(r.get(q + "_err")) or 0.0 for r in pts],
                            fmt="o", ms=4, capsize=2, label=name)
            r0 = next(r for r in rows if int(r["i_t"]) == j)
            if r0["t_lo"]:
                ax.set_title(f"{float(r0['t_lo']):.2f} < -t < {float(r0['t_hi']):.2f}")
            ax.set_xlabel(args.mass_label)
        axes[0][0].set_ylabel(q)
        axes[0][-1].legend(frameon=False, title="wave set")
        fig.tight_layout()
        fig.savefig(outdir / f"{q}.{args.format}")
        plt.close(fig)
    print(f"Wrote {len(quantities)} comparison plots to {outdir}/")


if __name__ == "__main__":
    main()
