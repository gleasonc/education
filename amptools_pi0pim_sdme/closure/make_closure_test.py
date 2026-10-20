#!/usr/bin/env python3
"""
Closure test: generate pi- pi0 events with known SDMEs and fit them back.

For every -t bin of your INI and every polarization orientation this writes a
gen_amp config whose TwoPiAngles amplitude has the SDMEs of a reference table
(by default the GlueX preliminary values in reference_sdmes_gluex_prelim.csv,
interpolated linearly to the bin centre). It also writes

  generate.sh        runs gen_amp for every config, plus flat (-f) MC per t bin,
                     and merges everything into data.root and flat.root in the
                     format amptools_sdme.py reads (tree kin, Proton Pi- Pi0,
                     Weight/F, pol/I)
  truth.csv          the SDMEs that were generated, per t bin
  closure_sdme.ini   your INI with the samples pointed at data.root / flat.root,
                     workdir = closure_fit and the plots overlaid with truth.csv

so the test is:

  python3 closure/make_closure_test.py -c my_sdme.ini --outdir closure_test
  bash closure_test/generate.sh            # needs gen_amp and root (GlueX env)
  python3 amptools_sdme.py all -c closure_test/closure_sdme.ini

The fitted SDMEs should agree with truth.csv (gray boxes in the plots).
This is a generator-level test (no detector): it checks the fit set-up, the
binning and the polarization bookkeeping. With --deut Neutron (and --fermi
FILE) the events are generated on a moving neutron in the deuteron, with a
gen_amp that has the -deut/-fermi options; the ROOT files then have the same
kinematics as gamma d -> pi- pi0 p (p_spectator) data. To include the
acceptance, pass the generated files (HDDM, with the spectator) through your
simulation and DSelector before fitting.
"""

import argparse
import configparser
import csv
import math
import os
import re
import shlex
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SDME_NAMES = ["rho000", "rho100", "rho1m10", "rho111", "rho001",
              "rho101", "rho1m11", "rho102", "rho1m12"]

CFG_TEMPLATE = """# gen_amp config written by make_closure_test.py
# -t bin {t_lo:.4f}-{t_hi:.4f} GeV^2, orientation {suf} (polAngle {ang} deg, P = {pval})
# SDMEs: {sdme_comment}
define rho {rho_mass} {rho_width}
define polAngle {ang}
define polFrac {pval}

fit closure_{tag}

reaction Pi0Pim Beam Proton Pi- Pi0

sum Pi0Pim xpol

amplitude Pi0Pim::xpol::rhoS TwoPiAngles {sdmes} polAngle polFrac
amplitude Pi0Pim::xpol::rhoS BreitWigner rho 1 2 3

initialize Pi0Pim::xpol::rhoS cartesian 1.0 0.0
"""

FLAT_CFG = """# gen_amp config for flat MC (used with -f), written by make_closure_test.py
define rho {rho_mass} {rho_width}
fit closure_flat
reaction Pi0Pim Beam Proton Pi- Pi0
sum Pi0Pim xpol
amplitude Pi0Pim::xpol::rhoS BreitWigner rho 1 2 3
initialize Pi0Pim::xpol::rhoS cartesian 1.0 0.0
"""


def read_reference(path):
    rows = []
    with open(path) as f:
        lines = [l for l in f if l.strip() and not l.lstrip().startswith("#")]
    for r in csv.DictReader(lines):
        rows.append({k.strip(): float(v) for k, v in r.items() if v not in (None, "")})
    missing = [n for n in SDME_NAMES if n not in rows[0]]
    if "t" not in rows[0] or missing:
        sys.exit(f"{path}: needs columns t and {' '.join(SDME_NAMES)}")
    return sorted(rows, key=lambda r: r["t"])


def interpolate(rows, t):
    """Linear interpolation in -t; constant beyond the first/last point."""
    if t <= rows[0]["t"]:
        return {n: rows[0][n] for n in SDME_NAMES}
    if t >= rows[-1]["t"]:
        return {n: rows[-1][n] for n in SDME_NAMES}
    for a, b in zip(rows, rows[1:]):
        if a["t"] <= t <= b["t"]:
            f = (t - a["t"]) / (b["t"] - a["t"])
            return {n: a[n] + f * (b[n] - a[n]) for n in SDME_NAMES}


def min_w(s, pol, n=24):
    """Smallest TwoPiAngles W(theta, phi, Phi) (without 3/4pi) on a grid.
    gen_amp uses |W|, so W < 0 somewhere means the sample does not follow
    these SDMEs there."""
    lo = float("inf")
    for i in range(n + 1):
        ct = -1 + 2 * i / n
        st2 = 1 - ct * ct
        s2t = 2 * ct * math.sqrt(max(st2, 0.0))
        for j in range(n):
            ph = 2 * math.pi * j / n
            w0 = (0.5 * (1 - s["rho000"]) + 0.5 * (3 * s["rho000"] - 1) * ct * ct
                  - math.sqrt(2) * s["rho100"] * s2t * math.cos(ph)
                  - s["rho1m10"] * st2 * math.cos(2 * ph))
            w1 = (s["rho111"] * st2 + s["rho001"] * ct * ct
                  - math.sqrt(2) * s["rho101"] * s2t * math.cos(ph)
                  - s["rho1m11"] * st2 * math.cos(2 * ph))
            w2 = math.sqrt(2) * s["rho102"] * s2t * math.sin(ph) + s["rho1m12"] * st2 * math.sin(2 * ph)
            # minimum over Phi of w0 - P (w1 cos2Phi + w2 sin2Phi)
            lo = min(lo, w0 - pol * math.hypot(w1, w2))
    return lo


def edges(lo, hi, n, log):
    n = int(n)
    if log:
        return [lo * (hi / lo) ** (i / n) for i in range(n + 1)]
    return [lo + (hi - lo) * i / n for i in range(n + 1)]


def floats(v):
    return [float(x) for x in re.split(r"[\s,]+", v.strip()) if x]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True,
                    help="your amptools_sdme INI (t binning, mass range, orientations)")
    ap.add_argument("--sdmes", default=str(HERE / "reference_sdmes_gluex_prelim.csv"),
                    help="CSV with columns t,rho000,...,rho1m12 (default: %(default)s)")
    ap.add_argument("--outdir", default="closure_test")
    ap.add_argument("--nevents", type=int, default=20000,
                    help="events per t bin and orientation (default %(default)s)")
    ap.add_argument("--nflat", type=int, default=200000,
                    help="flat MC events per t bin (default %(default)s)")
    ap.add_argument("--pol", nargs="+", metavar="SUF:ANGLE:P",
                    default=["00:0:0.3519", "45:45:0.3374", "90:90:0.3303", "135:135:0.3375"],
                    help="orientations: token suffix, polAngle (deg), polarization "
                         "(default: the values of the fit templates)")
    ap.add_argument("--ebeam", nargs=2, type=float, default=[8.2, 8.8], metavar=("LO", "HI"),
                    help="beam energy range for gen_amp -a/-b (default 8.2 8.8)")
    ap.add_argument("--slope", type=float, default=3.0,
                    help="gen_amp -t: exponential t slope inside each bin (default 3)")
    ap.add_argument("--rho", nargs=2, type=float, default=[0.7751, 0.1491],
                    metavar=("MASS", "WIDTH"), help="rho- Breit-Wigner (default 0.7751 0.1491)")
    ap.add_argument("--gen-amp", default="gen_amp", help="gen_amp executable")
    ap.add_argument("--root-exe", default="root")
    ap.add_argument("--seed", type=int, default=1000, help="first random seed")
    ap.add_argument("--deut", choices=["Neutron", "Proton"],
                    help="deuteron target (needs a gen_amp with -deut/-fermi): the nucleon the "
                         "photon hits, e.g. Neutron for gamma d -> pi- pi0 p (p)")
    ap.add_argument("--fermi", help="gen_amp -fermi: |p_Fermi| distribution file "
                                    "(default: gen_amp's Hulthen wave function)")
    ap.add_argument("--t-margin", type=float, default=0.03,
                    help="with --deut: generate each bin over [t_lo - M, t_hi + M] and keep "
                         "only events whose -t falls in the bin (default %(default)s)")
    ap.add_argument("--jobs", type=int, default=1, help="gen_amp jobs run in parallel by generate.sh")
    args = ap.parse_args()

    ini = Path(args.config).resolve()
    cp = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None)
    cp.read(ini)
    mass = floats(cp.get("bin", "mass", fallback="0.6 0.95 1"))
    tdef = cp.get("bin", "t", fallback="").strip()
    if not tdef:
        sys.exit("the INI has no t binning ([bin] t = LO HI N); the closure test needs t bins")
    tlo, thi, tn = floats(tdef)
    tlog = cp.get("bin", "t_log", fallback="false").strip().lower() in ("1", "true", "yes")
    tedges = edges(tlo, thi, tn, tlog)
    pols = []
    for p in args.pol:
        suf, ang, val = p.split(":")
        pols.append((suf, ang, val))
    pol_values = {}
    for item in cp.get("polarization", "orientations", fallback="00:0 45:45 90:90 135:135").split():
        suf, val = item.split(":")
        pol_values[suf] = int(val)
    unknown = [s for s, _, _ in pols if s not in pol_values]
    if unknown:
        sys.exit(f"--pol suffixes {unknown} are not in [polarization] orientations of {ini}")
    pol_branch = cp.get("polarization", "pol_branch", fallback="pol").strip() or "pol"

    ref = read_reference(args.sdmes)
    out = Path(args.outdir).resolve()
    (out / "cfg").mkdir(parents=True, exist_ok=True)
    (out / "raw").mkdir(exist_ok=True)

    truth, cmds, data_files, flat_files = [], [], [], []
    seed = args.seed
    warn = []
    (out / "cfg" / "gen_flat.cfg").write_text(FLAT_CFG.format(rho_mass=args.rho[0],
                                                              rho_width=args.rho[1]))
    common = (f"-l {mass[0]} -u {mass[1]} -a {args.ebeam[0]} -b {args.ebeam[1]} "
              f"-t {args.slope}")
    if args.deut:
        common += f" -deut {args.deut}"
        if args.fermi:
            common += f" -fermi {Path(args.fermi).resolve()}"
    # Fermi motion moves -t a little: generate a bit wider, cut back to the bin in the merge
    margin = args.t_margin if args.deut else 0.0

    def scale(n, a, b):
        """Events to generate so that about n remain in [a, b] after the margin cut."""
        return int(round(n * (b - a + 2 * margin) / (b - a)))
    for k, (a, b) in enumerate(zip(tedges, tedges[1:])):
        s = interpolate(ref, 0.5 * (a + b))
        truth.append({"t_lo": a, "t_hi": b, "t": 0.5 * (a + b), **s})
        trange = f"-tmin {max(a - margin, 0.0):.5f} -tmax {b + margin:.5f}"
        for suf, ang, val in pols:
            tag = f"t{k:02d}_pol{suf}"
            w = min_w(s, float(val))
            if w < 0:
                warn.append(f"  {tag}: W < 0 somewhere (min {w:.3f}); gen_amp uses |W|, "
                            "so this bin will not close exactly")
            cfg = out / "cfg" / f"gen_{tag}.cfg"
            cfg.write_text(CFG_TEMPLATE.format(
                t_lo=a, t_hi=b, suf=suf, ang=ang, pval=val, tag=tag,
                rho_mass=args.rho[0], rho_width=args.rho[1],
                sdmes=" ".join(f"{s[n]:.4f}" for n in SDME_NAMES),
                sdme_comment=" ".join(f"{n}={s[n]:.3f}" for n in SDME_NAMES)))
            raw = out / "raw" / f"{tag}.root"
            cmds.append(f"run {shlex.quote(str(cfg))} {shlex.quote(str(raw))} "
                        f"\"{common} {trange} -n {scale(args.nevents, a, b)} -s {seed}\"")
            data_files.append((raw, pol_values[suf], a, b))
            seed += 1
        raw = out / "raw" / f"t{k:02d}_flat.root"
        cmds.append(f"run {shlex.quote(str(out / 'cfg' / 'gen_flat.cfg'))} {shlex.quote(str(raw))} "
                    f"\"{common} {trange} -n {scale(args.nflat, a, b)} -s {seed} -f\"")
        flat_files.append((raw, -1, a, b))
        seed += 1

    with open(out / "truth.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["t_lo", "t_hi", "t"] + SDME_NAMES)
        w.writeheader()
        for r in truth:
            w.writerow({k: (f"{v:.5f}" if isinstance(v, float) else v) for k, v in r.items()})

    # file lists for the merge macro: "file pol tmin tmax" per line
    for name, files in (("data_files.txt", data_files), ("flat_files.txt", flat_files)):
        (out / name).write_text("".join(f"{p} {v} {a:.6f} {b:.6f}\n" for p, v, a, b in files))
    merge = HERE / "merge_closure.C"
    script = f"""#!/bin/bash
# Generate the closure-test samples (written by make_closure_test.py).
# Needs gen_amp and root in the path (the GlueX environment).
set -e
cd {shlex.quote(str(out / 'raw'))}
GEN_AMP=${{GEN_AMP:-{args.gen_amp}}}
JOBS=${{JOBS:-{args.jobs}}}
FAILED=$PWD/failed.txt; rm -f "$FAILED"
run() {{  # cfg output "options"; each job in its own folder (gen_amp writes local files)
  if [ -s "$2" ]; then echo "exists: $2"; return; fi
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 1; done
  echo "gen_amp -c $1 -o $2 $3"
  (d="${{2%.root}}.d"; mkdir -p "$d" && cd "$d" &&
   $GEN_AMP -c "$1" -o "$2" $3 > "${{2%.root}}.log" 2>&1 || echo "$2" >> "$FAILED") &
}}
{chr(10).join(cmds)}
wait
if [ -s "$FAILED" ]; then echo "gen_amp failed for:"; cat "$FAILED"; exit 1; fi
{args.root_exe} -l -b -q '{merge}("{out / 'data.root'}","{out / 'data_files.txt'}","{pol_branch}")'
{args.root_exe} -l -b -q '{merge}("{out / 'flat.root'}","{out / 'flat_files.txt'}","{pol_branch}")'
echo "Done: {out / 'data.root'} and {out / 'flat.root'}"
"""
    (out / "generate.sh").write_text(script)
    (out / "generate.sh").chmod(0o755)

    # INI for the fit of the generated samples
    cp2 = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None)
    cp2.optionxform = str
    cp2.read(ini)
    here_main = HERE.parent
    for sec in ("samples", "general", "plot", "fit", "bin", "plotfit"):
        if not cp2.has_section(sec):
            cp2.add_section(sec)
    for key in list(cp2["samples"].keys()):
        cp2.remove_option("samples", key)
    cp2.set("samples", "DATAFILE", str(out / "data.root"))
    cp2.set("samples", "ACCMCFILE", str(out / "flat.root"))
    cp2.set("samples", "GENMCFILE", str(out / "flat.root"))
    cp2.set("polarization", "split", "DATAFILE")
    cp2.set("general", "workdir", str(out / "closure_fit"))
    cp2.set("plot", "reference", str(out / "truth.csv"))
    cp2.set("plot", "reference_label", "generated")
    # templates and splitter: absolute paths, so the INI works from anywhere
    for sec, key in (("fit", "cfg"), ("fit", "cfg_single"), ("bin", "splitter"),
                     ("plotfit", "plotter")):
        v = cp2.get(sec, key, fallback="").strip()
        if v and not os.path.isabs(v):
            for base in (ini.parent, here_main):
                if (base / v).exists():
                    cp2.set(sec, key, str(base / v))
                    break
    if cp2.has_section("farm"):
        cp2.remove_section("farm")
    with open(out / "closure_sdme.ini", "w") as f:
        f.write("# Closure test: written by make_closure_test.py from "
                f"{ini}\n# truth.csv has the generated SDMEs\n")
        cp2.write(f)

    print(f"Wrote {len(data_files)} gen_amp configs + flat MC for {len(tedges) - 1} t bins in {out}")
    print(f"  truth:   {out / 'truth.csv'}")
    print(f"  next:    bash {out / 'generate.sh'}")
    print(f"           python3 {here_main / 'amptools_sdme.py'} all -c {out / 'closure_sdme.ini'}")
    if warn:
        print("Warning:\n" + "\n".join(warn))


if __name__ == "__main__":
    main()
