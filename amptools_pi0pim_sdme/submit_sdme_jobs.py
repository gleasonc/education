#!/usr/bin/env python3
"""
Run the SDME fits on the JLab farm: one swif2 job per bin.

Each job runs, for its bin, the fits (`amptools_sdme.py fit --bin ...`) and the
fit-quality plots (`amptools_sdme.py plotfit --bin ...`).  Binning, gathering
and the result plots are quick and run on ifarm itself:

  python3 amptools_sdme.py bin -c my_sdme.ini        # once, on ifarm
  python3 submit_sdme_jobs.py -c my_sdme.ini         # one job per bin
  swif2 status -workflow <workflow>                  # wait until all jobs are done
  python3 amptools_sdme.py gather  -c my_sdme.ini
  python3 amptools_sdme.py plot    -c my_sdme.ini
  python3 amptools_sdme.py plotfit -c my_sdme.ini    # only merges the per-bin pages

The job settings are in the [farm] section of the INI file (see
example_sdme.ini); options on the command line override them.  Run it from the
same folder as the other steps, since `workdir` is relative to it.
--dry-run writes the job scripts and prints the swif2 commands without
submitting anything.
"""

import argparse
import configparser
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

FARM_DEFAULTS = {
    "workflow": "sdme_pi0pim",
    "account": "halld",
    "partition": "production",
    "os": "el9",
    "cores": "5",               # the 5 fits of a bin (all, pol00 ... pol135) run in parallel
    "ram": "4GB",
    "disk": "5GB",
    "time": "8hrs",
    "env_setup": "source /group/halld/Software/build_scripts/gluex_env_boot_jlab.csh\n"
                 "gxenv $HALLD_VERSIONS/version.xml",
}

JOB_TEMPLATE = """#!/bin/tcsh -f
# SDME fits of {bin} (written by submit_sdme_jobs.py)
{env}
cd {launch}
python3 {driver} fit -c {ini} --workdir {workdir} --jobs {cores} --bin {bin}{fit_extra}
if ( $status != 0 ) exit $status
python3 {driver} plotfit -c {ini} --workdir {workdir} --jobs {cores} --bin {bin}
"""


def run(cmd, dry):
    print("$ " + cmd)
    if dry:
        return 0
    return subprocess.call(cmd, shell=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--config", required=True, help="the INI file used for the other steps")
    ap.add_argument("--workdir", help="override [general] workdir")
    ap.add_argument("--bin", nargs="+", metavar="NAME", help="submit only these bins")
    for key in FARM_DEFAULTS:
        if key != "env_setup":
            ap.add_argument("--" + key.replace("_", "-"), help=f"[farm] {key}")
    ap.add_argument("--fit-args", default="",
                    help='extra amptools_sdme.py fit options for the jobs, e.g. "-r 10"')
    ap.add_argument("--no-run", action="store_true",
                    help="add the jobs but do not start the workflow (swif2 run)")
    ap.add_argument("--dry-run", action="store_true",
                    help="write the job scripts, print the swif2 commands, submit nothing")
    args = ap.parse_args()

    ini = Path(args.config).resolve()
    cp = configparser.ConfigParser(inline_comment_prefixes=("#",), interpolation=None)
    cp.read(ini)
    farm = dict(FARM_DEFAULTS)
    if cp.has_section("farm"):
        farm.update({k: v.strip() for k, v in cp.items("farm") if v.strip()})
    for key in FARM_DEFAULTS:
        v = getattr(args, key, None)
        if v:
            farm[key] = v

    workdir = args.workdir or (cp.get("general", "workdir", fallback="").strip() or "SDME_fit")
    workdir = Path(workdir).resolve()
    bins_json = workdir / "bins.json"
    if not bins_json.exists():
        sys.exit(f"{bins_json} not found; run `amptools_sdme.py bin -c {args.config}` first")
    bins = [b["name"] for b in json.loads(bins_json.read_text())]
    if args.bin:
        unknown = [b for b in args.bin if b not in bins]
        if unknown:
            sys.exit(f"--bin: no bin {' '.join(unknown)} in {bins_json}")
        bins = [b for b in bins if b in args.bin]

    jobdir = workdir / "jobs"
    (jobdir / "log").mkdir(parents=True, exist_ok=True)
    wf = farm["workflow"]
    fit_extra = (" " + args.fit_args) if args.fit_args else ""

    if run(f"swif2 create -workflow {shlex.quote(wf)}", args.dry_run):
        print(f"(workflow {wf} probably exists already; adding the jobs to it)")
    failed = []
    for b in bins:
        script = jobdir / f"{b}.csh"
        script.write_text(JOB_TEMPLATE.format(
            bin=b, env=farm["env_setup"], launch=shlex.quote(os.getcwd()),
            driver=shlex.quote(str(HERE / "amptools_sdme.py")), ini=shlex.quote(str(ini)),
            workdir=shlex.quote(str(workdir)), cores=farm["cores"], fit_extra=fit_extra))
        script.chmod(0o755)
        cmd = (f"swif2 add-job -workflow {shlex.quote(wf)} -name {shlex.quote(wf + '_' + b)}"
               f" -account {farm['account']} -partition {farm['partition']} -os {farm['os']}"
               f" -cores {farm['cores']} -disk {farm['disk']} -ram {farm['ram']}"
               f" -time {farm['time']}"
               f" -stdout {jobdir / 'log' / (b + '.out')} -stderr {jobdir / 'log' / (b + '.err')}"
               f" {script}")
        if run(cmd, args.dry_run):
            failed.append(b)
    if failed:
        sys.exit(f"swif2 add-job failed for {', '.join(failed)}")
    if not args.no_run:
        run(f"swif2 run -workflow {shlex.quote(wf)}", args.dry_run)
    print(f"{len(bins)} jobs {'written' if args.dry_run else 'submitted'} "
          f"(scripts in {jobdir}). Check with: swif2 status -workflow {wf}")


if __name__ == "__main__":
    main()
