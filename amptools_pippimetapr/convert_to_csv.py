"""Convert eta' pi- (gamma p -> eta' pi- pi+ p) AmpTools .fit files, or their ROOT data files, into a csv.

1. .fit files -> one row per fit with likelihood, fit status, and for each wave in
   fit_etaprpim_SPD_TEMPLATE.cfg its intensity (+error) and production parameter
   (Re, Im), the reflectivity / total intensities and selected phase differences.
   Optionally also the covariance / correlation matrices (long format).
   Runs the compiled helper extract_fit_results_etapi (build it like
   plot_etapi_delta; see the top of extract_fit_results_etapi.cc).
2. .root files (ROOTDataReader format) -> one row per file with the number of
   events and the eta' pi- mass mean/rms/min/max. Runs extract_bin_info_etapi.C.

If the paths contain bin_<N> and --mass-range is given (default 0.8 2.0 20, same
as plot_etapi_delta), 'bin' and 'mass' columns are added.

Examples:
  python3 convert_to_csv.py -i "EtaPi_fit/bin_*/bin_*.fit" -o fits.csv
  python3 convert_to_csv.py -i EtaPi_fit/bin_*/bin_*.fit -a --phases P1+:D1+ --covariance
  python3 convert_to_csv.py -i "EtaPi_fit/bin_*/data_*.root" -o data.csv
"""

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FIT_EXE = os.path.join(SCRIPT_DIR, "extract_fit_results_etapi")
ROOT_MACRO = os.path.join(SCRIPT_DIR, "extract_bin_info_etapi.C")


def sort_input_files(files: list, sort_index: int = -1) -> list:
    """Sort files by the number at position sort_index among all numbers in the path
    (default: the last one), so bin_10 comes after bin_9."""

    def key(path):
        nums = re.findall(r"\d+", path)
        try:
            return (0, int(nums[sort_index]), path)
        except IndexError:
            return (1, 0, path)

    return sorted(files, key=key)


def expand_inputs(inputs: list) -> list:
    # a single non-.fit/.root file is treated as a list of files
    if (
        len(inputs) == 1
        and os.path.isfile(inputs[0])
        and not inputs[0].endswith((".fit", ".root"))
    ):
        with open(inputs[0]) as f:
            inputs = [line.strip() for line in f if line.strip()]

    files = []
    for item in inputs:
        if any(c in item for c in "*?["):
            matches = glob.glob(item)
            if not matches:
                raise FileNotFoundError(f"No files match {item}")
            files.extend(matches)
        else:
            if not os.path.exists(item):
                raise FileNotFoundError(f"The file {item} does not exist")
            files.append(item)
    return [os.path.abspath(f) for f in files]


def add_bin_mass_columns(csv_path: str, mass_range) -> None:
    """Insert 'bin' and 'mass' (bin center) columns after 'file', using bin_<N> in the path."""
    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))
    if not rows:
        return
    header, body = rows[0], rows[1:]
    low, high, nbins = mass_range
    step = (high - low) / nbins

    new_rows = [[header[0], "bin", "mass"] + header[1:]]
    for row in body:
        m = re.findall(r"bin_(\d+)", row[0])
        if m:
            b = int(m[-1])
            new_rows.append([row[0], str(b), f"{low + step * (b + 0.5):.6g}"] + row[1:])
        else:
            new_rows.append([row[0], "", ""] + row[1:])
    with open(csv_path, "w", newline="") as f:
        csv.writer(f).writerows(new_rows)


def run(cmd: list, verbose: bool) -> None:
    print("Running:", " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if verbose and proc.stdout:
        print(proc.stdout, end="")
    if proc.returncode != 0:
        print(proc.stdout, end="")
        print(proc.stderr, end="", file=sys.stderr)
        sys.exit(f"Command failed with return code {proc.returncode}")
    elif proc.stderr:  # warnings
        print(proc.stderr, end="", file=sys.stderr)


def main(args: dict) -> None:
    if args["output"] and not args["output"].endswith(".csv"):
        args["output"] += ".csv"

    input_files = expand_inputs(args["input"])

    if all(f.endswith(".fit") for f in input_files):
        file_type = "fit"
    elif all(f.endswith(".root") for f in input_files):
        file_type = "root"
    else:
        raise ValueError("All input files must be either .fit or .root files")

    if args["sorted"]:
        input_files = sort_input_files(input_files, args["sort_index"])

    if args["preview"]:
        print("Files that will be processed:")
        for f in input_files:
            print(f"\t{f}")
        return

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".txt") as tmp:
        tmp.write("\n".join(input_files) + "\n")
        list_path = tmp.name

    try:
        if file_type == "fit":
            output = args["output"] or "fits.csv"
            if not os.path.isfile(FIT_EXE):
                sys.exit(
                    f"{FIT_EXE} not found. Build it first, e.g.:\n"
                    "  g++ -O2 extract_fit_results_etapi.cc -I$AMPTOOLS"
                    " -L$AMPTOOLS/lib -lAmpTools $(root-config --cflags --libs)"
                    " -o extract_fit_results_etapi"
                )
            cmd = [
                FIT_EXE,
                list_path,
                output,
                "1" if args["acceptance_corrected"] else "0",
                args["phases"],
                output.replace(".csv", "_cov.csv") if args["covariance"] else "-",
                output.replace(".csv", "_corr.csv") if args["correlation"] else "-",
            ]
            run(cmd, args["verbose"])
        else:
            output = args["output"] or "data.csv"
            if not shutil.which("root"):
                sys.exit("'root' not found on PATH. Set up your ROOT / GlueX environment first.")
            macro_call = (
                f'{ROOT_MACRO}("{list_path}", "{output}",'
                f' "{args["tree_name"]}", "{args["indices"]}")'
            )
            run(["root", "-n", "-l", "-b", "-q", macro_call], args["verbose"])

        if args["mass_range"] and os.path.isfile(output):
            add_bin_mass_columns(output, args["mass_range"])
        print(f"Wrote {output}")
    finally:
        os.remove(list_path)


def parse_args() -> dict:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument(
        "-i", "--input", nargs="+", required=True,
        help="Input .fit or .root file(s); wildcards (quote them) or a text file listing files",
    )
    p.add_argument("-o", "--output", default="", help="Output csv (default fits.csv / data.csv)")
    p.add_argument(
        "--no-sort", dest="sorted", action="store_false",
        help="Keep input order instead of sorting by the number in the path",
    )
    p.add_argument(
        "--sort-index", type=int, default=-1,
        help="Which number in the path to sort by (default -1 = last)",
    )
    p.add_argument(
        "-a", "--acceptance-corrected", action="store_true",
        help="Acceptance-corrected intensities (default: reconstructed). Note:"
        " plot_etapi_delta uses AmpTools' default, which IS acceptance corrected",
    )
    p.add_argument(
        "--phases", default="auto",
        help="Phase differences to include, e.g. 'P1+:D1+,S0+:D1+'; 'auto' (default) ="
        " all pairs of non-zero waves within a reflectivity; 'none' to skip",
    )
    p.add_argument("--covariance", action="store_true", help="Also write <output>_cov.csv")
    p.add_argument("--correlation", action="store_true", help="Also write <output>_corr.csv")
    p.add_argument(
        "--mass-range", nargs=3, type=float, metavar=("LOW", "HIGH", "NBINS"),
        default=[0.8, 2.0, 20],
        help="Mass binning used to add bin/mass columns from bin_<N> in the path"
        " (default 0.8 2.0 20)",
    )
    p.add_argument(
        "--no-mass", dest="mass_range", action="store_const", const=None,
        help="Don't add bin/mass columns",
    )
    p.add_argument("-t", "--tree-name", default="kin", help="Tree name in ROOT files (default kin)")
    p.add_argument(
        "--indices", default="1,2",
        help="Final-state indices for the mass in ROOT files (default 1,2 = EtaPrime, Pi-; 0,3 = p pi+)",
    )
    p.add_argument("-p", "--preview", action="store_true", help="List files to process and exit")
    p.add_argument("-v", "--verbose", action="store_true", help="Print helper output")
    args = vars(p.parse_args())
    if args["mass_range"]:
        args["mass_range"][2] = int(args["mass_range"][2])
    return args


if __name__ == "__main__":
    main(parse_args())
