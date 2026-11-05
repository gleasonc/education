#!/usr/bin/env python3
"""
Pulls of a closure test: (fitted - generated) / error for every SDME.

  python3 closure/closure_pulls.py closure_test/closure_fit/results.csv closure_test/truth.csv

Writes closure_pulls.pdf next to results.csv:
  page 1  pull vs -t, a panel per SDME, one series per fit (all, pol00 ... pol135)
  page 2  pull histogram per fit, with mean and RMS, and a unit Gaussian
"""

import argparse
import csv
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from amptools_sdme import SDMES, FIT_STYLE, COMBINED  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="results.csv of the closure fit")
    ap.add_argument("truth", help="truth.csv written by make_closure_test.py")
    ap.add_argument("-o", "--output", help="output PDF (default: closure_pulls.pdf next to results.csv)")
    args = ap.parse_args()

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    truth = list(csv.DictReader(open(args.truth)))
    rows = [r for r in csv.DictReader(open(args.results)) if r.get("fit_ok") == "1"]

    def gen(t):
        return min(truth, key=lambda r: abs(float(r["t"]) - t))

    fits = [COMBINED] + sorted({r["fit"] for r in rows} - {COMBINED},
                               key=lambda s: int("".join(c for c in s if c.isdigit()) or 0))
    pulls = {f: {n: [] for n, _ in SDMES} for f in fits}     # (t, pull)
    for r in rows:
        t = float(r["t"])
        g = gen(t)
        for n, _ in SDMES:
            err = float(r[n + "_err"])
            if err > 0:
                pulls[r["fit"]][n].append((t, (float(r[n]) - float(g[n])) / err))

    out = Path(args.output) if args.output else Path(args.results).with_name("closure_pulls.pdf")
    with PdfPages(out) as pdf:
        fig, axes = plt.subplots(3, 3, figsize=(11, 10), sharex=True)
        nf = len(fits)
        for ax, (n, label) in zip(axes.flat, SDMES):
            for i, f in enumerate(fits):
                pts = pulls[f][n]
                if not pts:
                    continue
                color, name = FIT_STYLE.get(f, (None, f))
                dx = (i - (nf - 1) / 2) * 0.012            # spread the series a little
                ax.plot([t + dx for t, _ in pts], [p for _, p in pts], "o", ms=3.5,
                        color=color, label=name)
            ax.axhspan(-1, 1, color="gray", alpha=0.15, lw=0)
            ax.axhline(0, color="gray", lw=0.8)
            lim = max(4.0, max((abs(p) for f in fits for _, p in pulls[f][n]), default=0) * 1.1)
            ax.set_ylim(-lim, lim)
            ax.text(0.05, 0.92, label, transform=ax.transAxes, fontsize=13, va="top",
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.8, pad=1))
            ax.set_ylabel("pull")
            ax.tick_params(direction="in", top=True, right=True)
        for ax in axes[-1]:
            ax.set_xlabel(r"$-t$ [GeV$^2$]")
        h, l = axes.flat[0].get_legend_handles_labels()
        fig.legend(h, l, loc="upper center", ncol=len(l), frameon=False, fontsize=11)
        fig.suptitle("(fit $-$ generated) / error; band = $\\pm 1$", y=0.005, va="bottom", fontsize=10)
        fig.tight_layout(rect=(0, 0.02, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)

        fig, axes = plt.subplots(1, nf, figsize=(3.3 * nf, 3.4), squeeze=False)
        bins = [-5 + 0.5 * k for k in range(21)]
        for ax, f in zip(axes[0], fits):
            vals = [p for n, _ in SDMES for _, p in pulls[f][n]]
            if not vals:
                ax.axis("off")
                continue
            color, name = FIT_STYLE.get(f, (None, f))
            clipped = [max(min(v, 4.999), -4.999) for v in vals]
            ax.hist(clipped, bins=bins, color=color, alpha=0.6)
            mean = sum(vals) / len(vals)
            rms = math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals))
            xs = [-5 + 0.05 * k for k in range(201)]
            norm = len(vals) * 0.5 / math.sqrt(2 * math.pi)
            ax.plot(xs, [norm * math.exp(-0.5 * x * x) for x in xs], color="black", lw=1)
            ax.set_title(f"{name}: mean {mean:.2f}, RMS {rms:.2f}", fontsize=10)
            ax.set_xlabel("pull (overflow in the outer bins)")
            ax.tick_params(direction="in", top=True, right=True)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
