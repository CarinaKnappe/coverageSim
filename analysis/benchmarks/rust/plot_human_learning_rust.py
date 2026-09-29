"""Render and summarize unchanged upstream RUST output for all staged simulations."""

import csv
import importlib
import json
import os
from pathlib import Path
import sys

import numpy as np


REPO = Path(__file__).resolve().parents[3]
ROOT = REPO / "analysis/benchmarks/rust/runs/2026-09-11_human-learning"
os.environ["MPLCONFIGDIR"] = str(ROOT / "matplotlib_cache")
sys.path.insert(0, str(REPO / "analysis/benchmarks/rust/upstream/RUST"))
sys.path.insert(0, str(REPO / "analysis/benchmarks/rust"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rust_plot_compat import compatible_metafootprint


SCENARIOS = [
    "baseline", "counts_only", "geometry_only", "codon_only", "ends_only",
    "counts_geometry", "counts_geometry_codon", "all_supported",
]
TITLES = {
    "baseline": "Baseline",
    "counts_only": "Learned counts",
    "geometry_only": "Learned geometry",
    "codon_only": "Learned codon profile",
    "ends_only": "Learned end profiles",
    "counts_geometry": "Counts + geometry",
    "counts_geometry_codon": "Counts + geometry + codons",
    "all_supported": "All supported learned inputs",
}
LEFT_YLIM = (-3.0, 2.0)
RIGHT_YLIM = (0.0, 1.0)


def numerical_path(scenario, stratum, mode):
    folder = ROOT / "simulations" / scenario / "rust" / stratum / mode
    paths = list(folder.glob(f"RUST_{mode}_file_*"))
    if len(paths) != 1:
        raise RuntimeError(f"Expected one numerical output in {folder}, found {paths}")
    return paths[0]


def read_profile(scenario, stratum, mode):
    path = numerical_path(scenario, stratum, mode)
    lines = path.read_text().splitlines()
    split = lines.index("")
    rows = list(csv.reader(lines[:split], skipinitialspace=True))
    positions = np.array([int(value) for value in rows[0][2:]])
    features = [row[0] for row in rows[1:]]
    expected = np.array([float(row[1]) for row in rows[1:]])
    observed = np.array([[float(value) for value in row[2:]] for row in rows[1:]])
    divergence = np.array([
        float(value) if value.strip() != "NA" else np.nan
        for value in next(csv.reader([lines[split + 1]], skipinitialspace=True))[2:]
    ])
    return {
        "path": path, "positions": positions, "features": features,
        "ratios": observed / expected[:, None], "divergence": divergence,
    }


def draw(ax, profile, mode, title):
    before = set(ax.figure.axes)
    module = importlib.import_module(f"RUST.{mode}")
    with profile["path"].open() as handle:
        compatible_metafootprint(module)(handle, ax)
    added = set(ax.figure.axes) - before
    if len(added) == 1:
        right = added.pop()
        np.testing.assert_allclose(
            right.lines[0].get_ydata(), profile["divergence"], atol=1e-12
        )
    elif len(added) == 0 and np.isnan(profile["divergence"]).any():
        # Upstream suppresses the entire divergence axis when any position is NA.
        # Draw the unchanged numerical values with gaps so sparse outputs remain visible.
        right = ax.twinx()
        right.plot(profile["divergence"], color="blue")
    else:
        raise RuntimeError("Unexpected upstream RUST divergence plot state")
    with np.errstate(divide="ignore"):
        log_ratios = np.log2(profile["ratios"])
    for line, expected in zip(ax.lines[:len(profile["features"])], log_ratios):
        np.testing.assert_allclose(line.get_ydata(), expected, atol=1e-12)
    ax.set_ylim(LEFT_YLIM)
    right.set_ylim(RIGHT_YLIM)
    np.testing.assert_allclose(ax.get_ylim(), LEFT_YLIM)
    np.testing.assert_allclose(right.get_ylim(), RIGHT_YLIM)
    ticks = np.linspace(*RIGHT_YLIM, 5)
    right.set_yticks(ticks)
    right.set_yticklabels([f"{value:.2f}" for value in ticks])
    right.set_ylabel("RUST divergence", color="blue")
    ax.set_title(title)
    unit = "codons" if mode == "codon" else "nt"
    ax.set_xlabel(f"Distance from assigned A-site [{unit}]")
    ax.set_ylabel(f"{mode.capitalize()} RUST ratio, log2(observed/expected)")
    ax.axhline(0, color="grey", linestyle=":", linewidth=.5)
    ax.tick_params(labelsize=8)
    right.tick_params(labelsize=8)

def summarize(data, manifest):
    output = []
    for scenario in SCENARIOS:
        for stratum in ("length28", "offset15_pooled"):
            nt = data[scenario, stratum, "nucleotide"]
            codon = data[scenario, stratum, "codon"]
            npos = nt["positions"]
            cpos = codon["positions"]
            validation = manifest["validations"][scenario]
            selected = validation["full"] if stratum == "length28" else validation["offset15"]
            accepted = (sum(
                count for key, count in selected["records_by_length_offset"].items()
                if key.startswith("28:")
            ) if stratum == "length28" else selected["records"])
            retained = (sum(
                count for key, count in selected["retained_by_length_offset"].items()
                if key.startswith("28:")
            ) if stratum == "length28" else sum(selected["retained_by_length_offset"].values()))
            output.append({
                "scenario": scenario,
                "stratum": stratum,
                "bam_records": selected["records"],
                "accepted_records": accepted,
                "retained_records": retained,
                "a_site_nt_divergence": float(nt["divergence"][npos == 0][0]),
                "five_window_nt_peak": float(np.nanmax(nt["divergence"][(npos >= -16) & (npos <= -11)])),
                "three_window_nt_peak": float(np.nanmax(nt["divergence"][(npos >= 9) & (npos <= 18)])),
                "a_site_codon_divergence": float(codon["divergence"][cpos == 0][0]),
                "max_nt_divergence_position": int(npos[np.nanargmax(nt["divergence"])]),
            })
    with (ROOT / "rust_signal_summary.tsv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(output)
    return output


def main():
    data = {
        (scenario, stratum, mode): read_profile(scenario, stratum, mode)
        for scenario in SCENARIOS
        for stratum in ("length28", "offset15_pooled")
        for mode in ("codon", "nucleotide")
    }
    manifest = json.loads((ROOT / "rust_run_manifest.json").read_text())
    summary = summarize(data, manifest)

    with PdfPages(ROOT / "human_learning_RUST_metafootprints.pdf") as pdf:
        for stratum in ("length28", "offset15_pooled"):
            label = ("28 nt, offset 15" if stratum == "length28"
                     else "all reads with actual offset 15")
            for start in (0, 4):
                scenarios = SCENARIOS[start:start + 4]
                fig, axes = plt.subplots(2, 2, figsize=(15, 10))
                for ax, scenario in zip(axes.flat, scenarios):
                    draw(ax, data[scenario, stratum, "codon"], "codon",
                         TITLES[scenario])
                fig.suptitle(f"Original RUST codon footprint — {label}", fontsize=16)
                fig.tight_layout(rect=(0, 0, 1, .96))
                pdf.savefig(fig)
                plt.close(fig)

    with PdfPages(ROOT / "human_learning_RUST_nucleotide_diagnostics.pdf") as pdf:
        for stratum in ("length28", "offset15_pooled"):
            label = ("28 nt, offset 15" if stratum == "length28"
                     else "all reads with actual offset 15")
            for start in (0, 4):
                scenarios = SCENARIOS[start:start + 4]
                fig, axes = plt.subplots(2, 2, figsize=(15, 10))
                for ax, scenario in zip(axes.flat, scenarios):
                    draw(ax, data[scenario, stratum, "nucleotide"], "nucleotide",
                         TITLES[scenario])
                fig.suptitle(f"RUST nucleotide diagnostic — {label}", fontsize=16)
                fig.tight_layout(rect=(0, 0, 1, .96))
                pdf.savefig(fig)
                plt.close(fig)

    for stratum in ("length28", "offset15_pooled"):
        for mode in ("codon", "nucleotide"):
            fig, axes = plt.subplots(4, 2, figsize=(15, 18))
            for ax, scenario in zip(axes.flat, SCENARIOS):
                draw(ax, data[scenario, stratum, mode], mode,
                     TITLES[scenario])
            label = "28nt" if stratum == "length28" else "offset15_pool"
            fig.tight_layout()
            fig.savefig(ROOT / f"RUST_{label}_{mode}_comparison.png", dpi=140)
            plt.close(fig)

    (ROOT / "rust_plot_validation.json").write_text(json.dumps({
        "profiles": 32,
        "plotted_curves_match_upstream_numerical_output": True,
        "upstream_numerical_files_modified": False,
        "axis_limits": {
            "log2_observed_expected": LEFT_YLIM,
            "divergence": RIGHT_YLIM,
        },
        "axis_limits_verified_on_all_panels": True,
        "primary_pdf_mode": "codon",
        "nucleotide_output_role": "separate diagnostic",
        "profiles_with_sparse_NA_divergence": int(sum(
            np.isnan(profile["divergence"]).any() for profile in data.values()
        )),
        "zero_observed_expected_ratios": int(sum(
            np.count_nonzero(profile["ratios"] == 0) for profile in data.values()
        )),
        "NA_divergence_plot_handling": (
            "Upstream omits the divergence axis when any value is NA; the comparison "
            "renderer plots unchanged finite values with gaps at NA positions."
        ),
        "summary_rows": len(summary),
    }, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
