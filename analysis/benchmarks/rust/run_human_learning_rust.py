"""Convert all staged human-learning simulations to upstream RUST inputs and run them."""

import csv
import json
import os
from collections import Counter
from pathlib import Path
import subprocess
import sys

import pysam


REPO = Path(__file__).resolve().parents[3]
ROOT = REPO / "analysis/benchmarks/rust/runs/2026-09-11_human-learning"
UPSTREAM = REPO / "analysis/benchmarks/rust/upstream/RUST"
SCENARIOS = [
    "baseline", "counts_only", "geometry_only", "codon_only", "ends_only",
    "counts_geometry", "counts_geometry_codon", "all_supported",
]


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_transcript_bam(scenario, rows, annotation, offset15_only=False):
    name = "offset15.bam" if offset15_only else "transcript_fragments.bam"
    output_path = ROOT / "simulations" / scenario / name
    selected = [r for r in rows if not offset15_only or int(r["site_offset"]) == 15]
    refs = {r["transcript_id"]: i for i, r in enumerate(annotation)}
    selected.sort(key=lambda r: (refs[r["transcript_id"]], int(r["tx_start"])))
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": r["transcript_id"], "LN": int(r["length"])} for r in annotation],
    }
    records = 0
    with pysam.AlignmentFile(output_path, "wb", header=header) as target:
        for row in selected:
            read = pysam.AlignedSegment(target.header)
            read.reference_id = refs[row["transcript_id"]]
            read.reference_start = int(row["tx_start"]) - 1
            read.flag = 0
            read.mapping_quality = 60
            read.cigarstring = f'{row["fragment_length"]}M'
            read.query_sequence = row["sequence"]
            read.set_tag("NH", 1)
            for _ in range(int(row["score"])):
                records += 1
                read.query_name = f"{scenario}_{records}"
                target.write(read)
    pysam.index(str(output_path))
    return output_path, selected, records


def validate_bam(path, rows, annotation, expected_records):
    expected = Counter()
    retained = Counter()
    by_geometry = Counter()
    cds_start = {r["transcript_id"]: int(r["cds_start"]) for r in annotation}
    cds_end = {r["transcript_id"]: int(r["cds_end"]) for r in annotation}
    for row in rows:
        count = int(row["score"])
        key = (row["transcript_id"], int(row["tx_start"]) - 1,
               int(row["fragment_length"]), row["sequence"])
        expected[key] += count
        by_geometry[(int(row["fragment_length"]), int(row["site_offset"]))] += count
        site0 = int(row["site_tx"]) - 1
        assert site0 == key[1] + int(row["site_offset"])
        assert (site0 - cds_start[row["transcript_id"]]) % 3 == 0
        if cds_start[row["transcript_id"]] + 120 <= site0 < cds_end[row["transcript_id"]] - 60:
            retained[(int(row["fragment_length"]), int(row["site_offset"]))] += count
    actual = Counter()
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam:
            key = (read.reference_name, read.reference_start, read.query_length,
                   read.query_sequence)
            actual[key] += 1
    assert actual == expected
    assert sum(actual.values()) == expected_records
    return {
        "records": expected_records,
        "truth_match": True,
        "frame0_fraction": 1.0,
        "records_by_length_offset": {
            f"{length}:{offset}": count
            for (length, offset), count in sorted(by_geometry.items())
        },
        "retained_by_length_offset": {
            f"{length}:{offset}": count
            for (length, offset), count in sorted(retained.items())
        },
    }


def run_rust(scenario, stratum, bam, lengths):
    suffix = lengths.replace(":", "_")
    for mode in ("codon", "nucleotide"):
        folder = ROOT / "simulations" / scenario / "rust" / stratum / mode
        folder.mkdir(parents=True, exist_ok=True)
        numerical = folder / f"RUST_{mode}_file_{bam.name}_15_{suffix}"
        command = [
            sys.executable, "-c", "from RUST.__main__ import main; main()", mode,
            "-t", str(ROOT / "rust_transcripts.fa"), "-a", str(bam),
            "-o", "15", "-l", lengths, "-P", str(folder),
        ]
        (folder / "command.json").write_text(json.dumps(command, indent=2))
        env = os.environ.copy()
        env["PYTHONPATH"] = str(UPSTREAM)
        env["MPLCONFIGDIR"] = str(ROOT / "matplotlib_cache")
        if not numerical.exists():
            with (folder / "run.log").open("w") as log:
                subprocess.run(
                    command, env=env, stdout=log,
                    stderr=subprocess.STDOUT, check=True
                )
        if not numerical.exists() or numerical.stat().st_size == 0:
            raise RuntimeError(f"Missing RUST numerical output: {numerical}")


def main():
    annotation = read_table(ROOT / "rust_transcript_annotation.tsv")
    pysam.faidx(str(ROOT / "rust_transcripts.fa"))
    validations = {}
    for scenario in SCENARIOS:
        print(f"Preparing {scenario}", flush=True)
        folder = ROOT / "simulations" / scenario
        rows = read_table(folder / "rust_fragment_truth.tsv")
        full_bam, full_rows, full_n = write_transcript_bam(
            scenario, rows, annotation, offset15_only=False
        )
        offset_bam, offset_rows, offset_n = write_transcript_bam(
            scenario, rows, annotation, offset15_only=True
        )
        full_qc = validate_bam(full_bam, full_rows, annotation, full_n)
        offset_qc = validate_bam(offset_bam, offset_rows, annotation, offset_n)
        lengths = sorted({int(r["fragment_length"]) for r in offset_rows})
        assert 28 in lengths
        validations[scenario] = {"full": full_qc, "offset15": offset_qc,
                                 "offset15_lengths": lengths}
        (folder / "rust_input_validation.json").write_text(
            json.dumps(validations[scenario], indent=2)
        )
        print(f"RUST 28 nt: {scenario}", flush=True)
        run_rust(scenario, "length28", full_bam, "28")
        pooled_range = f"{min(lengths)}:{max(lengths)}"
        print(f"RUST offset-15 pool ({pooled_range}): {scenario}", flush=True)
        run_rust(scenario, "offset15_pooled", offset_bam, pooled_range)
    commit = subprocess.check_output(
        ["git", "-C", str(UPSTREAM), "rev-parse", "HEAD"], text=True
    ).strip()
    (ROOT / "rust_run_manifest.json").write_text(json.dumps({
        "upstream": "https://github.com/JackCurragh/RUST",
        "commit": commit,
        "version": "1.3.0 upstream Python-3 port",
        "source_modified": False,
        "scenarios": SCENARIOS,
        "analyses": ["28 nt / offset 15", "all actual-offset-15 reads"],
        "validations": validations,
    }, indent=2))
    print("All upstream RUST runs complete", flush=True)


if __name__ == "__main__":
    main()
