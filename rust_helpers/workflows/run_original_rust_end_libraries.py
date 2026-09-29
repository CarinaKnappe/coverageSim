"""Export CovSim library candidates in transcript coordinates; run upstream RUST."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from collections import Counter
import pysam

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / 'rust_helpers/runs/2026-09-10_original-rust_read-ends'
UPSTREAM = REPO / 'rust_helpers/upstream/RUST'
SCENARIOS = ['control', 'five_prime', 'three_prime', 'both_ends']


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def export_bam(scenario, annotation):
    folder = ROOT / scenario
    bam = folder / 'footprints.bam'
    if bam.exists():
        return bam, read_table(folder / 'fragment_truth.tsv')
    header = {'HD': {'VN': '1.6', 'SO': 'coordinate'},
              'SQ': [{'SN': row['transcript_id'], 'LN': int(row['length'])} for row in annotation]}
    refs = {row['transcript_id']: i for i, row in enumerate(annotation)}
    rows = read_table(folder / 'fragment_truth.tsv')
    rows.sort(key=lambda r: (refs[r['transcript_id']], int(r['tx_start'])))
    counter = 0
    with pysam.AlignmentFile(bam, 'wb', header=header) as output:
        for row in rows:
            count = int(row['score'])
            if not count:
                continue
            read = pysam.AlignedSegment(output.header)
            read.reference_id = refs[row['transcript_id']]
            read.reference_start = int(row['tx_start']) - 1
            read.flag = 0  # transcript sequences and reads are in biological orientation
            read.mapping_quality = 60
            read.cigarstring = f"{row['fragment_length']}M"
            read.query_sequence = row['sequence']
            read.set_tag('NH', 1)
            for _ in range(count):
                counter += 1
                read.query_name = f'{scenario}_{counter}'
                output.write(read)
    pysam.index(str(bam))
    assert counter == 1500000
    return bam, rows


def validate_bam(bam, rows, annotation):
    expected = Counter()
    for row in rows:
        if int(row['score']):
            expected[(row['transcript_id'], int(row['tx_start']) - 1, row['sequence'])] += int(row['score'])
    actual = Counter()
    retained = Counter()
    starts = {r['transcript_id']: int(r['cds_start']) for r in annotation}
    ends = {r['transcript_id']: int(r['cds_end']) for r in annotation}
    with pysam.AlignmentFile(bam, 'rb') as reads:
        for read in reads:
            assert read.flag == 0 and read.cigarstring == '28M'
            assert read.query_length == 28
            tx = read.reference_name
            actual[(tx, read.reference_start, read.query_sequence)] += 1
            site = read.reference_start + 15
            assert (site - starts[tx]) % 3 == 0
            if starts[tx] + 120 <= site < ends[tx] - 60:
                retained[tx] += 1
    assert actual == expected, 'BAM differs from simulated library counts/sequences/positions'
    assert len(retained) == 300
    return {'records': sum(actual.values()), 'retained_reads': sum(retained.values()),
            'retained_transcripts': len(retained), 'minimum_retained_per_transcript': min(retained.values()),
            'exact_truth_match': True, 'frame0_fraction': 1.0}


def run_rust(scenario, bam, mode):
    folder = ROOT / scenario / f'rust_{mode}'
    numerical = folder / f'RUST_{mode}_file_footprints.bam_15_28'
    env = os.environ.copy()
    env['PYTHONPATH'] = str(UPSTREAM)
    env['MPLCONFIGDIR'] = str(ROOT / 'matplotlib_cache')
    command = [sys.executable, '-c', 'from RUST.__main__ import main; main()', mode,
               '-t', str(ROOT / 'transcripts.fa'), '-a', str(bam), '-o', '15',
               '-l', '28', '-P', str(folder)]
    (ROOT / scenario / f'{mode}_command.json').write_text(json.dumps(command, indent=2))
    if not numerical.exists():
        with (ROOT / scenario / f'{mode}.log').open('w') as log:
            subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    prefix = 'codon' if mode == 'codon' else 'nucleotide'
    assert list(folder.glob(f'RUST_{prefix}_file_*')), 'Missing numerical output'
    # Upstream port passes a map iterator to matplotlib; recover its plot separately.
    import importlib
    from rust_plot_compat import render_original_plot
    sys.path.insert(0, str(UPSTREAM))
    os.environ['MPLCONFIGDIR'] = str(ROOT / 'matplotlib_cache')
    module = importlib.import_module(f'RUST.{mode}')
    if not (folder / f'RUST_{prefix}_metafootprint_footprints.bam_15_28.png').exists():
        render_original_plot(module, numerical, folder / f'RUST_{prefix}_metafootprint_footprints.bam_15_28_compat.png')


def main():
    annotation = read_table(ROOT / 'transcript_annotation.tsv')
    pysam.faidx(str(ROOT / 'transcripts.fa'))
    commit = subprocess.check_output(['git', '-C', str(UPSTREAM), 'rev-parse', 'HEAD'], text=True).strip()
    manifest = {'upstream': 'https://github.com/JackCurragh/RUST', 'commit': commit,
                'version': '1.3.0 (upstream Python 3 port)', 'source_modified': False, 'plot_adapter': 'map iterator converted to list in copied plotting function; numerical output unmodified'}
    (ROOT / 'rust_version.json').write_text(json.dumps(manifest, indent=2))
    for scenario in SCENARIOS:
        print(f'Exporting {scenario}', flush=True)
        bam, rows = export_bam(scenario, annotation)
        qc = validate_bam(bam, rows, annotation)
        (ROOT / scenario / 'bam_validation.json').write_text(json.dumps(qc, indent=2))
        for mode in ['codon', 'nucleotide']:
            print(f'Original RUST {mode}: {scenario}', flush=True)
            run_rust(scenario, bam, mode)
        print(f'Completed {scenario}: {qc}', flush=True)
    print('All original RUST runs complete', flush=True)


if __name__ == '__main__':
    main()
