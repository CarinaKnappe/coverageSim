"""Summarise upstream RUST output and render its metafootprint plotting function."""
import csv
import importlib
import json
import os
from pathlib import Path
import sys
import numpy as np

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / 'rust_helpers/runs/2026-09-10_original-rust_read-ends'
os.environ['MPLCONFIGDIR'] = str(ROOT / 'matplotlib_cache')
sys.path.insert(0, str(REPO / 'rust_helpers/upstream/RUST'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rust_plot_compat import compatible_metafootprint

SCENARIOS = ['control', 'five_prime', 'three_prime', 'both_ends']
TITLES = {'control': 'Neutral control', 'five_prime': "5-prime bias: G preferred",
          'three_prime': "3-prime bias: C preferred", 'both_ends': 'Both ends: G ... C preferred'}


def read_profile(scenario, mode):
    path = ROOT / scenario / f'rust_{mode}' / f'RUST_{mode}_file_footprints.bam_15_28'
    lines = path.read_text().splitlines()
    split = lines.index('')
    rows = list(csv.reader(lines[:split], skipinitialspace=True))
    positions = np.array([int(v) for v in rows[0][2:]])
    features = [r[0] for r in rows[1:]]
    expected = np.array([float(r[1]) for r in rows[1:]])
    observed = np.array([[float(v) for v in r[2:]] for r in rows[1:]])
    divergence = np.array([float(v) if v.strip() != 'NA' else np.nan
                           for v in next(csv.reader([lines[split+1]], skipinitialspace=True))[2:]])
    return {'path': path, 'positions': positions, 'features': features,
            'ratios': observed / expected[:, None], 'divergence': divergence}


def draw_profile(ax, scenario, mode, profile, limits, title=True):
    module = importlib.import_module(f'RUST.{mode}')
    before = set(ax.figure.axes)
    with profile['path'].open() as handle:
        compatible_metafootprint(module)(handle, ax)
    new_axes = set(ax.figure.axes) - before
    if len(new_axes) != 1:
        raise RuntimeError('Original RUST did not draw a divergence curve')
    right = new_axes.pop()
    # Check the actual drawn curves against unchanged upstream numerical output.
    for line, expected in zip(ax.lines[:len(profile['features'])], np.log2(profile['ratios'])):
        np.testing.assert_allclose(line.get_ydata(), expected, atol=1e-12)
    np.testing.assert_allclose(right.lines[0].get_ydata(), profile['divergence'], atol=1e-12)
    right.set_ylim(0, limits['divergence'])
    right.set_yticks(np.linspace(0, limits['divergence'], 5))
    right.set_yticklabels([f'{v:.2f}' for v in np.linspace(0, limits['divergence'], 5)])
    right.set_ylabel('RUST divergence', color='blue')
    ax.set_ylim(limits['log2'])
    offset = 40 if mode == 'codon' else 120
    ends = [-5, 4] if mode == 'codon' else [-15, 12]
    for value, label in zip(ends, ["5-prime", "3-prime"]):
        ax.axvline(offset + value, color='black', linestyle=':', alpha=.65)
        ax.text(offset + value, .96, label, transform=ax.get_xaxis_transform(),
                rotation=90, va='top', ha='right', fontsize=8)
    ax.axhline(0, color='grey', linestyle=':', linewidth=.6)
    ax.set_xlabel('Distance from A-site [' + ('codons' if mode == 'codon' else 'nt') + ']')
    # The upstream nucleotide plot mistakenly says codon; only the label is corrected.
    ax.set_ylabel(('Codon' if mode == 'codon' else 'Nucleotide') + ' RUST ratio, log2(observed/expected)')
    ax.tick_params(labelsize=8)
    right.tick_params(labelsize=8)
    if title:
        ax.set_title(TITLES[scenario], fontsize=12)


def main():
    profiles = {(s,m): read_profile(s,m) for s in SCENARIOS for m in ['codon','nucleotide']}
    limits = {}
    for mode in ['codon','nucleotide']:
        values = np.concatenate([np.log2(profiles[s,mode]['ratios']).ravel() for s in SCENARIOS])
        if not np.isfinite(values).all():
            raise RuntimeError('Non-finite ratio: inspect before rendering')
        div = np.concatenate([profiles[s,mode]['divergence'] for s in SCENARIOS])
        if not np.isfinite(div).all():
            raise RuntimeError('Original RUST returned NA divergence: inspect before rendering')
        limits[mode] = {'log2': (float(np.floor(values.min()-.1)), float(np.ceil(values.max()+.1))),
                        'divergence': float(np.ceil(div.max()*11)/10)}
    summary = []
    for s in SCENARIOS:
        nt = profiles[s,'nucleotide']; codon = profiles[s,'codon']
        def at(p, position): return int(np.where(p['positions'] == position)[0][0])
        qc = json.loads((ROOT / s / 'bam_validation.json').read_text())
        summary.append({'scenario': s, 'retained_reads': qc['retained_reads'],
          'five_G_ratio': nt['ratios'][nt['features'].index('G'),at(nt,-15)],
          'three_C_ratio': nt['ratios'][nt['features'].index('C'),at(nt,12)],
          'five_divergence': codon['divergence'][at(codon,-5)],
          'three_divergence': codon['divergence'][at(codon,4)],
          'a_site_divergence': codon['divergence'][at(codon,0)],
          'max_nt_divergence_position': nt['positions'][np.argmax(nt['divergence'])]})
    with (ROOT / 'rust_signal_summary.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(summary)
    with PdfPages(ROOT / 'original_RUST_metafootprints.pdf') as pdf:
        for mode in ['codon','nucleotide']:
            for s in SCENARIOS:
                fig, ax = plt.subplots(figsize=(9,6))
                draw_profile(ax,s,mode,profiles[s,mode],limits[mode])
                fig.tight_layout()
                pdf.savefig(fig)
                fig.savefig(ROOT / s / f'metafootprint_{mode}.png', dpi=180)
                plt.close(fig)
    for mode in ['codon','nucleotide']:
        fig, axes = plt.subplots(2,2,figsize=(15,10))
        for ax,s in zip(axes.flat,SCENARIOS):
            draw_profile(ax,s,mode,profiles[s,mode],limits[mode])
        fig.suptitle('Original RUST: read-end sequence bias (' + mode + ')', fontsize=16)
        fig.tight_layout(rect=(0,0,1,.96))
        fig.savefig(ROOT / f'original_RUST_{mode}_comparison.png', dpi=160)
        plt.close(fig)
    (ROOT / 'plot_validation.json').write_text(json.dumps({
        'all_plotted_feature_curves_match_upstream': True,
        'all_plotted_divergence_curves_match_upstream': True,
        'profiles_validated': 8, 'source_numerical_files_modified': False}, indent=2))
    print(json.dumps(summary, indent=2, default=float))


if __name__ == '__main__':
    main()
