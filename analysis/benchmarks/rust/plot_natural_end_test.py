"""Plot unchanged upstream RUST results for the heterogeneous end-bias benchmark."""
import csv, importlib, json, os, sys
from pathlib import Path
import numpy as np
REPO=Path(__file__).resolve().parents[3]
ROOT=REPO/'analysis/benchmarks/rust/runs/2026-09-10_natural-read-ends'
os.environ['MPLCONFIGDIR']=str(ROOT/'matplotlib_cache')
sys.path.insert(0,str(REPO/'analysis/benchmarks/rust/upstream/RUST'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from rust_plot_compat import compatible_metafootprint
from natural_end_model import SCENARIOS
TITLES={'control':'Control','five_prime':'5-prime bias','three_prime':'3-prime bias','both_ends':'Both ends'}


def profile(rep,scenario,stratum,mode):
    suffix='27_31' if stratum=='pooled' else '28'
    path=ROOT/rep/scenario/f'rust_{stratum}_{mode}'/f'RUST_{mode}_file_footprints.bam_15_{suffix}'
    lines=path.read_text().splitlines(); split=lines.index('')
    rows=list(csv.reader(lines[:split],skipinitialspace=True))
    positions=np.array([int(x) for x in rows[0][2:]])
    expected=np.array([float(r[1]) for r in rows[1:]])
    observed=np.array([[float(x) for x in r[2:]] for r in rows[1:]])
    divergence=np.array([float(x) if x.strip()!='NA' else np.nan for x in next(csv.reader([lines[split+1]],skipinitialspace=True))[2:]])
    return dict(path=path,positions=positions,ratios=observed/expected[:,None],divergence=divergence)


def draw(ax,p,mode,limits,title):
    before=set(ax.figure.axes)
    with p['path'].open() as f:compatible_metafootprint(importlib.import_module('RUST.'+mode))(f,ax)
    axes=set(ax.figure.axes)-before
    if len(axes)!=1:raise RuntimeError('Missing upstream divergence curve')
    right=axes.pop()
    for line,expected in zip(ax.lines[:len(p['ratios'])],np.log2(p['ratios'])):
        np.testing.assert_allclose(line.get_ydata(),expected,atol=1e-12)
    np.testing.assert_allclose(right.lines[0].get_ydata(),p['divergence'],atol=1e-12)
    ax.set_ylim(limits[0]);right.set_ylim(0,limits[1])
    ticks=np.linspace(0,limits[1],5);right.set_yticks(ticks)
    right.set_yticklabels([f'{x:.2f}' for x in ticks])
    right.set_ylabel('RUST divergence',color='blue')
    ax.set_title(title);ax.set_xlabel('Distance from assigned A-site ['+('codons' if mode=='codon' else 'nt')+']')
    ax.tick_params(labelsize=8);right.tick_params(labelsize=8)
    ax.axhline(0,color='grey',ls=':',lw=.5)


def main():
    data={(r,s,t,m):profile(r,s,t,m) for r in ['rep1','rep2'] for s in SCENARIOS
          for t in ['pooled','length28'] for m in ['codon','nucleotide']}
    limits={}
    for mode in ['codon','nucleotide']:
        ps=[v for k,v in data.items() if k[-1]==mode]
        ratios=np.concatenate([np.log2(p['ratios']).ravel() for p in ps])
        div=np.concatenate([p['divergence'] for p in ps])
        if not np.isfinite(ratios).all() or not np.isfinite(div).all():
            raise RuntimeError('Non-finite upstream values need inspection')
        limits[mode]=((float(np.floor(ratios.min()*2-.1)/2),float(np.ceil(ratios.max()*2+.1)/2)),
                      float(np.ceil(div.max()*110)/100))
    summary=[]
    for rep in ['rep1','rep2']:
        for stratum in ['pooled','length28']:
            for scenario in SCENARIOS:
                p=data[rep,scenario,stratum,'nucleotide'];pos=p['positions']
                q=json.loads((ROOT/rep/scenario/'validation.json').read_text())
                background=p['divergence'][(pos < -25)|(pos > 25)]
                control=data[rep,'control',stratum,'nucleotide']['divergence']
                delta=p['divergence']-control
                summary.append(dict(replicate=rep,scenario=scenario,stratum=stratum,
                    retained_reads=sum(q['retained_pooled' if stratum=='pooled' else 'retained_28'].values()),
                    retained_transcripts=len(q['retained_pooled' if stratum=='pooled' else 'retained_28']),
                    five_peak=float(p['divergence'][(pos>=-16)&(pos<=-11)].max()),
                    three_peak=float(p['divergence'][(pos>=9)&(pos<=18)].max()),
                    five_excess=float(delta[(pos>=-16)&(pos<=-11)].sum()),
                    three_excess=float(delta[(pos>=9)&(pos<=18)].sum()),
                    background_mean=float(background.mean()),a_site=float(p['divergence'][pos==0][0])))
    with (ROOT/'signal_summary.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t');w.writeheader();w.writerows(summary)
    with PdfPages(ROOT/'natural_RUST_metafootprints.pdf') as pdf:
        for stratum in ['pooled','length28']:
            for mode in ['codon','nucleotide']:
                for rep in ['rep1','rep2']:
                    fig,axes=plt.subplots(2,2,figsize=(15,10))
                    for ax,s in zip(axes.flat,SCENARIOS):draw(ax,data[rep,s,stratum,mode],mode,limits[mode],TITLES[s])
                    label='27-31 nt, common offset 15' if stratum=='pooled' else '28 nt only, offset 15'
                    fig.suptitle(f'Upstream RUST - {rep} - {label} - {mode}',fontsize=16)
                    fig.tight_layout(rect=(0,0,1,.96));pdf.savefig(fig)
                    fig.savefig(ROOT/f'{rep}_{stratum}_{mode}.png',dpi=140);plt.close(fig)
    (ROOT/'plot_validation.json').write_text(json.dumps(dict(profiles=32,
        plotted_curves_match_upstream=True,numerical_files_changed=False),indent=2))
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
