"""Generate mixed-geometry libraries from fresh CovSim counts/coverage weights."""
from collections import Counter, defaultdict
import csv
import importlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
import pysam
from natural_end_model import *
from rust_plot_compat import render_original_plot

REPO = Path(__file__).resolve().parents[3]
ROOT = REPO / 'analysis/benchmarks/rust/runs/2026-09-10_natural-read-ends'
UPSTREAM = REPO / 'analysis/benchmarks/rust/upstream/RUST'
os.environ['MPLCONFIGDIR'] = str(ROOT / 'matplotlib_cache')
sys.path.insert(0,str(UPSTREAM))


def read_table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def generate_libraries(rep, annotation, counts, site_data):
    folder = ROOT / rep
    if folder.exists():
        raise RuntimeError(f'Refusing to overwrite {folder}')
    folder.mkdir()
    fasta=pysam.FastaFile(str(ROOT/'transcripts.fa'))
    header={'HD':{'VN':'1.6','SO':'coordinate'},
            'SQ':[{'SN':a['transcript_id'],'LN':int(a['length'])} for a in annotation]}
    outputs={}; handles={}; writers={}; total=Counter(); qc={}
    for s in SCENARIOS:
        (folder/s).mkdir()
        outputs[s]=pysam.AlignmentFile(folder/s/'footprints.bam','wb',header=header)
        handles[s]=(folder/s/'fragment_truth.tsv').open('w')
        writers[s]=csv.writer(handles[s],delimiter='\t')
        writers[s].writerow(['transcript_id','start0','site0','fragment_length','site_offset',
                            'sequence','five_kmer','three_kmer','five_weight','three_weight','score'])
        qc[s]={'lengths':Counter(),'offsets':Counter(),'frame_common15':Counter(),
               'five_kmers':Counter(),'three_kmers':Counter(),'transcripts':Counter(),
               'retained_pooled':Counter(),'retained_28':Counter()}
    try:
        for tx_index, a in enumerate(annotation):
            tx=a['transcript_id']; sequence=fasta.fetch(tx)
            rows=site_data[tx]
            sites=np.array([int(r['site0']) for r in rows])
            base=np.array([float(r['weight']) for r in rows])
            seed=202609102 + int(rep[-1])*100000 + tx_index
            latent=noisy_site_weights(base,np.random.default_rng(seed))
            candidates=fragment_candidates(sequence,sites,latent)
            order=np.lexsort((candidates['offset'],candidates['length'],candidates['start']))
            n=int(counts[tx][rep])
            for s_index,s in enumerate(SCENARIOS):
                selected=sample_library(candidates,n,s,np.random.default_rng(seed+1000000*(s_index+1)))
                assert selected.sum()==n
                qc[s]['transcripts'][tx]=n
                for k in order:
                    score=int(selected[k])
                    if not score: continue
                    start=int(candidates['start'][k]); site=int(candidates['site'][k])
                    width=int(candidates['length'][k]); offset=int(candidates['offset'][k])
                    seq=candidates['sequence'][k]
                    assert start+offset==site and (site-int(a['cds_start']))%3==0
                    assert seq==sequence[start:start+width] and len(seq)==width
                    writers[s].writerow([tx,start,site,width,offset,seq,candidates['five_kmer'][k],
                        candidates['three_kmer'][k],candidates['five_weight'][k],candidates['three_weight'][k],score])
                    read=pysam.AlignedSegment(outputs[s].header)
                    read.reference_id=tx_index; read.reference_start=start; read.flag=0
                    read.mapping_quality=60; read.cigarstring=f'{width}M'; read.query_sequence=seq
                    read.set_tag('NH',1)
                    for _ in range(score):
                        total[s]+=1; read.query_name=f'{rep}_{s}_{total[s]}'
                        outputs[s].write(read)
                    q=qc[s]
                    q['lengths'][width]+=score; q['offsets'][f'{width}:{offset}']+=score
                    q['frame_common15'][(start+15-int(a['cds_start']))%3]+=score
                    q['five_kmers'][candidates['five_kmer'][k]]+=score
                    q['three_kmers'][candidates['three_kmer'][k]]+=score
                    if int(a['cds_start'])+120 <= start+15 < int(a['cds_end'])-60:
                        q['retained_pooled'][tx]+=score
                        if width==28:q['retained_28'][tx]+=score
            if tx_index%50==49: print(f'{rep}: {tx_index+1}/300 transcripts exported',flush=True)
    finally:
        fasta.close()
        for s in SCENARIOS: outputs[s].close(); handles[s].close()
    for s in SCENARIOS:
        assert total[s]==1500000
        bam=folder/s/'footprints.bam'; pysam.index(str(bam))
        validate_export(bam,folder/s/'fragment_truth.tsv')
        record={'records':total[s], 'exact_truth_match':True, **qc[s]}
        (folder/s/'validation.json').write_text(json.dumps(record,indent=2))


def validate_export(bam,truth):
    expected=Counter()
    for row in read_table(truth):
        expected[(row['transcript_id'],int(row['start0']),row['sequence'])]+=int(row['score'])
    actual=Counter()
    with pysam.AlignmentFile(bam,'rb') as reads:
        for read in reads:
            assert read.flag==0 and read.query_length in LENGTHS
            assert read.cigarstring==f'{read.query_length}M'
            actual[(read.reference_name,read.reference_start,read.query_sequence)]+=1
    assert actual==expected,'Exported BAM differs from fragment truth'


def run_analysis(rep,scenario,stratum,mode):
    folder=ROOT/rep/scenario/f'rust_{stratum}_{mode}'
    length_arg='27:31' if stratum=='pooled' else '28'
    suffix='27_31' if stratum=='pooled' else '28'
    raw=folder/f'RUST_{mode}_file_footprints.bam_15_{suffix}'
    command=[sys.executable,'-c','from RUST.__main__ import main; main()',mode,
             '-t',str(ROOT/'transcripts.fa'),'-a',str(ROOT/rep/scenario/'footprints.bam'),
             '-o','15','-l',length_arg,'-P',str(folder)]
    (ROOT/rep/scenario/f'{stratum}_{mode}_command.json').write_text(json.dumps(command,indent=2))
    if not raw.exists():
        env=os.environ.copy();env['PYTHONPATH']=str(UPSTREAM)
        with (ROOT/rep/scenario/f'{stratum}_{mode}.log').open('w') as log:
            subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    assert raw.exists()
    original_plot=folder/f'RUST_{mode}_metafootprint_footprints.bam_15_{suffix}.png'
    if not original_plot.exists():
        render_original_plot(importlib.import_module(f'RUST.{mode}'),raw,
            folder/f'RUST_{mode}_metafootprint_footprints.bam_15_{suffix}_compat.png')


def main():
    annotation=read_table(ROOT/'transcript_annotation.tsv')
    counts={r['transcript_id']:r for r in read_table(ROOT/'transcript_counts.tsv')}
    site_data=defaultdict(list)
    for r in read_table(ROOT/'covsim_site_weights.tsv'): site_data[r['transcript_id']].append(r)
    if not (ROOT/'transcripts.fa.fai').exists():pysam.faidx(str(ROOT/'transcripts.fa'))
    for rep in ['rep1','rep2']:
        if not all((ROOT/rep/s/'validation.json').exists() for s in SCENARIOS):
            generate_libraries(rep,annotation,counts,site_data)
        for s in SCENARIOS:
            for stratum in ['pooled','length28']:
                for mode in ['codon','nucleotide']:
                    print(f'Upstream RUST: {rep}/{s}/{stratum}/{mode}',flush=True)
                    run_analysis(rep,s,stratum,mode)
    print('Naturalistic synthetic test completed',flush=True)


if __name__=='__main__':main()
