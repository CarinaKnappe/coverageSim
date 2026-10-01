# choros benchmark: real human Ribo-seq against simulated coverage

Two workflows that answer the same question from opposite directions. The first
runs a real human Ribo-seq library, where sequence-dependent technical bias is
expected and serves as a positive control. The second runs a fully simulated
sample generated without any deliberate 5'- or 3'-end sequence bias, where the
correct answer is that there is none to find.

The upstream method is attributed in `../README.md`.

## The two workflows

| | Real human | Simulated |
|---|---|---|
| Prepare | `prepare_choros_real_human.R` | `prepare_choros_covsim_human_genome_only.R` |
| Fit | `choros_real_human.R` | `choros_covsim_human_genome_only.R` |
| Input | `SRR32491292.sorted.unique_nh1.ofst` | `RFP_WT_1.bam` from a coverageSim run |

All four scripts source `choros_utils.R` from their own directory, so they work
wherever this folder is placed. Preparation has to finish before the fitting
script is started.

## Running them

Nothing is hardcoded to a particular machine; the locations come from four
environment variables.

```r
Sys.setenv(
  CHOROS_REAL_BASE  = "<dataset with reads/ and genome/>",   # real human
  CHOROS_COVSIM_BASE = "<coverageSim run directory>",        # simulated
  CHOROS_OUTPUT_DIR = "<where results should go>",
  CHOROS_RUN_ID     = "RFP_WT_1"                             # which sample
)
source("prepare_choros_covsim_human_genome_only.R")
source("choros_covsim_human_genome_only.R")
```

A simulated dataset is expected to hold `reads/<run_id>.bam`,
`genome/human_flavoured_sim.fasta` and `genome/human_flavoured_sim.gtf.db`.
Runs whose genome files are named differently need
`resolve_covsim_config()` adjusted.

The preparation script repairs the zero chromosome lengths that a simulated BAM
header can carry, taking the real lengths from the FASTA index, and validates
every alignment against the reference. It reads the BAM and does not modify it.

## What the inputs contribute

The reference FASTA supplies the sequences for codons and for the 5'/3'
fragment-end features; its index supplies the authoritative chromosome lengths.
The transcript annotation supplies transcripts, CDSs and 5' leaders. The
workflow keeps transcripts with at least 1 nt of leader and at least 231 nt of
CDS, then the longest per gene.

## Stage 1 — preparation

1. load and align transcript, CDS and leader annotations;
2. import the mapped reads;
3. map each 5' read end into transcript coordinates;
4. select read lengths — for real data by abundance and clear 3-nt frame
   enrichment, for simulated data by keeping the configured lengths 28, 29 and
   30 nt, since filtering those on observed frame would only confirm what was
   put in;
5. infer length-specific offsets from the observed peak at the start codon;
6. construct frame-specific `d5`, `d3` and A-site codon indices;
7. write a validated input RDS and QC tables.

`prepared/` then holds `*_choros_input.rds` (counts, sequences, geometry),
`*_offsets.csv`, `*_periodicity_qc.csv` and `*_tis_profile.csv`.

## Stage 2 — model and correction

1. build the zero-filled transcript/codon/geometry universe;
2. select the 250 highest-density training transcripts;
3. fit a negative-binomial base model with codon effects, GC and geometry;
4. fit a full model that also includes the 5'/3' end sequences and their
   interactions with `d5` and `d3`;
5. keep the full model when `BIC_full < BIC_base`, otherwise the base model — on
   data without end bias the base model wins, which is the correct negative
   result;
6. compute correction factors and corrected counts;
7. write summary tables and diagnostic plots.

## Results

- `*_choros_summary.csv` — read count, both BICs, selected model, absolute bias
- `*_choros_corrected_counts.csv.gz` — raw and corrected footprint counts
- `*_choros_diagnostics.pdf` — length and frame, geometry, start and stop
  profiles, positional bias, metagene

Each run writes below its own `CHOROS_OUTPUT_DIR`, so one analysis cannot
overwrite another.
