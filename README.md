# coverageSim
### Simulating Ribosome profiling data for tool validation
![](inst/images/coverageSim_overview.png)


## How to install
```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("Roleren/coverageSim")
```

## Simple example

Lets make a genome with 6 genes, each being coding (they have a CDS) and
having 1 translated uORF:

```r
library(coverageSim)
library(ORFik)
## Simple example
# 6 genes on 6 chromosomes (1 active uORF each)
simGenome6 <- simGenome(n = 6, max_uorfs = 1)
# Simulate Ribo-seq only
gene_count_table <-simCountTables(loadRegion(simGenome6["txdb"], "cds"),
libtypes = "RFP", print_statistics = FALSE)
region_count_table <- simCountTablesRegions(gene_count_table,
     regionsToSample = c("leader", "cds", "trailer"))
df <- simNGScoverage(simGenome6, region_count_table)
```

The ORFik experiment object (df), now contains all linkers to resulting files.

## Convert to other track formats

```r
# Convert to bigwig
ORFik::convert_to_bigWig(df)

```
ORFik supports multiple other formats to convert to, also coverageSim internally has
an option to export simulated libraries directly as `.ofst`, `sam`, or `bam`
through `simNGScoverage(libFormats = ...)`.

## Relevant paths for other programs

```r
# Transcript annotation (.gtf)
gtf_path <- ORFik:::getGtfPathFromTxdb(loadTxdb(df))
# Genome sequences (.fasta)
fasta_genome_path <- df@fafile
# NGS track files (.bigwig, pairwise for forward and reverse strands)
bigwig_paths <- filepath(df, "bigwig")
```

## Export to IGV

- Open IGV,
- Press "Genomes" tab (top left), "Create .genome file"
- Unique Identifier: "Sim genome", Descriptive name: "Sim genome", 
- FASTA file: input fasta_genome_path above
- Gene file: input gtf_path above
- Press OK

Now load bigwig files by: 

- Press "File" (top left), "Load from File"
- load all paths from bigwig_paths above

You now have genome, gtf and tracks loaded in IGV

## Learn end biases independently of the simulation inputs

`learn_end_bias()` estimates 5-prime and 3-prime sequence preferences from a
complete-fragment genomic Ribo-seq BAM. Supply its matching reference and
annotation, and an explicit length-specific offset table. These learning inputs
can differ from the genome, count table and region proportions used for simulation.

```r
devtools::load_all(".")

# Select matching, unambiguous transcript/CDS annotations for the real sample.
transcripts <- ORFik::loadRegion(real_txdb, "mrna")
cds <- ORFik::loadRegion(real_txdb, "cds", names.keep = names(transcripts))

# Learn supported lengths, their empirical probabilities and one A-site offset
# per length from frame evidence in the same BAM.
geometry <- learn_fragment_geometry(
  bam = "real_sample.bam", transcripts = transcripts, cds = cds
)
learned <- learn_end_bias(
  bam = "real_sample.bam", fasta = "matching_reference.fa",
  transcripts = transcripts, cds = cds,
  fragment_geometry = geometry, k = 1, by_length = TRUE
)
saveRDS(learned, "learned_end_bias.rds")

# No original BAM is needed for later simulations. Choose the simulation's
# own geometry, genome and region counts independently.
simulation_geometry <- list(source = "default", site_reference = "a_site")
simulation_geometry$five_prime_bias <- learned$five_prime_bias
simulation_geometry$three_prime_bias <- learned$three_prime_bias
simulation_geometry$codon_bias <- learned$codon_bias
simulation_geometry$frame_bias <- learned$frame_bias
simulation_geometry$five_prime_bias$strength <- 0.5
simulation_geometry$three_prime_bias$strength <- 1
simNGScoverage(
  simGenome = simulated_genome,
  count_table = simulated_region_counts,
  seq_bias = learned$sequence_bias,
  fragment_geometry = simulation_geometry
)
```

The fit compares observed reads with possible CDS fragments, including positions
with zero reads. It jointly estimates both end effects and nuisance codon effects,
conditioning on transcript-by-length totals. A positive ridge penalty stabilizes
sparse motifs. `k = 2` or `k = 3` allows longer motifs; `by_length = TRUE` estimates
separate end, codon and frame profiles per supported length. Missing length/motif combinations in a
profile retain the simulator's neutral fallback. The fitted codon effects are
returned as `learned$sequence_bias`. This table also stores the robustly estimated
Dirichlet-multinomial concentration. When it is passed as `seq_bias`, the default
`dmn_alpha_scale = NULL` uses that learned concentration automatically. An
explicit numeric `dmn_alpha_scale` always overrides it; profiles without a
learned value fall back to `1`.

Inspect `learned$diagnostics` for used/excluded read counts, motif support and fit
convergence. The BAM is read in chunks, retaining counts per distinct alignment.
The end learner requires one supplied or learned offset per length and complete CDS
annotations. It excludes clips, indels, paired reads, low-quality/secondary/
supplementary mappings, NH>1, ambiguous transcript assignments, and reads whose
assigned site is not at a CDS codon boundary. Duplicate flags alone are retained.
It estimates transferable sequence preferences under these assumptions, rather
than identifying an enzyme-specific effect or providing confidence intervals.
Length probabilities and offsets come from `learn_fragment_geometry()`; region
proportions come from `learn_region_proportions()`. Coverage roughness is estimated from transcript-level overdispersion
after accounting for fitted codon and end preferences. Strong position-specific
biological effects or unmodelled mapping biases can still affect the estimate.
The fit also returns a non-negative codon autocorrelation kernel and auditable
zero, spike, peak and gap summaries. Pass the kernel as the RFP value in
`auto_correlation`; use the QC tables to compare real and simulated libraries.

Region allocation can be learned independently from the same or another BAM:

```r
regions <- list(leader = leader_ranges, cds = cds_ranges,
                trailer = trailer_ranges, uorf = uorf_ranges)
region_fit <- learn_region_proportions(
  bam = "sample.bam",
  transcripts = learning_transcripts,
  regions = regions,
  fragment_geometry = learned_geometry
)

region_counts <- simCountTablesRegions(
  simulated_counts,
  region_proportion = region_fit$region_proportion
)
```

Reads compatible with several transcripts are excluded. A-sites falling in
overlapping regions are assigned by the documented `uorf`, `cds`, `leader`,
`trailer` priority by default; fractional allocation and exclusion are also
available. The fit reports both global simulator-ready proportions and
per-transcript estimates, with a configurable pseudocount for sparse regions.
