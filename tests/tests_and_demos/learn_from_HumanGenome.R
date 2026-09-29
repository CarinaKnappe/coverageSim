### CoverageSim with given input
#
#
# input is: BAMs + FASTA + GTF/TxDb + sample metadata
#
# The genome gives coordinates and sequences and the Ribo-seq data is what coverageSim learns coverage behavior from.
#
# Genome data:
# - FASTA: the reference genome sequence, e.g. human chromosomes.
# - GTF: genome annotation, e.g. genes, transcripts, exons, CDS coordinates, UTRs.
# - TxDb: an R/Bioconductor database version of the GTF annotation. ORFik uses this to load regions like cds, mrna, leader, trailer.
# So genome data answers: “Where are the genes/CDSs/uORFs in the genome, and what sequence do they have?”
#
# Ribo-seq data:
# - BAMs: aligned ribosome profiling reads. These are the actual Ribo-seq footprints mapped to the genome.
# - sample metadata: table describing which BAM belongs to which sample, condition, replicate, library type, etc.
# So Ribo-seq data answers: “Where did ribosome footprints occur in this experiment, and how many reads are there per sample?”

# For coverageSim learning, the workflow is roughly:
# use real Ribo-seq FASTQ/BAM + matching human FASTA + matching GTF/TxDb
# -> ORFik experiment object
# -> estimate coverage/codon/read-length/frame/region features
# -> feed learned parameters into simNGScoverage()
#
# In coverageSim terms, the learned pieces can become:
# - seq_bias: learned codon or amino-acid bias table
# - read_lengths_per: learned footprint length distribution
# - region_proportion: learned leader/CDS/trailer/uORF proportions
# - auto_correlation: learned local coverage structure
# - count table parameters: learned gene/sample count distributions

# RWA FILES FOR PROCESSING BEFORE INPUT HERE IS:
# RiboSeq sequences from Human HeLa control: https://www.ncbi.nlm.nih.gov/sra/SRX27804664[accn] (fastq)
# Human Genome: https://www.gencodegenes.org/human/ GRCh.38.p14
  # Basic gene annotation - PRI - GTF
  # Genome Sequence primary assembly GRCh38 - PRI - Fasta


