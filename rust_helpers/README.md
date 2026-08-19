# RUST helper scripts

This directory contains analysis helpers and workflows that are useful for
feeding coverageSim outputs into RUST-style analyses or for making technical
bias test datasets.

These files are intentionally kept outside `R/` so they are not part of the
core coverageSim package API and do not change the simulator output by being
loaded with the package.

## Contents

- `R/Learn_from_input.R`: helper functions for learning read lengths, CDS
  counts, and codon bias from real BAM input.
- `R/G_end_dropout_helpers.R`: helper functions for G-ending read dropout
  datasets.
- `workflows/Workflow coverageSim learn from real BAM.R`: real-BAM learning
  workflow.
- `workflows/Workflow coverageSim human genome only.R`: human-genome-only
  coverageSim workflow.
- `workflows/create_G_end_dropout_datasets.R`: technical G-ending dropout
  dataset workflow.

## Alignment-coordinate semantics

The RUST learning helpers treat imported alignments as complete physical
fragments. `point_reads()` anchors each alignment at its strand-aware biological
5-prime end and then applies the optional `focal_offset` to obtain an A-/P-site
or another focal position. RUST therefore expects true fragment ends in OFST or
BAM input, not coverage points that were already shifted to a ribosome site.

The physical-fragment mode in coverageSim follows that expectation. The
explicit `legacy_point` mode does not: its alignment start is the simulated
signal position and should not be interpreted as a physical 5-prime end.

The read-end helpers prefer BAM `SEQ`. Their historical FASTA fallback infers a
single contiguous reference interval from POS and reference-consuming CIGAR
width; that fallback is not suitable for spliced CIGARs containing `N`.
CoverageSim does not modify these downstream RUST helpers as part of physical
fragment generation.
