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
