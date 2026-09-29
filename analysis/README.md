# Analysis code

Everything here *uses* coverageSim; none of it is part of the package. The
simulator lives in `../R/` and is what gets installed and published. These files
are deliberately kept outside it, so that loading the package cannot pull in an
analysis adapter and so that a built package stays limited to the simulator.

They are sourced by path, never loaded as a library.

## Layout

```
analysis/
├── R/              shared helpers, sourced by the scripts below
├── workflows/      scripts that produce a simulation run or a test dataset
└── benchmarks/
    └── rust/       comparison against the original RUST tool
        ├── *.R, *.py   our scripts
        ├── upstream/   a clone of JackCurragh/RUST (GPL-3, not redistributed
        │               here; each run manifest records the exact commit)
        ├── runs/       stored outputs, not in version control
        └── .venv-rust-original/   Python environment, rebuild from
                                   upstream/RUST/requirements.txt
```

The choros comparison lives in `benchmarks/choros/`.

## Contents

- `R/Learn_from_input.R`: learns read lengths, CDS counts and codon bias from a
  real BAM. Despite where it used to live, this has nothing to do with RUST --
  RUST is only one of the things that consumes what it learns.
- `R/Read_end_bias_helpers.R`, `R/G_end_dropout_helpers.R`,
  `R/end_bias_library_helpers.R`: measuring and simulating read-end bias.
- `workflows/Workflow coverageSim learn from real BAM.R`: learns its settings
  from a real library, then simulates.
- `workflows/Workflow coverageSim human genome only.R`: uses the real human
  genome and annotation but invents the expression; reads nothing real.
- `workflows/create_G_end_dropout_datasets.R`: builds datasets in which reads
  ending in G are dropped on purpose, to test bias correction.
- `benchmarks/rust/`: builds benchmark libraries with a known injected bias and
  checks whether the original RUST recovers it.

Tests for this code are in `../tests/testthat/`, so they run with the rest of
the suite.

## Alignment-coordinate semantics

The RUST learning helpers treat imported alignments as complete simulated RPF
fragments. `point_reads()` anchors each alignment at its strand-aware biological
5-prime end and then applies the optional `focal_offset` to obtain an A-/P-site
or another focal position. RUST therefore expects true fragment ends in OFST or
BAM input, not coverage points that were already shifted to a ribosome site.

The simulated-RPF mode in coverageSim follows that expectation. The
explicit `legacy_point` mode does not: its alignment start is the simulated
signal position and should not be interpreted as a simulated RPF 5-prime end.

The read-end helpers prefer BAM `SEQ`. Their historical FASTA fallback infers a
single contiguous reference interval from POS and reference-consuming CIGAR
width; that fallback is not suitable for spliced CIGARs containing `N`.
CoverageSim does not modify these downstream RUST helpers as part of simulated-RPF
fragment generation.
