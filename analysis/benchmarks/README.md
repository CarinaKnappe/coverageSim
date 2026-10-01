# Benchmarks against external analysis tools

Each directory here answers the same question for a different tool: if a known
bias is built into simulated coverage on purpose, does an established analysis
tool recover it -- and does it correctly report nothing where nothing was put?

## What is done with each tool

For both tools, two things happen, and it matters which is which.

**The original is run.** Under `rust/`, `run_human_learning_rust.py`,
`run_natural_end_test.py` and `run_original_rust_end_libraries.py` execute the
original RUST as a subprocess, against the clone kept locally in
`rust/upstream/`. For choros, `run_choros.R` (kept outside this repository)
calls the installed package through `choros::`. Running the original is the
stronger comparison, since it removes any doubt about whether a difference
comes from the tool or from our version of it.

**The method is also reimplemented.** `choros/choros_utils.R` reimplements the
choros bias model in R; a separate R reimplementation of RUST exists outside
this repository. These exist because the originals cannot always be applied to
the data in the form we have it, and because a reimplementation can be
instrumented in ways a subprocess cannot.

## What that means for provenance

No code is copied. Comparing every line of 25 characters or more between
`choros/choros_utils.R` and the choros package, in the same language, gives
zero matches; the RUST reimplementation is in R while the original is Python,
so the question does not arise there.

The reimplementations are close, however, and were written with the original
source in view -- they are not independent work. The GC term in
`compute_gc_dt()` computes the same quantity as `compute_rpf_gc()` in choros,
including the constant `- 7L`, which is that function's
`- 1 - 3 * num_omit_codons` evaluated for `omit = "APE"`, and including its
choice of dividing by the full fragment length while omitting nine bases from
the numerator. The regression is the same model, fitted with `fixest::fenegbin`
absorbing the transcript effect where the original uses `MASS::glm.nb` with
transcript as a term.

Anyone reusing this should cite the upstream methods and keep the notices
below. The absence of copied text does not by itself make a reimplementation
independent of the work it reimplements.

## The two upstream tools

**RUST** — <https://github.com/JackCurragh/RUST>, GPL-3.

A clone is kept at `rust/upstream/` so the benchmarks can run the original
rather than a reimplementation of it, and every stored run records the commit
it was run against. That directory is excluded from version control and is not
redistributed here; to reproduce a run, clone the repository at the commit the
run's manifest names. `rust/.venv-rust-original/` is the Python environment it
needs, also excluded, rebuildable from `rust/upstream/RUST/requirements.txt`.

**choros** — <https://github.com/lareaulab/choros>, MIT.

    Copyright (c) 2023 Liana Lareau

The choros bias model is reimplemented in `choros/choros_utils.R`, so that
reimplementation is treated as derived from their work and carries the notice
above, as the MIT licence requires. No choros code is copied: the package is
installed separately and called through `choros::`.

## What is not tracked, and why

`rust/upstream/`, `rust/.venv-rust-original/` and `choros/results/` are
excluded. The first two are someone else's software and a rebuildable
environment; the third is output, which belongs with the other stored runs
under `coverageSim_data/runs/` rather than inside the repository.
