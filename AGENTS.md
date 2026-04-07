# AGENTS.md

## Project overview
Simulation project for ribo-seq coverage analysis.

## Guidelines
- Keep functions small and readable
- Do not modify data files
- Add tests for new features

## General workflow
- First, catch up on the latest local changes and understand the repository.
- Always inspect related code before making assumptions or edits.
- If necessary, follow dependencies deeper (e.g. ORFik / Bioconductor) before making changes.

## Plan mode
- Before implementing, propose 2–3 alternative approaches for this change.
- Focus on correctness vs performance trade-offs.
- Do not modify files yet.

## Before retrying any command outside the sandbox:
- explain why the sandbox execution failed
- identify which capability is missing (filesystem, network, system libs, etc.)
- propose a sandbox-compatible alternative if possible
- ask for confirmation before proceeding

## Running code
- Always use `devtools::load_all()` when working with coverageSim.
- Do NOT rely on an installed version of the package.

## Running R commands
- Run R commands using:

  env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); <COMMAND>'

- This avoids `LC_ALL='C.UTF-8'` warnings from testthat/withr.

## Testing (REQUIRED)
- Try every command in the sandbox first.
- Only escalate if a specific command is still blocked and that block is clearly required for the task.
- Tell exactly which command hit the limit.
- Ensure tests are run in the same environment using `devtools::load_all()`.
- Do not run tests against an installed version of the package.
- Every Codex task MUST include a test that proves the implementation works.
- Tests must validate the actual behavior being implemented or changed (not just superficial checks).
- Tests must be placed in the `tests/testthat/` directory.
- Do not consider a task complete unless the new or updated tests pass.


## Example commands
Run tests:
env -u LC_ALL R --vanilla -q -e 'devtools::load_all("."); testthat::test_dir("tests/testthat")'

## Expectations
- Prefer correctness over speed.
- Verify assumptions by reading code before modifying anything.
- Do not declare a task complete without tests that demonstrate correctness.
