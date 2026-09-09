# Code Review Checklist

Apply relevant items to the requested scope; do not manufacture findings.

## Contract and design

- The review includes all requested changes and governing acceptance criteria.
- The implementation neither silently narrows nor expands the agreed scope.
- Control flow and procedure responsibilities are understandable in context.
- New helpers, options, validation, and dependencies serve current needs.
- Remove unnecessary concepts when correctness, compatibility, and clarity
  would be preserved; do not generalize for hypothetical future needs.
- No unrelated edits, dead paths, duplicated logic, or accidental artifacts.

## Numerical and interface correctness

- Units, precision, normalization, indexing, periodic images, and gradient
  signs remain correct, including degenerate and boundary cases.
- Analytical, derivative, or invariance checks provide independent evidence
  where applicable. Tolerances have numerical justification.
- Array dimensions, bounds, initialization, allocation lifetimes, and aliases
  are valid; failure paths do not hide errors or leak resources.
- C bindings preserve interoperable types, layout, strings, and ownership.
- CLI, input/output formats, persisted models, and installed names remain
  compatible or have an intentional documented migration.
- MPI paths agree on collective ordering and failure handling; reductions
  and reproducibility assumptions match the supported contract.
- Compiler features, preprocessing, BLAS/LAPACK kinds, and link requirements
  work for affected supported configurations.

## Test integrity and evidence

- Changed functionality has meaningful unit or integration coverage.
- Regressions exercise the production path and would fail without the fix.
- Assertions establish outcomes, not merely successful execution or internal
  self-consistency. Tests do not duplicate the algorithm as their oracle.
- Skips, retries, conditional assertions, fallbacks, and relaxed tolerances
  do not conceal defects. Randomness and output locations are controlled.
- Parallel tests cannot race on shared files; MPI runs use appropriate rank
  counts and have bounded execution when hangs are possible.
- Relevant CMake/CTest checks ran, with compiler, build type, and backend
  choices recorded. Unavailable configurations are stated as gaps.

## Documentation and final report

- Procedure help, CLI documentation, examples, and format descriptions match
  the implemented contract and identify consequential units and assumptions.
- Maintained sources are edited, not just their generated exports.
- New dependencies justify portability, ABI, installation, and license costs.
- Private paths, backup files, and generated binaries do not enter the diff.
- Findings identify concrete scenarios and impact with precise locations.
- Task closure and final claims match actual review and validation evidence.
