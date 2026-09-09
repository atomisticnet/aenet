# Shared Engineering Standards

Apply relevant standards during implementation and review. These guide
judgment; they do not expand the approved task into unrelated cleanup.

## Numerical correctness and compatibility

- Preserve precision, units, normalization, indexing, periodic-boundary
  conventions, and force/gradient signs. Document consequential assumptions.
- Treat Fortran interfaces, C bindings, CLI options and output consumed by
  tools, input keywords, model/training-data formats, and installed names as
  compatibility surfaces. Describe deliberate changes and migration needs.
- At C boundaries, check interoperable kinds, shapes and array layout,
  string conventions, allocation ownership, lifetimes, and error reporting.
- Check numerical routines against meaningful independent evidence where
  applicable: analytical cases, finite-difference derivatives, and physical
  invariants such as rotation, translation, or permutation invariance.
- Justify absolute/relative tolerances using scale, conditioning, and expected
  floating-point error. Do not loosen tolerances or add fallback results to
  conceal production defects. Do not promise bitwise agreement across
  compilers, BLAS libraries, or MPI reductions without evidence.
- For MPI changes, reason about collective ordering, rank-dependent paths,
  reductions, initialization/finalization, and failures that could hang peers.
- Check BLAS/LAPACK integer kinds, interfaces, and linkage when affected.
  Do not assume LP64 and ILP64 libraries are interchangeable.

## Fortran style and design

- Follow surrounding Fortran conventions; prefer local consistency over
  restyling. Much of the repository follows Fortran 90-era conventions.
- Modern features are acceptable when supported by target GNU and Intel
  compilers from roughly the last five years; verify support when uncertain.
- Avoid new goto statements and line labels. Do not rewrite legacy control
  flow unless the task calls for it.
- Aim for 72 columns, keep source within 80 where practical, and allow minor
  comment overshoot for clarity. Respect language-specific line limits.
- Give procedures concise help text describing purpose, arguments, outputs,
  units, and consequential side effects or ownership. Explain non-obvious
  numerical rationale rather than restating syntax.
- Use the smallest direct design satisfying current requirements. Avoid
  speculative options, duplicated algorithms, test-only production branches,
  broad error suppression, and unrelated cleanup.
- Start new source files with the license text in src/license-header.txt,
  formatted as valid comments for the language. Preserve required shebangs
  and format directives. Do not prepend raw license text to structured data,
  Markdown, or skill YAML frontmatter; retain any existing license notices.

## Tests and dependencies

- Add focused unit or integration tests for changed behavior and meaningful
  failure paths. A regression check should fail without the fix and exercise
  the real supported path, not a copy of the implementation.
- Control random seeds where practical. Review skips, retries, expected
  failures, conditional assertions, and numerical tolerances critically.
- Keep tests independent and isolate generated files. Use fixtures and
  cleanup compatible with parallel execution before running tests in parallel.
- Prefer existing dependencies. Justify new libraries through concrete
  benefit and account for compiler/platform support, ABI, licensing,
  installation, transitive dependencies, and maintenance costs.
- Keep MPI and vendor BLAS choices optional for unrelated serial workflows.
  Support performance claims with representative measurements.

## Documentation, hygiene, and validation

- Update maintained sources when behavior, API, CLI, defaults, formats,
  units, dependencies, or supported workflows change. Qualify scientific and
  performance claims to the available evidence.
- Preserve unrelated working-tree changes and inspect the final scoped diff.
  Do not track private paths, credentials, backups, or accidental build output.
- Keep build options in CMake; local tool locations belong in ignored config.
  Do not invent configuration files without a defined consumer.
- Use the build-test skill for backend validation. Broaden checks in
  proportion to regression risk; unavailable checks are gaps, not passes.
