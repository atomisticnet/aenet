# Issue 2: Revisit canonical versioning

**Status:** Done
**Parent:** [Issue 1](../issues/1-binary-distribution.md)
**Legacy ID:** REL001

## Problem

Release metadata and executable version queries need a minimal coherent
policy for the first macOS arm64 and Linux x86_64 serial binary release.
The previous uncommitted attempt was reverted; the feature is not complete.
That attempt also introduced uncoordinated MPI shutdown on input errors.

## Acceptance criteria

- Agree canonical version, tag, and shared-library ABI conventions.
- Add automated checks for executable version reporting.
- Preserve coordinated MPI shutdown on missing arguments and input files.
- Validate affected build paths and update maintained documentation sources.

## Scope and order

Follow the release contract required by
[issue 1](../issues/1-binary-distribution.md). Agree the release version/tag
and library metadata needed for packaging; elaborate
pre-release/development-version machinery is deferred unless a concrete
release requirement needs it. Do not reinstate the reverted patch wholesale.

This is the first implementation unit after initial platform feasibility.
It precedes final archive naming and publication in issues 3 and 4. Shared
startup changes still require relevant MPI regression checks even though
MPI binaries are outside the first release.

## Resolution

`src/VERSION` remains the canonical plain `MAJOR.MINOR.PATCH` value, and Git
release tags use `vMAJOR.MINOR.PATCH`. CMake validates and consumes that value,
reconfigures when it changes, and applies the full version and major ABI
version to the shared library. The three main executables now provide exact,
side-effect-free `--version` output. Legacy Makefiles use the same generated
version module without changing their existing artifact names.

Startup validation now exits nonzero before application resources are created.
For MPI builds, rank zero decides version and input-error paths, broadcasts the
decision, and all ranks finalize together. The CMake MPI option now defines the
existing `PARALLEL` preprocessor path so its MPI linkage exercises actual MPI
behavior. Release preparation validates its version argument before editing
files and prints the corresponding prefixed tag commands.

## Validation

- GNU Fortran 14.2, Release, macOS arm64, Accelerate: complete 30-test serial
  suite passed, including a disposable version-change rebuild and Mach-O
  compatibility/current-version inspection.
- GNU Fortran 14.2, Release, macOS arm64, OpenBLAS 0.3.29: complete 30-test
  suite passed, including the same rebuild/metadata integration check.
- GNU Fortran 14.2 with Open MPI 5.0.7 and Accelerate: complete 35-test suite
  passed; two-rank version, no-argument, and nonexistent-input tests passed
  with bounded execution for `train.x` and `predict.x`.
- A normal generate/train/predict smoke workflow and a C API initialization,
  atom-type conversion, and finalization check passed with Accelerate.
- The GNU/macOS legacy Makefile built the version module and reported
  `generate.x 2.0.4` from its existing suffixed executable.
- `prepare-release.sh` syntax and disposable valid/invalid version fixtures
  passed. `git diff --check` passed.

Intel Fortran and a current Linux build were unavailable locally and were not
rerun for this change. The implementation uses standard CMake and Fortran 2008
features, but those configurations remain validation gaps. The maintained Org
manual source was updated. Emacs 29.4 introduced broad unrelated formatting
churn in generated text output, and Pandoc was unavailable on the host PATH,
so generated manual exports were intentionally left unchanged.
