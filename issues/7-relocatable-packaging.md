# Issue 7: Package and validate relocatable binaries

**Status:** Active
**Parent:** [Issue 1](1-binary-distribution.md)
**Dependencies:** Approved release contract under issue 1;
[2](../closed-issues/2-canonical-versioning.md),
[3](../closed-issues/3-installed-artifacts.md)

## Problem

Installed backend artifacts need reproducible archive packaging and portable
runtime dependencies before CI can deliver usable release candidates.

## Current Linux evidence

Local issue L8 validated the versioned backend on GitHub's Ubuntu 22.04
x86_64 runner using GNU Fortran 11.4, CMake 3.31.6, and OpenBLAS 0.3.20. The
fresh Release build passed all 30 CTest entries, including version reporting,
startup errors, release-input validation, and rebuilding after a disposable
version change. The installed shared library was `libaenet.so.2.0.4` with
SONAME `libaenet.so.2` and the expected symlink chain.

An experimental relocation bundled `libgfortran.so.5`, `libquadmath.so.0`,
and `libgcc_s.so.1` and set relative ELF RPATHs; OpenBLAS was linked
statically. With no compiler or `LD_LIBRARY_PATH`, a clean Ubuntu 22.04
container passed exact version queries, the C API smoke check, and the
generate/train/predict workflow. This demonstrates a viable runtime strategy,
not a final archive or approved dependency/notices implementation. Preserve
and productionize the behavior under this issue.

Evidence: [successful workflow
run](https://github.com/atomisticnet/aenet/actions/runs/35636177866) at
temporary source commit `75681eb`; committed log `922fe23` on
`codex/l8-linux-validation`. The first run exposed a CMake script-policy
portability defect, addressed separately by local issue L9.

## Production Linux packaging evidence

The maintained issue 7 entry points passed on GitHub's Ubuntu 22.04 x86_64
runner in [workflow run
35815467080](https://github.com/atomisticnet/aenet/actions/runs/35815467080).
The run used the repository's common build command with GNU Fortran and static
OpenBLAS, passed the full 31-entry CTest suite, bundled the GNU Fortran,
quadmath, and GCC support runtimes, and created the named archive plus SHA-256
sidecar. The production inspection rejected non-x86_64 objects, unexpected
dependencies, dynamic BLAS/LAPACK, and non-relative runtime search paths.

The final gate validated the compressed archive itself in a fresh Ubuntu
22.04 container with no Fortran compiler and no `LD_LIBRARY_PATH`. All bundled
dependencies resolved from the extracted `lib/` directory. Exact CLI version
checks, the prebuilt native C API test, and the generate/train/predict numerical
smoke workflow passed. The temporary test branch and workflow were deleted
after the run. This establishes the Linux implementation; the candidate is not
publishable until local issue L15 adds and validates the complete dependency
notices and paired-platform metadata.

## Acceptance criteria

- Implement repeatable packaging for both platforms using the approved
  runtime strategy and issue 3's layout, without build-machine-specific paths.
- Include dependency/license notices, version/platform identification, and
  archive checksums; preserve required symlinks and executable permissions.
- Add automated extracted-archive validation in a new prefix and independent
  runtime environment, covering architecture/files, generate/train/predict
  numerical smoke tests, and native-library loading/API behavior.
- Ensure test harnesses do not supply dependencies from the build toolchain.
- Provide documented build/package/validate entry points consumed unchanged
  by issue 4. Do not duplicate packaging logic in CI YAML.
- Demonstrate passing candidates for native macOS arm64 and Linux x86_64.

Issue 4 owns CI orchestration; issue 8 owns user installation guidance;
local release work owns publication. No publication is part of this issue.

## Feasibility handoff

See the [proposed release contract](../doc/binary-release-contract.md) for
experimental evidence and limitations. Use the demonstrated platform-specific
runtime strategy as the starting point, and validate actual archives and the
full minimum-OS contract. The
first macOS candidate uses GNU 14; GNU 16 C exports remain unverified.

The extended feasibility report proposes system Accelerate on macOS and
OpenBLAS on Linux. Both macOS variants passed runtime checks; the proposal
simplifies dependencies and does not claim a benchmarked performance gain.
