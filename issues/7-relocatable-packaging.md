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
