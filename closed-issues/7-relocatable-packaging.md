# Issue 7: Package and validate relocatable binaries

**Status:** Done
**Parent:** [Issue 1](../issues/1-binary-distribution.md)
**Dependencies:** Approved release contract under issue 1;
[2](2-canonical-versioning.md), [3](3-installed-artifacts.md)

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
after the run. This established the Linux runtime implementation before L15
added and validated the common notices and metadata contract.

## Production macOS packaging evidence

The maintained issue 7 entry points passed on GitHub's native macOS 14 arm64
runner in [workflow run
36026156581](https://github.com/atomisticnet/aenet/actions/runs/36026156581).
The run used GNU Fortran 14 with system Accelerate, passed the full 31-entry
CTest suite, bundled the GNU Fortran, quadmath, and GCC support runtimes, and
created the named archive plus SHA-256 sidecar. Preparation rewrote Mach-O
dependencies to relative loader paths, set relative library IDs, applied
ad-hoc signatures after modification, and rejected non-arm64 objects,
deployment targets newer than macOS 14.0, OpenBLAS/OpenMP, and unresolved or
toolchain dependencies.

The final gate validated the compressed archive after extraction with common
Homebrew and active Xcode developer paths denied by `sandbox-exec`. Exact CLI
version checks, the separately built native C API test, and the
generate/train/predict numerical smoke workflow passed. The temporary test
branch and workflow were deleted after the run. L15 subsequently completed
the redistribution notices and paired-platform metadata.

## Paired final validation

[Workflow run
36027659670](https://github.com/atomisticnet/aenet/actions/runs/36027659670)
validated the final contract from one source revision on both Ubuntu 22.04
x86_64 and native macOS 14 arm64. Both jobs passed 21 packaging regressions,
the full 31-entry CTest suite, platform runtime preparation, deterministic
metadata and notice installation, archive and checksum creation, and the
independent extracted-archive gate. The Linux container had no Fortran
compiler or `LD_LIBRARY_PATH`; the macOS smoke tests denied the common
Homebrew and active Xcode developer paths. Both native C API tests and
generate/train/predict numerical workflows passed. The temporary branch and
workflow were deleted after recording this evidence.

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

## Resolution

The repository now provides maintained commands to build and stage GNU serial
release artifacts, relocate their platform-specific runtime dependencies,
install deterministic metadata and required redistribution notices, create
reproducible archives and checksums, and validate the extracted archives in
independent runtime environments. Linux uses static OpenBLAS with bundled GNU
runtimes; macOS uses system Accelerate with bundled GNU runtimes and ad-hoc
signatures. The archive metadata and hashed manifest identify every candidate
and its payload without build-machine paths.

## Validation

Focused implementation evidence is recorded above for Linux run 35815467080,
macOS run 36026156581, and final paired run 36027659670. The paired run is the
completion gate because it exercised the notice and metadata contract on both
platforms from the same revision. No archive was published. Issue 4 owns
permanent release CI, issue 8 owns installation guidance, and local release
work owns publication and downloaded-artifact verification.
