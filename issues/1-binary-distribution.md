# Issue 1: Umbrella — first backend binary release

**Status:** Active
**Type:** Umbrella (no implementation work)
**Legacy ID:** DIST001

## Problem

Missing ready-to-run Apple silicon macOS and Linux backend binaries cause
installation difficulties and repeated user questions. This issue coordinates
the first usable GitHub binary release. Maintained deliverables belong to
shared children; investigations and release operations use local issues.

## First-release scope

- Native macOS arm64 and Linux x86_64 archives, built with GNU Fortran in
  serial mode. Minimum macOS and Linux compatibility baselines must be
  selected and documented before the release contract is finalized.
- Main executables, tools, native libraries, and C header in a stable layout.
- Users can download, extract, and run a documented example without installing
  a Fortran compiler or relying on the build machine's library paths.
- A documented runtime dependency strategy, checksums, and required license
  notices accompany each archive.
- MPI/Intel binaries, Linux arm64, Intel macOS binaries, Conda/PyPI binary
  distribution, and automatic Python downloads are outside this milestone.
  Existing supported source-build configurations must not regress.

## Child issues and ownership

| Issue | Deliverable | Depends on |
| --- | --- | --- |
| [2](2-canonical-versioning.md) | Minimal versioning and regression checks | Release contract |
| [3](3-installed-artifacts.md) | Stable installed names and layout | Release contract, 2 |
| [7](7-relocatable-packaging.md) | Reproducible archives and runtime validation | Release contract, 2, 3 |
| [4](4-binary-release-automation.md) | CI orchestration and release-candidate artifacts | 7 |
| [8](8-binary-installation-docs.md) | Tested installation/troubleshooting guidance | 7 |

[Issue 5](5-python-binary-integration.md) is a deferred follow-up, not a
first-release dependency. New concrete work discovered during delivery must
be assigned to an appropriate shared or local issue before implementation.

## Acceptance criteria

- Shared children 2, 3, 7, 4, and 8 are complete against their criteria.
- Feasibility evidence and the approved OS/runtime, archive, and validation
  contract are recorded in maintained documentation before implementation
  relies on them. This investigation is tracked locally.
- Both validated archives are published with authorization; public checksums,
  links, and the downloaded user workflow are verified. Release execution is
  tracked locally and remains required for umbrella closure.
- The approved scope works end to end: download, extract, and run documented
  examples without a Fortran compiler or build-machine library paths.
- Any change in milestone scope is explicit and reflected in the child
  issues; issue 5 remains independently deferred.

The umbrella remains open while any required child deliverable is missing.
Its first practical task is local platform feasibility; approval of this
issue structure is not
approval to execute every child or publish a release.
