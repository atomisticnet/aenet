# Issue 4: Automate GNU serial binary releases

**Status:** Done
**Legacy ID:** CI001
**Parent:** [Issue 1](../issues/1-binary-distribution.md)
**Dependency:** [Issue 7](../closed-issues/7-relocatable-packaging.md)

## Problem

The agreed build, packaging, and runtime checks must run reproducibly for
both supported platforms and produce reviewable release candidates.

## Acceptance criteria

- Add CMake-based GNU serial CI jobs for native macOS arm64 and Linux x86_64
  using the environments approved in the release contract.
- Invoke issue 7's packaging and validation entry points; include relevant
  CTest checks and independent runtime environments for extracted artifacts.
- Retain validated archives, checksums, and sufficient build/test provenance
  as candidate artifacts for review; report failures before publication.
- Implement the authorized tagged-release publication path so it promotes
  tested artifacts rather than rebuilding different binaries without testing.
- Exercise candidate generation without publishing a release and document
  the release trigger/approval process for local release execution.

Platform feasibility is tracked locally under issue 1, packaging/tests by issue 7, and
user documentation by issue 8. Actual first-release publication and public
asset verification are tracked locally under issue 1; workflow development does not itself
authorize a push, remote dispatch, or publication.

## Feasibility handoff

See the [proposed release contract](../doc/binary-release-contract.md) for
experimental evidence and limitations. Before enabling parallel builds,
resolve the observed static/shared Fortran module-output race or retain single-job builds as an explicit limitation.

## Resolution

Added commit-pinned GitHub Actions workflows that build GNU serial candidates
on Ubuntu 22.04 x86_64 and macOS 15 arm64, invoke the maintained packaging and
independent runtime checks, and retain a verified paired candidate set. The
manual release workflow requires an existing annotated version tag, matching
source and artifact metadata, an unpublished tag, and approval through the
protected `release` environment. Only its final publication job receives
`contents: write`, and that job downloads the already validated artifacts
from the same run instead of rebuilding them.

The `release` environment requires review by `alexurba`, disallows
administrator bypass, and accepts deployments only from `master`.

## Validation

Thirty packaging tests pass locally, including candidate-set and release-tag
regressions. Candidate run 36039685568 passed both native build, CTest,
packaging, independent-runtime, and paired-artifact jobs. Full dry run
36042038039 additionally passed the annotated-tag, absence-of-release, and
prepublication gates; its publication job was skipped and no release was
created. The temporary validation branch and tag were deleted. An earlier
unchanged Linux attempt intermittently failed `symmfunc:derivatives`; its
successful rerun passed all 31 CTest entries, and L18 tracks that test issue.
