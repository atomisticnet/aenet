# Issue 4: Automate GNU serial binary releases

**Status:** Pending
**Legacy ID:** CI001
**Parent:** [Issue 1](1-binary-distribution.md)
**Dependency:** [Issue 7](7-relocatable-packaging.md)

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
