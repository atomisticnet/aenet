# Issue 2: Revisit canonical versioning

**Status:** Pending
**Parent:** [Issue 1](1-binary-distribution.md)
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

Follow the release contract required by [issue 1](1-binary-distribution.md). Agree the
release version/tag and library metadata needed for packaging; elaborate
pre-release/development-version machinery is deferred unless a concrete
release requirement needs it. Do not reinstate the reverted patch wholesale.

This is the first implementation unit after initial platform feasibility.
It precedes final archive naming and publication in issues 3 and 4. Shared
startup changes still require relevant MPI regression checks even though
MPI binaries are outside the first release.
