# Issue 2: Revisit canonical versioning

**Status:** Deferred
**Legacy ID:** REL001

## Problem

Release metadata and executable version queries need a coherent policy.
The previous uncommitted attempt was reverted; the feature is not complete.
That attempt also introduced uncoordinated MPI shutdown on input errors.

## Acceptance criteria

- Agree canonical version, tag, and shared-library ABI conventions.
- Add automated checks for executable version reporting.
- Preserve coordinated MPI shutdown on missing arguments and input files.
- Validate affected build paths and update maintained documentation sources.
