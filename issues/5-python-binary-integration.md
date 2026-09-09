# Issue 5: Integrate GitHub binaries with aenet-python

**Status:** Deferred
**Legacy ID:** PY001

## Problem

The companion Python project needs a convenient way to discover or explicitly
install compatible GitHub-hosted backend binaries. Implementation belongs in
that project; this issue tracks the backend distribution dependency.

## Acceptance criteria

- Prefer discovery of an existing local installation.
- Provide an explicit download/install command, with no implicit network
  activity during pip install.
- Verify checksums, unpack into a managed location, and configure executable
  and library paths without manual entry.
- Coordinate with [binary release automation](4-binary-release-automation.md)
  and record the companion implementation or an explicit scope decision.
