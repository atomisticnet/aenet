# Issue 4: Automate GNU serial binary releases

**Status:** Pending
**Legacy ID:** CI001

## Problem

Supported platforms need reproducible builds and downloadable release assets.

## Acceptance criteria

- Add Linux and macOS GNU serial CI builds using CMake for main executables,
  tools, and the shared library.
- Run agreed release-gating tests.
- Package archives with checksums and publish assets for tagged releases.

## Dependencies

Agree [versioning](2-canonical-versioning.md) and
[installed layout](3-installed-artifacts.md) before finalizing packaging.
Publication still requires the applicable authorization.
