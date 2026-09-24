# Issue 3: Normalize installed artifact naming and layout

**Status:** Done
**Parent:** [Issue 1](../issues/1-binary-distribution.md)
**Legacy ID:** REL002

## Problem

External consumers need stable installed names and a predictable layout;
build-variant information should be represented in release archive names.

## Acceptance criteria

- Install main executables as generate.x, train.x, and predict.x.
- Keep tool names stable under tools/ and install aenet.h with native libraries.
- Define versioned GNU serial archive names identifying macOS arm64 or
  Linux x86_64; record the agreed minimum OS/runtime requirements.
- Decide whether legacy suffixes remain an optional developer build setting.
- Validate installation into a temporary prefix, stable names, included
  tools/header/libraries, and discovery by companion tools.

## Dependencies and boundary

Use the release contract required by
[issue 1](../issues/1-binary-distribution.md) and
[versioning policy](2-canonical-versioning.md). Settle runtime-library
placement with [packaging](7-relocatable-packaging.md) so extracted
archives can be relocated. Issue 7 owns dependency bundling and clean-
environment execution checks. Do not add developer suffix options without a
demonstrated need.

## Resolution

CMake now installs the main executables as `generate.x`, `train.x`, and
`predict.x`, the existing utilities under `tools/`, the public C header under
`include/`, and the static and versioned shared libraries under `lib/`.
Developer build trees retain their variant suffixes; no additional suffix
option was introduced.

The release contract fixes the archive names as
`aenet-<version>-macos-arm64-gnu-serial.tar.gz` and
`aenet-<version>-linux-x86_64-gnu-serial.tar.gz`. The first-release runtime
baselines are macOS 14 or newer on arm64 and Ubuntu 22.04/glibc 2.35 or newer
on x86_64. Issue 7 still owns runtime bundling and final extracted-archive
validation.

## Validation

The install-layout test first failed against the previous suffixed install.
After implementation, a local macOS GNU 13.2/Accelerate Release build installed
the exact layout into an isolated prefix, checked all three CLI versions, and
compiled and ran a minimal C consumer using only the installed header and
shared library. The focused test and 29 other CTest entries passed.
`version:rebuild` exceeded its five-minute harness timeout while still
compiling; its direct rerun completed and produced the expected 7.8.9
executables and shared-library metadata. The companion `aenet-python`
configuration command's actual discovery logic found all three main
executables, `trnset2ASCII.x`, and the shared library in the installed tree.

[Hosted run 35775370837](https://github.com/atomisticnet/aenet/actions/runs/35775370837)
then passed the complete configure, serial build, and CTest step on Ubuntu
22.04 with GNU Fortran and OpenBLAS. The workflow existed only on the temporary
`test/issue-3-linux-install` branch, was not added to the product history, and
the branch was deleted after the result was recorded.
