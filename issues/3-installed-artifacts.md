# Issue 3: Normalize installed artifact naming and layout

**Status:** Pending
**Parent:** [Issue 1](1-binary-distribution.md)
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

Use the release contract required by [issue 1](1-binary-distribution.md) and
[versioning policy](2-canonical-versioning.md). Settle runtime-library placement
with [packaging](7-relocatable-packaging.md) so extracted archives can be
relocated. Issue 7 owns dependency bundling and clean-environment execution
checks. Do not add developer suffix options without a demonstrated need.
