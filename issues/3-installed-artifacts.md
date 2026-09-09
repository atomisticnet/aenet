# Issue 3: Normalize installed artifact naming and layout

**Status:** Pending
**Legacy ID:** REL002

## Problem

External consumers need stable installed names and a predictable layout;
build-variant information should be represented in release archive names.

## Acceptance criteria

- Install main executables as generate.x, train.x, and predict.x.
- Keep tool names stable under tools/ and install aenet.h with native libraries.
- Define GNU serial archive names for Linux and macOS.
- Decide whether legacy suffixes remain an optional developer build setting.
- Validate the installed layout for discovery by companion tools.
