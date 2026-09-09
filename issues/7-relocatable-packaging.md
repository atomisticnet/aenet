# Issue 7: Package and validate relocatable binaries

**Status:** Pending
**Parent:** [Issue 1](1-binary-distribution.md)
**Dependencies:** Approved release contract under issue 1;
[2](2-canonical-versioning.md), [3](3-installed-artifacts.md)

## Problem

Installed backend artifacts need reproducible archive packaging and portable
runtime dependencies before CI can deliver usable release candidates.

## Acceptance criteria

- Implement repeatable packaging for both platforms using the approved
  runtime strategy and issue 3's layout, without build-machine-specific paths.
- Include dependency/license notices, version/platform identification, and
  archive checksums; preserve required symlinks and executable permissions.
- Add automated extracted-archive validation in a new prefix and independent
  runtime environment, covering architecture/files, generate/train/predict
  numerical smoke tests, and native-library loading/API behavior.
- Ensure test harnesses do not supply dependencies from the build toolchain.
- Provide documented build/package/validate entry points consumed unchanged
  by issue 4. Do not duplicate packaging logic in CI YAML.
- Demonstrate passing candidates for native macOS arm64 and Linux x86_64.

Issue 4 owns CI orchestration; issue 8 owns user installation guidance;
local release work owns publication. No publication is part of this issue.
