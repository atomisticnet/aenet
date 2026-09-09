# Issue 8: Document binary installation

**Status:** Pending
**Parent:** [Issue 1](1-binary-distribution.md)
**Dependency:** [Issue 7](7-relocatable-packaging.md)

## Problem

Users need a short, tested path from selecting the correct archive to running
an example, with answers to common platform/runtime questions.

## Acceptance criteria

- Update maintained documentation sources with architecture/OS support,
  download selection, checksum verification, extraction, PATH/library usage,
  a small runnable example, and troubleshooting.
- Exercise instructions against issue 7's candidates; use actual filenames,
  dependencies, and minimum OS claims supported by the approved release-contract evidence.
- Assess macOS downloaded-artifact security behavior using an authorized
  candidate download when available. Record any signing/notarization work as
  an explicit issue if required; do not infer it from local execution alone.
- Clearly identify any public-download steps awaiting verification during local release execution.
  Draft instructions can use the agreed release URL convention until then.
- Regenerate and inspect affected maintained exports under the documentation
  skill. No automatic Python downloader is required.

Local release execution verifies the final public links and downloaded user workflow.
