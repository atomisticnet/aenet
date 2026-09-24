# Issue 8: Document binary installation

**Status:** Done
**Parent:** [Issue 1](../issues/1-binary-distribution.md)
**Dependency:** [Issue 7](../closed-issues/7-relocatable-packaging.md)

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

## Current evidence

The modular binary-installation source and maintained manual exports document
platform selection, checksums, extraction, paths, a runnable generate example,
and troubleshooting. The documented checksum, extraction, version, and
generate commands pass on the packaged macOS candidate after security-policy
isolation.

A Firefox download of candidate run 36042038039 carried macOS quarantine and
GitHub provenance attributes through the artifact ZIP, tar archive, and
extracted Mach-O files. Strict code-signature verification passed, but
Gatekeeper rejected the ad-hoc signatures and all three main executables were
killed with signal 9. After checksum verification, removing quarantine from
the downloaded tar archive before extraction produced a tree without
quarantined Mach-O files; all version checks and the generate example passed.
Approving =generate.x= through System Settings allowed that executable to
start but did not approve the quarantined GNU runtime libraries; the dynamic
loader then rejected =libgfortran.5.dylib=. The other main executables also
remained separately blocked, so per-executable approval is not the documented
installation path.
The installation page documents this bounded approval step. Issue 10 tracks
Developer ID signing/notarization as a deferred usability improvement.

## Resolution

Added a modular binary-installation page to the maintained Org manual and
regenerated its Markdown, text, and PDF exports. The page documents platform
selection, checksum verification, the checksum-first macOS quarantine
approval step, extraction, paths, a runnable example, and troubleshooting.
Issue 10 records the deferred signing, notarization, and graphical-installer
options.

## Validation

- Verified checksums, extraction, all three version commands, and the runnable
  generate example with the packaged macOS candidate.
- Repeated the checks with a Firefox-downloaded hosted candidate. Confirmed
  that approving one executable in System Settings did not authorize its GNU
  runtime libraries, while removing quarantine from the verified archive
  before extraction produced a working installation.
- Inspected the regenerated PDF pages containing the complete installation
  guide. The Linux archive contract and runtime behavior retain the hosted
  Ubuntu validation recorded by issues 4 and 7.
