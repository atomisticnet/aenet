# Issue 10: Improve macOS binary download security

**Status:** Deferred
**Related:** [Issue 1](1-binary-distribution.md)
**Dependencies:** [Issue 4](../closed-issues/4-binary-release-automation.md),
[Issue 7](../closed-issues/7-relocatable-packaging.md)

## Problem

The macOS archive is ad-hoc signed after runtime relocation. Files extracted
from a Firefox-downloaded candidate retain quarantine metadata, and Gatekeeper
rejects all three main executables. Users would have to approve or remove
quarantine manually. The first release may document removal of quarantine
from the checksum-verified archive before extraction, but Developer ID signing
and notarization would provide a smoother and more conventional experience.
System Settings approval of one executable does not approve the quarantined
GNU runtime libraries or the other executables.

A small macOS launcher or installer may provide a familiar graphical entry
point, but it must cover the complete runtime dependency tree rather than
merely launching each command once. Evaluate an app bundle, disk image, or
flat installer package alongside the current archive. Apple supports
notarization of all three container types; a signed and notarized installer
package is likely to fit command-line tools better than a launcher whose only
purpose is to provoke separate Gatekeeper prompts.

This route requires Apple Developer Program membership, currently listed by
Apple at USD 99 per membership year with regional pricing. Eligible nonprofit
organizations, accredited educational institutions, and government entities
may request a fee waiver. Apple does not list a separate per-release
notarization charge.

## Acceptance criteria

- Select a secure Developer ID signing and Apple notarization procedure for
  every executable and dynamic library in the final relocated tree, after all
  load-command changes.
- Compare the archive with a minimal graphical installer or launcher. If a
  GUI is retained, verify that one deliberate user action authorizes the
  executables and bundled runtime libraries and that installation paths and
  upgrades remain predictable.
- Store signing/notarization credentials in protected GitHub configuration;
  do not expose them to pull requests or read-only candidate workflows.
- Preserve the relative runtime dependency closure and strict code-signature
  validity. If notarization requires a different release container, update the
  archive contract, automation, checksums, and installation documentation
  deliberately.
- When signed releases are enabled, make signing and notarization failures
  stop publication with useful diagnostics before release creation.
- Download the resulting candidate through a normal macOS browser and verify
  that quarantine propagates but Gatekeeper permits the documented version
  checks and runnable example without manual security exceptions.
- Record credential prerequisites, validation evidence, and any remaining
  online-verification or minimum-OS limitations. Do not publish a release
  without the separate authorization required by the release workflow.

## References

- [Apple Developer Program enrollment and fees](https://developer.apple.com/programs/enroll/)
- [Developer ID certificates](https://developer.apple.com/help/account/certificates/create-developer-id-certificates)
- [Notarizing macOS software before distribution](https://developer.apple.com/documentation/security/notarizing-macos-software-before-distribution)
