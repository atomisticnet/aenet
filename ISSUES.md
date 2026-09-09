# Shared Issues

Repository issue IDs are independent of GitHub issue numbers.
Details live in issues/; completed issue files move to closed-issues/.

**Last issue ID assigned:** 9

| ID | Status | Issue |
| --- | --- | --- |
| 1 | Active | [Umbrella: first backend binary release](issues/1-binary-distribution.md) |
| 2 | Pending | [Revisit canonical versioning](issues/2-canonical-versioning.md) |
| 3 | Pending | [Normalize installed artifacts](issues/3-installed-artifacts.md) |
| 4 | Pending | [Automate GNU serial binary releases](issues/4-binary-release-automation.md) |
| 5 | Deferred | [Integrate binaries with aenet-python](issues/5-python-binary-integration.md) |
| 7 | Pending | [Package and validate relocatable binaries](issues/7-relocatable-packaging.md) |
| 8 | Pending | [Document binary installation](issues/8-binary-installation-docs.md) |

Issue 1 is an umbrella only. Start with local platform-feasibility work,
then implement shared issues 2 → 3 → 7 → 4 → 8, followed by local release
publication and verification. Issue 8 may start after issue 7.
Issue 5 remains deferred and does not block the first-release milestone.

Retired IDs: [6](closed-issues/6-platform-feasibility.md) and
[9](closed-issues/9-first-binary-release.md) moved to local tracking.
Their work remains outstanding; IDs will not be reused.
