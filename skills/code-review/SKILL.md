---
name: code-review
description: Review aenet commits, branches, diffs, or working-tree changes for correctness, numerical validity, compatibility, maintainability, and test quality. Reviews are read-only unless fixes are requested.
---

# Code Review

Follow AGENTS.md and read
[shared engineering standards](../references/engineering-standards.md).

Resolve the requested review boundary: commits, branch relative to merge
base, staged/unstaged changes, or another explicit range. Include requested
refinement commits. Read the complete scoped diff, governing tasks, relevant
neighboring code, callers, tests, documentation, and build configuration.
Do not silently reduce a history review to the final snapshot.

Apply the [review checklist](references/review-checklist.md). Prioritize
observable failures, numerical validity, API/ABI and format compatibility,
MPI behavior, and meaningful test evidence. Evaluate simplicity: each new
abstraction, option, fallback, or dependency should serve an actual contract.
Do not report stylistic preferences as defects without concrete impact.

Use the [build-test skill](../build-test/SKILL.md) for relevant validation.
Run checks in isolated build directories without changing reviewed sources
or tests. Distinguish unavailable configurations from failed validation.
Consult the [documentation skill](../documentation/SKILL.md) when maintained
documentation or public interface comments are materially affected.

Report actionable findings first, ordered by severity. Give a precise file
and line location, triggering scenario, and consequence for each finding:

- P0: immediate catastrophic or security-critical impact.
- P1: release-blocking correctness, data-loss, or major compatibility bug.
- P2: substantive defect, missing requirement, or concrete regression risk.
- P3: low-risk improvement with a demonstrated maintenance or usability benefit.

Separate findings from questions and assumptions. Summarize checks and
remaining coverage gaps. State explicitly when there are no actionable
findings; do not manufacture minor issues. Do not implement fixes or add a
fix plan unless requested.
