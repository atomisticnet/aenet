---
name: issue-workflow
description: Plan, implement, review, validate, and close aenet work tracked by shared or local issues, or explicitly linked GitHub issues. Use for addressing tracked work, not read-only exploration or triage.
---

# Issue Workflow

Follow AGENTS.md for approval, task storage, and external-action permissions.
Read [shared engineering standards](../references/engineering-standards.md)
before implementation or review. Use ISSUES.md for shared work and
LOCAL_ISSUES.md for implementation tasks; a small task needs only a local
issue. Keep records proportional to the work.

## Establish the contract and plan

Read the complete task, related records, relevant implementation, tests,
public callers, and documentation. Verify that the task still matches actual
behavior. Identify acceptance criteria, compatibility constraints, numerical
assumptions, required configurations, dependencies, and out-of-scope work.

Split work only when it contains independently reviewable units or requires
separate contracts or validation. Record dependencies and keep an umbrella
open until its children and combined acceptance criteria are complete.

Present the behavioral outcome, implementation units, regression checks,
compatibility decisions, documentation changes, and validation plan. During
planning, write only private development notes; leave code, tests,
maintained documentation, and task status unchanged. Obtain approval under
AGENTS.md before implementation. An approved detailed child sequence does
not need repeated approval unless its scope or contracts materially change.

## Implement and review

Apply the shared engineering standards' test-first policy, including
characterization coverage and documented exceptions. Implement one approved
logical unit, updating relevant tests,
interface comments, documentation, and task evidence together. Use the
[build-test skill](../build-test/SKILL.md) for backend checks and the
[documentation skill](../documentation/SKILL.md) for substantive docs work.

Before declaring the unit ready to commit, inspect its complete diff and
apply the [review checklist](../code-review/references/review-checklist.md).
Passing tests do not substitute for reviewing numerical correctness,
compatibility, simplicity, and test integrity. Preserve unrelated edits.

Resolve known P0-P2 findings within scope. Record proportionate P3 follow-ups.
When a required fix materially expands the approved contract, describe the
finding and obtain direction rather than silently broadening the task.

## Validate, close, and hand off

Run focused checks, then broaden as regression risk requires. Compare every
acceptance criterion with the final implementation and actual validation.
Distinguish environment failures, product failures, and unrun checks.

Close local work after validation, when ready to commit, by moving its entry
into CLOSED_LOCAL_ISSUES.md with a concise outcome and validation receipt.
For completed shared work, move its issue file to closed-issues/ and remove
its active index entry before merging, following AGENTS.md. Neither closure
requires a final commit hash or a date in the filename. Preserve ID counters.
Leave deferred issues indexed and keep an issue open when a promised
deliverable or required validation is missing. Do not close a parent merely
because one child is complete. Only update GitHub with user authorization.

Summarize the result and evidence, identify remaining work, and propose a
commit message referencing the relevant task or explicitly identified GitHub
issue. Do not commit without confirmation.
