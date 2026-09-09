---
name: issue-workflow
description: Plan, implement, review, validate, and close aenet work tracked by ROADMAP.md or GitHub issues. Use for addressing tracked work, not read-only exploration or triage.
---

# Issue Workflow

Follow AGENTS.md for approval, task storage, and external-action permissions.
Read [shared engineering standards](../references/engineering-standards.md)
before implementation or review. The current tracker is ROADMAP.md; do not
import aenet-python's global/local issue migration implicitly.

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

Establish a failing test or reproducible check before changing behavior when
practical. Implement one approved logical unit, updating relevant tests,
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

Archive a resolved roadmap task only after validation, when review-ready,
using AGENTS.md's location and status rules. The receipt records outcome,
commands/results, limitations, related tasks, and follow-up work. Remove the
resolved entry after archiving; do not require a commit hash in advance.
Leave work open if a promised deliverable or required validation is missing.
Only update public issue records with appropriate user authorization.

Summarize the result and evidence, identify remaining work, and propose a
commit message referencing the relevant task or explicitly identified GitHub
issue. Do not commit without confirmation.
