# Context

- This repository contains ænet's Fortran backend for machine-learned
  interatomic potentials: descriptors, training, prediction, command-line
  tools, and native libraries.
- The companion `aenet-python` project consumes backend executables and
  library interfaces. Preserve those contracts deliberately; this repository
  must not require a sibling checkout for its development policies.

# Development environment

- Use CMake for builds and CTest for tests. Maintain legacy build Makefiles
  only when the task requires them.
- Target GNU Fortran (`gfortran`) and Intel Fortran (`ifort`, `ifx`). Keep
  compiler flags and build defaults authoritative in CMake.
- Read `config/local.toml` when present for machine-specific tool commands;
  use `config/local.toml.example` as its template. Local configuration is
  untracked and is guidance for the agent, not a CMake input.
- Do not hard-code user paths or environment locations in shared policy.
- Current documentation is maintained under `doc/` using Org mode and its
  export Makefile. A Sphinx/Read the Docs migration is separate work.
- Release preparation uses `src/prepare-release.sh`; inspect it and the
  current version/build configuration before release work. Preparing a
  change does not authorize tagging, pushing, or publishing a release.

# Repository skills

Read each applicable `SKILL.md` completely before using its workflow:

- [issue-workflow](skills/issue-workflow/SKILL.md): implement and close work
  tracked by shared or local issues, or explicitly linked GitHub issues.
- [code-review](skills/code-review/SKILL.md): review commits, branches, or
  working-tree changes; includes the pre-commit review checklist.
- [build-test](skills/build-test/SKILL.md): configure, build, and validate
  backend changes with CMake and CTest.
- [documentation](skills/documentation/SKILL.md): substantially revise
  maintained documentation, interface comments, or runnable examples.

Use issue-workflow to coordinate tracked work and specialized skills for
its deliverables. Apply
[shared engineering standards](skills/references/engineering-standards.md)
during implementation and review.

# Development workflow

- Begin development with planning. Do not modify code, tests, or maintained
  documentation during planning; private notes may be created in dev-notes/.
- Obtain approval before implementation unless the user has already approved
  the proposed plan. Approval covers its described implementation and
  validation steps. Ask again when discoveries materially change scope,
  behavior, risk, or public contracts, not for routine implementation choices.
- Follow the practical test-first policy in the shared engineering standards.
  Add meaningful automated coverage for new or changed behavior; explain
  exceptions to test-first development and validate before completion.
- Keep changes scoped to the approved logical unit; preserve unrelated user
  edits. Ask when requirements are ambiguous or progress is blocked.
- Keep procedure comments, help text, and maintained documentation current
  with changed behavior and interfaces.
- Complete the applicable review and validation before declaring work done.
  Report unavailable checks and unresolved limitations explicitly.
- Summarize changes, testing, task status, and follow-up work. Leave handoff
  notes for phased work and propose a focused commit message.
- Do not commit without user confirmation. Do not infer authorization to
  send messages, update GitHub, push, merge, or publish from local task work.

# Shared issues

- Use tracked `ISSUES.md` as the concise index of active and deferred shared
  work. Store substantial descriptions in `issues/<id>-<description>.md`.
- Assign stable integer IDs. Retain a last-assigned-ID counter in ISSUES.md;
  never reuse IDs. Categories belong in titles or optional metadata.
- Require only a problem, acceptance criteria, and status. Add plans,
  dependencies, and evidence when useful. Small work may live entirely in
  the local tracker; do not require a shared/local pair for every change.
- Substantive shared issues should normally use separate branches.
- After validation and review, before merging the issue branch, move resolved
  issue files to `closed-issues/<id>-<description>.md` and remove their active
  index entries. Record resolution, meaningful validation, and limitations.
  Dates and commit hashes are optional; Git records history. Keep deferred
  issues in the active index with an explicit deferred status.
- Refer to repository issues as `issue 3` and GitHub issues explicitly as
  `GitHub #3` or by URL; their number spaces are independent. Commit messages
  should reference the relevant shared issue when applicable.
- Public GitHub updates still require authorization. Do not copy private
  notes into shared files without reviewing them for contributor relevance.
  Closing a child does not imply that its parent is complete.

# Local issues

- Use ignored `LOCAL_ISSUES.md` for current implementation tasks, following
  `LOCAL_ISSUES.md.example`. Assign IDs L1, L2, etc., and retain its
  last-assigned-ID counter even when completed history is purged.
- Each local issue is a coherent unit of work with a problem, acceptance
  criteria, and status; a shared-issue reference is optional.
- After validation, when ready to commit, move completed entries into ignored
  `CLOSED_LOCAL_ISSUES.md` with a short resolution and validation receipt.
  A final commit hash is not required. Purge these disposable records only
  after work is committed and durable findings are promoted as needed.
- ROADMAP.md is retired. Existing dev-archive/ records remain private legacy
  history; do not rename them or use them as the active tracker. Migrated
  issues retain their legacy ID once for traceability.

# Development notes

- Use untracked `dev-notes/` for diagnostics, experiments, design sketches,
  and intermediate findings; create it when needed.
- Notes are non-authoritative and may become stale. Promote durable findings
  into maintained documentation, tests, source comments, or task contracts.
