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
  tracked by the private roadmap or GitHub issues.
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
- Use test-driven development: establish a failing regression test or other
  reproducible check before changing behavior when practical. Add or update
  meaningful tests for new or modified functionality.
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

# Planning and task tracking

- Keep `ROADMAP.md` as the private authoritative list of active, pending,
  and explicitly deferred tasks. Use short stable IDs such as `T012` or
  `BUG007`; check both ROADMAP.md and dev-archive/ before assigning an ID.
- Give work a problem, scope, acceptance criteria, and relevant dependencies.
  One task should normally be a coherent reviewable unit, not an arbitrary
  implementation step.
- After validation, when resolved work is ready to commit, archive it in
  `dev-archive/<yyyymmdd>-<task-id>-<description>.md` with an explicit status
  (`done`, `deferred`, `dropped`, or `moved-to-github`), resolution, validation,
  limitations, and follow-up references. A commit hash is not required yet.
- Remove resolved tasks from ROADMAP.md after archiving. Deferred work may
  remain on the roadmap when it is still intended future work.
- GitHub issues are public, contributor-facing records; private roadmap,
  notes, and archives are not automatically suitable for publication.
  Update or close authorized GitHub issues only when the governing criteria
  and validation are satisfied. Closing a child does not close its parent.
- Migration to global/local issue files is pending separate review. Do not
  create competing trackers, renumber old tasks, or publish private records
  as part of ordinary issue work.

# Development notes

- Use untracked `dev-notes/` for diagnostics, experiments, design sketches,
  and intermediate findings; create it when needed.
- Notes are non-authoritative and may become stale. Promote durable findings
  into maintained documentation, tests, source comments, or task contracts.
- Keep local backups under untracked `bak/`; preserve existing backups.
