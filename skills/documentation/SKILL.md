---
name: documentation
description: Create or substantially revise aenet maintained documentation, Fortran or C interface comments, CLI help, and runnable examples using the repository's current documentation sources.
---

# Documentation

Follow AGENTS.md for planning and approval. Inspect the affected interface,
tests, neighboring prose, doc/README.md, and doc/Makefile before editing.

## Choose the maintained source

The current manual source is doc/aenet-doc.org. Inspect its includes and the
files under doc/input-files/ to locate the source of a particular section.
README.org is also exported to README.md by doc/Makefile. Do not overwrite
independent README edits: reconcile source and output within approved scope.
Some input examples may be maintained directly as text; inspect the actual
file and rule rather than assuming every file is generated.

Update maintained source, then regenerate relevant tracked exports when tools
are available. Current exports include Markdown, text, and PDF. Do not treat
Sphinx or Read the Docs as configured here; introducing them is separate work.

Keep procedure comments and C interface documentation close to their code.
Describe purpose, arguments, outputs, units, shapes/layout, allocation
ownership, side effects, and failure behavior where consequential. Do not
impose Python docstring conventions on Fortran.

## Write and validate

Document changed input keywords, CLI behavior, defaults, units, file formats,
and compatibility implications. Distinguish executables from the native
library and companion Python interfaces. Bound scientific and performance
claims to evidence and state assumptions needed to reproduce examples.

Use small deterministic examples through supported public interfaces. Test
changed executable examples with appropriate fixtures; label pseudocode and
intentionally partial examples. Avoid machine-specific paths.

Inspect export rules before running them. `make -C doc md` can update both
manual Markdown and the root README; select the individual target when only
one export is in scope. Text exports use Emacs/Org, Markdown also uses Pandoc,
and PDF export requires the relevant LaTeX toolchain. Inspect tools actually
available; do not claim an export passed when a prerequisite is absent.

For changed rendered documentation, regenerate and inspect affected output
for equations, code, links, tables, and layout. Review the diff for unrelated
export churn. Shared policy/skill Markdown does not require manual export.
Record intentionally untested examples, unavailable rendering checks, and
follow-up work in the handoff.
