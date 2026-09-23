# Repository instructions

## Commit messages

Every commit subject must use `<type>: <short imperative summary>`.
Choose the prefix that describes the primary purpose of the change:

- `feat:` — Add functionality or intentionally change supported behavior.
- `refact:` — Restructure or simplify implementation while preserving behavior.
- `docs:` — Update documentation, examples of usage, or explanatory reports.
- `chores:` — Maintain the repository, tooling, dependencies, or validation records.

Use these spellings consistently. Keep the subject concise; use the commit body
for rationale and relevant verification details when needed.

Examples:

```text
feat: Add periodic boundary support
refact: Unify energy evaluation across neighbor builders
docs: Explain caller-supplied neighbor arrays
chores: Record LAMMPS verification results
```

Use `master` as the main branch.
