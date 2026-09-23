# Repository instructions

## Commit messages

Every commit subject must use `<type>: <short imperative summary>`.
Choose the prefix that describes the primary purpose of the change:

- `feat:` — Add functionality or intentionally change supported behavior.
- `refact:` — Restructure or simplify implementation while preserving behavior.
- `docs:` — Update documentation, examples of usage, or explanatory reports.
- `format:` — Apply code formatting or configure formatting conventions and tools.
- `chores:` — Maintain the repository, tooling, dependencies, or validation records.

Use these spellings consistently. Keep the subject concise; use the commit body
for rationale and relevant verification details when needed.

Examples:

```text
feat: Add periodic boundary support
refact: Unify energy evaluation across neighbor builders
docs: Explain caller-supplied neighbor arrays
format: Format Python code with Ruff at 119 columns
chores: Record LAMMPS verification results
```

Use `master` as the main branch.

## Python formatting

Format Python code with Ruff using the repository's `pyproject.toml` settings:
line length 119 and Python 3.10 syntax. The `dev` extra pins the formatter version.
Use the existing `catorch3` environment:

```sh
mamba run -n catorch3 python -m pip install -e '.[dev]'
mamba run -n catorch3 python -m ruff format .
mamba run -n catorch3 python -m ruff format --check .
```

Run the formatting check before committing Python changes. Generated documentation
under `docs/_build` is excluded.
