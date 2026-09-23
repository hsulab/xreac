# Building and editing the documentation

## Build locally

All narrative pages are Markdown parsed by MyST. API sections embed Sphinx
autodoc directives in `eval-rst` fences, so signatures and docstrings come from the installed
package. The docs extra accepts compatible Sphinx/MyST versions so the existing
`catorch3` Markdown-conversion tools can keep their dependencies.

From the repository root:

```sh
mamba run -n catorch3 python -m pip install -e '.[docs]'
mamba run -n catorch3 python -m sphinx -n -W --keep-going -b html docs docs/_build/html
mamba run -n catorch3 python -m sphinx -n -W --keep-going -b doctest docs docs/_build/doctest
```

Open `docs/_build/html/index.html`. The HTML build fails on unresolved references
and other warnings. The doctest builder executes the marked quickstart, small-cell,
and ASE examples. Ordinary Python code blocks are illustrative and are not
executed. LAMMPS, MPI, LaTeX, and Matplotlib are not needed for either build.

Alternatively, run `mamba run -n catorch3 make -C docs html` or replace `html`
with `doctest` or `linkcheck`. The optional link checker needs network access and
may be affected by upstream rate limits. HTML output and doctest logs are ignored
by Git; commit Markdown sources and configuration instead.

The release displayed in the site comes from installed `xreac` package metadata.
Reinstall the editable package after changing `pyproject.toml`'s version.

## Editing pages

Add new `.md` pages to the toctree in `index.md`. Use relative Markdown links for
guide pages, Python-domain roles for API cross-references, and Sphinx download
roles for repository artifacts that should be copied into the built site.
Do not link generated pages back to local filesystem paths or assume a GitHub
repository URL. The [MyST documentation](https://myst-parser.readthedocs.io/)
describes Markdown directives and cross-references.

When changing calculation behavior, update its guide and API documentation and
run the relevant scientific tests. Routine documentation edits need the strict
HTML build and any affected doctests, rather than rerunning LAMMPS benchmarks.
