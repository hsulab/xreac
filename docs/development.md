# Building and publishing the documentation

## Build locally

All narrative pages are Markdown parsed by MyST. API sections embed Sphinx
autodoc directives in `eval-rst` fences, so signatures and docstrings come from the installed
package. The docs extra accepts compatible Sphinx/MyST versions so the existing
`catorch3` Markdown-conversion tools can keep their dependencies. RTD's clean
environment uses the pinned stack in `docs/requirements.txt`.

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

## Read the Docs deployment

The root `.readthedocs.yaml` uses the
[Read the Docs v2 configuration](https://docs.readthedocs.com/platform/stable/config-file/v2.html):

```{literalinclude} ../.readthedocs.yaml
:language: yaml
```

Read the Docs installs the pinned docs stack and the checkout with its `docs`
extra, then builds HTML using
`docs/conf.py`, failing on warnings. It uses Ubuntu 24.04 and Python 3.12. The
documentation imports the real package, including the ASE adapter, but does not
launch LAMMPS or execute the marked calculations during the HTML build.
When provided by RTD, `READTHEDOCS_CANONICAL_URL` sets the canonical site URL.
No project slug or hosted URL is hard-coded.

To publish:

1. Push the repository, including `.readthedocs.yaml`, `docs/`, `pyproject.toml`,
   `src/`, `data/`, and the linked `validation/` artifacts, to the chosen Git host.
2. Import that repository into Read the Docs and authorize the Git integration.
3. Set the project's default branch to **`master`** and build its `latest` version.
4. Inspect the build log, then enable additional release versions as needed.

The project has no Git remote configured in the checkout used to prepare this
documentation. No hosted RTD project or live URL is assumed. Configuration and
local builds prepare the repository for deployment; importing the project and
triggering the first remote build happen in your RTD account.

The provider's [project-import guide](https://docs.readthedocs.com/platform/stable/intro/add-project.html)
covers account and integration setup. No credentials belong in this repository.

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
