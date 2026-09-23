# Documentation source

Pages are Markdown with MyST/Sphinx directives. Run from the repository root:

```sh
mamba run -n catorch3 python -m pip install -e '.[docs]'
mamba run -n catorch3 python -m sphinx -n -W --keep-going -b html docs docs/_build/html
mamba run -n catorch3 python -m sphinx -n -W --keep-going -b doctest docs docs/_build/doctest
```

The local site starts at `docs/_build/html/index.html`.
Generated HTML is intentionally not committed.
