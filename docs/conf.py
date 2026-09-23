"""Sphinx configuration shared by local builds and Read the Docs."""

from importlib.metadata import version as package_version
import os

project = "xreac"
author = "xreac contributors"
copyright = "2026, xreac contributors"
release = package_version("xreac")
version = release

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.doctest",
    "sphinx.ext.mathjax",
]
# Autodoc emits reStructuredText internally; keep its parser enabled even
# though every hand-written page is Markdown.
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = ["_build", "README.md"]
language = "en"
nitpicky = True
# External annotation types are rendered as text without fetching intersphinx
# inventories. All project/API cross-references remain strict and work offline.
nitpick_ignore = [("py:class", "numpy.ndarray"), ("py:class", "pathlib.Path")]
autodoc_member_order = "bysource"
autodoc_typehints = "none"
myst_enable_extensions = ["colon_fence", "dollarmath"]
myst_heading_anchors = 3

html_theme = "sphinx_rtd_theme"
html_title = f"xreac {release}"
html_theme_options = {"navigation_depth": 3, "collapse_navigation": False}
html_baseurl = os.environ.get("READTHEDOCS_CANONICAL_URL", "")
html_show_sourcelink = True

# Executed only by the doctest builder; HTML builds do not run calculations.
doctest_global_setup = """
import numpy as np
from xreac import Calculator, ForceField
"""
