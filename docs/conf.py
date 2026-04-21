# Configuration file for the Sphinx documentation builder.

project = "epsic"
copyright = "2026, Willem van Straten"
author = "Willem van Straten"

extensions = [
    "breathe",
    "exhale",
    "myst_parser",
    "sphinx_click",
    "sphinx_copybutton",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

# -- Breathe (Doxygen XML) configuration ------------------------------------

breathe_projects = {
    "epsic": "_build/doxygen_xml",
}
breathe_default_project = "epsic"

# -- Exhale configuration ---------------------------------------------------

exhale_args = {
    "verboseBuild": True,
    "containmentFolder": "./api",
    "rootFileName": "index.rst",
    "doxygenStripFromPath": "..",
}

# -- General configuration ---------------------------------------------------

templates_path = ["_templates"]
exclude_patterns = ["_build", ".venv"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"
# html_static_path = ["_static"]

# -- MyST (Markdown) support -------------------------------------------------

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

myst_enable_extensions = [
    "attrs_inline",
    "colon_fence",
    "amsmath",
    "dollarmath",
]

# -- Copy button configuration -----------------------------------------------

copybutton_prompt_text = "$ "
