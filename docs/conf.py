# Configuration file for the Sphinx documentation builder.

project = "epsic"
copyright = "2026, Willem van Straten"
author = "Willem van Straten"

extensions = [
    "breathe",
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
    "epsic": "_doxygen/xml/epsic",
    "util": "_doxygen/xml/util",
}
breathe_default_project = "epsic"

# -- General configuration ---------------------------------------------------

templates_path = ["_templates"]
exclude_patterns = ["_build", ".venv"]

# -- Options for HTML output -------------------------------------------------

html_theme = "sphinx_rtd_theme"

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

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'

# -- Copy button configuration -----------------------------------------------

copybutton_prompt_text = "$ "
