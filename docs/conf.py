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

mathjax3_config = {
    'loader': {'load': ['[tex]/ams']},
    'tex': {
        'macros': {
            'C': r'{\mathbb{C}}',
            'R': r'{\mathbb{R}}',
            'Ci': r'{\mathrm{i}}',
            'bm': [r'\boldsymbol{#1}', 1],
            'boost': [r'{\bf H}{#1}_{\hat{\boldsymbol{m}}}(\beta)', 1],
            'rotat': [r'{\bf U}{#1}_{\hat{\boldsymbol{n}}}(\phi)', 1],
            'pauli': [r'{\boldsymbol{\sigma}}_{#1}', 1],
        }
    }
}

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
