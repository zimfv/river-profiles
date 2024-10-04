# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'BGU-River'
copyright = '2024, Fedor Zimin'
author = 'Fedor Zimin'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os, sys

sys.path.insert(0, os.path.abspath('../src/bguriver'))

# Extensions list
extensions = ["sphinx.ext.todo", 
			  "sphinx.ext.viewcode", 
			  "sphinx.ext.autodoc", 
              #"sphinx.ext.autosummary",
			  "sphinx.ext.mathjax",
]

# MathJax 3 configuration
mathjax3_config = {
    'tex': {
        #'tags': 'all',  # Automatically number all equations
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
    }
}

html_js_files = [
    'MathJax/es5/tex-mml-chtml.js',  # Adjust the path as necessary
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'sphinx_rtd_theme'
html_theme = 'pydata_sphinx_theme'
html_static_path = ['_static']
