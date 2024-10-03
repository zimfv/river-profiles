# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys

sys.path.insert(0, os.path.abspath("src\bguriver"))

project = 'BGU-River'
copyright = '2024, Fedor Zimin'
author = 'Fedor Zimin'
release = '0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.todo", # Support for todo items
			  "sphinx.ext.viewcode", # Add links to highlighted source code
			  "sphinx.ext.autodoc", # Include documentation from docstrings
			  #"sphinx.ext.napoleon", # Support for NumPy and Google style docstrings
			  #"sphinx.ext.imgmath", # Render math as images
			  #"sphinx.ext.mathjax", # Render math via JavaScript
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'sphinx_rtd_theme'
html_theme = 'furo'

html_static_path = ['_static']
