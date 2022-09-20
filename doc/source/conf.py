# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
from sphinx.ext import autodoc

sys.path.insert(0, os.path.abspath("../source"))

from create_docs_for_mspypeline_config import main as create_docs


sys.path.insert(0, os.path.abspath("./_ext"))  # for custom extensions
sys.path.insert(0, os.path.abspath("../.."))  # to find the package


# -- Project information -----------------------------------------------------

project = 'MSPypeline'
copyright = '2020, Simon Heming'
author = 'Simon Heming'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',
    'IPython.sphinxext.ipython_console_highlighting',
    'IPython.sphinxext.ipython_directive',
    'sphinx.ext.napoleon',
    'sphinx_rtd_theme',
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates', "_autoapi_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# autodoc
autodoc_typehints = "description"
autoclass_content = "class"


# napoleon options
napoleon_google_docstring = False
napoleon_include_init_with_doc = True

# create documentor for description only
# from https://stackoverflow.com/questions/7825263/including-docstring-in-sphinx-documentation

# create docs for mspypeline conf
create_docs()

class DescriptionOnlyDocumenter(autodoc.MethodDocumenter):
    objtype = "descriptiononly"
    priority = 0

    # do not indent the content
    content_indent = ""

    def add_directive_header(self, sig):
        pass

    def get_doc(self, **kwargs):
        res = super(DescriptionOnlyDocumenter, self).get_doc(**kwargs)
        assert len(res) == 1
        try:
            res = [res[0][:res[0].index("Parameters")]]
        except ValueError:
            pass
        return res


def setup(app):
    app.add_autodocumenter(DescriptionOnlyDocumenter)
