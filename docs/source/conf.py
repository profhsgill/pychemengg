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
sys.path.insert(0, os.path.abspath('../..'))


# -- Project information -----------------------------------------------------

project = 'PyChemEngg'
copyright = '2021, Harvinder Singh Gill'
author = 'Harvinder Singh Gill'

# The full version, including alpha/beta/rc tags
__location__ = os.path.realpath(
                    os.path.join(os.getcwd(), os.path.dirname(__file__)))
path_to_version_file =  os.path.dirname(os.path.dirname(__location__))
# path_to_version_file = os.path.dirname(os.path.dirname(os.getcwd()))
version_file = "pychemengg/Version.txt"
quotes = ('"', "'")
with open(os.path.join(path_to_version_file, version_file),"r") as file:
    for line in file:
        if line.startswith('__version__'):
            _, _, current_version = line.split()
            for quote in quotes:
                if current_version.startswith(quote):
                    current_version = current_version.replace(quote,'')            
            break

version = current_version
release = current_version


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinxcontrib.katex",
    "sphinx.ext.autodoc",    
    "sphinx.ext.autosummary",    
    "sphinx.ext.viewcode",
    "nbsphinx",   
    "sphinx.ext.coverage",
    "numpydoc",
    "sphinx.ext.intersphinx"
    ]
    

autosummary_generate = True  # Turn on sphinx.ext.autosummary
numpydoc_show_class_members = False
# as per the katex documentation
# add the following to render math equations using katex
# It will add appropriate javascript file locations to 
# generated html files so that math equations can be parsed
katex_css_path = \
    'https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.css'
katex_js_path = \
    'https://cdn.jsdelivr.net/npm/katex@0.12.0/dist/katex.min.js'
katex_autorender_path = \
    'https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.12.0/contrib/auto-render.min.js'
    
    
# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

# RTD THEME CUSTOMIZATION START
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
html_static_path = ['_static']
def setup(app):
    app.add_css_file('my_theme.css')
html_style = 'css/my_theme.css' 
# RTD END

# html_theme = 'sphinx_rtd_theme'

# OTHER RTD theme options START
html_theme_options = {
    #'analytics_id': 'UA-XXXXXXX-1',  #  Provided by Google in your dashboard
    #'analytics_anonymize_ip': False,
    #'logo_only': False,
    #'display_version': True,
    'prev_next_buttons_location': 'both', #  = bottom, top, or both
    # 'style_external_links': False,
    # 'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'white',
    # Toc options
    # 'collapse_navigation': True,
    # 'sticky_navigation': True,
    # 'navigation_depth': 4,
    # 'includehidden': False,
    # 'titles_only': False
}
# OTHER RTD theme options END


# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
# html_logo = 'artwork/shap_logo_white.png'


# The name of an image file (relative to this directory) to use as a favicon of
# the docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
#
# html_favicon = 'artwork/favicon.ico'
