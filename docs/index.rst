.. BGU-River documentation master file, created by
   sphinx-quickstart on Fri Oct  4 23:19:41 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

BGU-River documentation
=======================

Here is a package, with tools to describe the river profile evolution.

The module ``slope_patches`` contains the realisation of the concept of slope patches, which is described in the `article <https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrf.20031>`_ by Leigh Royden and J. Taylor Perron.

The module ``approximation`` contains the difference schemes approximating the euqation (8) from the mentioned article.


Installation
============

The first you should have `python <https://www.python.org/>`_, `pip <https://pypi.org/project/pip/>`_ and `git <https://git-scm.com/>`_ installed.

Then run in terminal or CMD:

.. code-block:: bash

   git clone --depth=1 --branch main https://github.com/zimfv/river-profiles
   cd river-profiles
   pip install -r requirements.txt
   python setup.py bdist_wheel sdist
   pip install .


Full Documentation
==================

.. toctree::
   :maxdepth: 4

   modules
   

