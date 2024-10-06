BGU-River: Tools for Describing River Profile Evolution
=======================================================

This package provides tools for modeling and analyzing the evolution of river longitudinal profiles.

The :mod:`bguriver.slope_patches` module implements the concept of slope patches, as described in the article `article <https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrf.20031>`_ by Leigh Royden and J. Taylor Perron.

The :mod:`bguriver.approximation` module contains numerical schemes for approximating equation (8) from the mentioned article:

.. math::
   \frac{\partial \lambda}{\partial \tau} + \left(\frac{\partial \lambda}{\partial \chi}\right)^n = \nu(\tau, \chi)

This package is developed by Fedor Zimin for `Liran Goren <https://sites.google.com/site/gorenliran/home>`_'s research group at `Ben Gurion University of the Negev <https://www.bgu.ac.il/en>`_.

The project repository can be found `here <https://github.com/zimfv/river-profiles>`_.

You can read the project documentation `here <https://zimfv.github.io/river-profiles/>`_.


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


Full Documentation Contents
===========================

.. toctree::
   :maxdepth: 4

   guide
   modules
   

