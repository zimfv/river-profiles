# BGU-River: Tools for Describing River Profile Evolution

This package provides tools for modeling and analyzing the evolution of river longitudinal profiles.

The `bguriver.slope_patches` module implements the concept of slope patches, as described in the article [article](https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrf.20031) by Leigh Royden and J. Taylor Perron.

The `bguriver.approximation` module contains numerical schemes for approximating equation (8) from the mentioned article:

![Image](pics/eq8.png)


This package is developed by Fedor Zimin for [Liran Goren](https://sites.google.com/site/gorenliran/home)'s research group at [Ben Gurion University of the Negev](https://www.bgu.ac.il/en).

The project repository can be found [here](https://github.com/zimfv/river-profiles).

You can read the project documentation [here](https://zimfv.github.io/river-profiles).


# Installation 
The first you should have [python](https://www.python.org/), [pip](https://pypi.org/project/pip/) and [git](https://git-scm.com/) installed.

Then run in terminal or CMD:
```bash
git clone --depth=1 --branch main https://github.com/zimfv/river-profiles
cd river-profiles
pip install -r requirements.txt
python setup.py bdist_wheel sdist
pip install .
```