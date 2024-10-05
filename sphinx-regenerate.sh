# reinstall the package
pip uninstall bguriver
python setup.py bdist_wheel sdist
pip install .

# regenerate the docs
sphinx-build -M clean docs docs/html
sphinx-apidoc -o docs src/bguriver
sphinx-build -M html docs docs
