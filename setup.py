from setuptools import find_packages, setup

with open("README.md", "r") as file:
	long_description = file.read()

setup(
	name="bguriver", 
	version="0.1.0", 
	description="Research: River Profile Evolution", 
	package_dir={"": "src"}, 
	packages=find_packages(where="src"), 
	long_description=long_description,
	long_description_content_type="text/markdown", 
	url="https://github.com/zimfv/river-profiles", 
#	author="", 
#	author_email="", 
#	license="", 
#	classifiers=[],
	install_requiers=["numpy", "scipy"], 
)
