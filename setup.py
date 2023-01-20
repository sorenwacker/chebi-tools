from setuptools import setup, find_packages

import versioneer

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    "Levenshtein",
    "pandas",
    "rdkit",
]

config = {
    "description": long_description,
    "author": "Soren Wacker",
    "url": "https://github.com/sorenwacker",
    "download_url": "https://github.com/sorenwacker/chebi_tools",
    "author_email": "swacker@ucalgary.ca",
    "version": versioneer.get_version(),
    "cmdclass": versioneer.get_cmdclass(),
    "install_requires": install_requires,
    "packages": find_packages(),
    "scripts": [],
    "name": "chebi_tools",
}

setup(**config)
