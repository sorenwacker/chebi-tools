from setuptools import setup, find_packages

try:
    import versioneer
except ImportError:
    versioneer = None

with open("README.md", "r") as fh:
    long_description = fh.read()

install_requires = [
    "Levenshtein",
    "pandas",
    "rdkit",
    "wget",
    "obonet",
    "pyvis",
    "tqdm",
    "pyarrow"
]

config = {
    "description": long_description,
    "author": "Soren Wacker",
    "url": "https://github.com/sorenwacker",
    "download_url": "https://github.com/sorenwacker/chebi_tools",
    "author_email": "swacker@ucalgary.ca",
    "version": versioneer.get_version() if versioneer is not None else None,
    "cmdclass": versioneer.get_cmdclass() if versioneer is not None else None,
    "install_requires": install_requires,
    "packages": find_packages(),
    "scripts": [],
    "name": "chebi_tools",
}

setup(**config)
