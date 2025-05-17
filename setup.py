"""
DOLFINx Modal Structural Solver - Setup script

Author: Sreekar Reddy, Sajjala
Contact: sreekar2858.tech
License: GPL-3.0
"""
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setup(
    name="dolfinx-modal-solver",
    version="0.1.0",
    author="Sreekar Reddy, Sajjala",
    author_email="sreekar2858.tech",
    description="A modal structural solver using DOLFINx",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/sreekar-reddy/dolfinx",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
)
