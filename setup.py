# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011-2022 PyCoast Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""Setuptools-based packaging and installation configuration."""

from setuptools import setup

import versioneer

requires = ["aggdraw", "pyshp", "numpy", "pyproj", "pillow"]

extras_require = {
    "docs": [
        "sphinx",
        "pyresample",
        "pytest",
        "pytest-lazy-fixture",
        "sphinx_rtd_theme",
        "sphinxcontrib-apidoc",
    ],
    "tests": ["pyresample", "pytest", "pytest-cov", "coverage", "coveralls", "pytest-lazy-fixture"],
}

with open("README", "r") as readme_file:
    long_description = readme_file.read()

setup(
    name="pycoast",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Writing of coastlines, borders and rivers to images in Python",
    long_description=long_description,
    author="Esben S. Nielsen",
    author_email="esn@dmi.dk",
    packages=["pycoast", "pycoast.tests"],
    include_package_data=True,
    install_requires=requires,
    extras_require=extras_require,
    python_requires=">3.7",
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
)
