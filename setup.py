"""Setuptools-based packaging and installation configuration."""

from setuptools import find_packages, setup

import versioneer

requires = ["aggdraw", "pyshp", "numpy", "pyproj", "pillow"]

extras_require = {
    "docs": [
        "sphinx",
        "pyresample",
        "pytest",
        "pytest-lazy-fixtures",
        "sphinx_rtd_theme",
        "sphinxcontrib-apidoc",
    ],
    "tests": ["pyresample", "pytest", "pytest-cov", "coverage", "coveralls", "pytest-lazy-fixtures"],
}

with open("README", "r") as readme_file:
    long_description = readme_file.read()

setup(
    name="pycoast",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Writing of coastlines, borders and rivers to images in Python",
    long_description=long_description,
    license="Apache-2.0",
    license_files=["LICENSE.txt"],
    author="Esben S. Nielsen",
    author_email="esn@dmi.dk",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requires,
    extras_require=extras_require,
    python_requires=">3.9",
    zip_safe=False,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Programming Language :: Python",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
    ],
)
