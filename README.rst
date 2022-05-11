PyCoast
=======

.. image:: https://github.com/pytroll/pycoast/workflows/CI/badge.svg?branch=main
    :target: https://github.com/pytroll/pycoast/actions?query=workflow%3A%22CI%22

.. image:: https://coveralls.io/repos/github/pytroll/pycoast/badge.svg?branch=main
    :target: https://coveralls.io/github/pytroll/pycoast?branch=main

.. image:: https://img.shields.io/pypi/v/pycoast.svg
        :target: https://pypi.python.org/pypi/pycoast

.. image:: https://results.pre-commit.ci/badge/github/pytroll/pycoast/main.svg
   :target: https://results.pre-commit.ci/latest/github/pytroll/pycoast/main
   :alt: pre-commit.ci status

Python package for adding coastlines, borders, rivers, lakes, cities, and other
overlays to raster images.

Installation
------------

PyCoast can be installed from PyPI using pip::

    pip install pycoast

Or with conda using the conda-forge channel::

    conda install -c conda-forge pycoast

Example
-------

::

    >>> from PIL import Image
    >>> from pycoast import ContourWriterAGG
    >>> img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriterAGG('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, area_def, resolution='l', level=4)
    >>> cw.add_rivers(img, area_def, level=5, outline='blue')
    >>> cw.add_borders(img, area_def, outline=(255, 0, 0))
    >>> img.show()
