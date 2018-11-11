PyCoast
=======

.. image:: https://travis-ci.org/pytroll/pycoast.svg?branch=master
    :target: https://travis-ci.org/pytroll/pycoast

.. image:: https://img.shields.io/pypi/v/pycoast.svg
        :target: https://pypi.python.org/pypi/pycoast

Python package for adding coastlines and borders to raster images.

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
