.. pycoast documentation master file, created by
   sphinx-quickstart on Thu Sep 22 15:38:22 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pycoast
=======
Pycost is a Python package to add coastlines, borders and rivers to raster images using data from the GSHHS and WDBII datasets

Installation
------------
Pycoast depends on pyshp_ and PIL_.

Install pycoast and dependencies.

Download the zipped GSHHS and WDBII Shapefiles from SOEST_.
Unzip the files to a data directory (hereafter called *GSHHS_DATA_ROOT*).
The structure of *GSHHS_DATA_ROOT* should now be::

    .
    ├── GSHHS_shp
    │   ├── c
    │   ├── f
    │   ├── h
    │   ├── i
    │   └── l
    └── WDBII_shp
        ├── c
        ├── f
        ├── h
        ├── i
        └── l

Where each dir on the lowest level contains Shapefiles like *GSHHS_shp/c/GSHHS_c_L1.shp*

Usage
-----
Pycoast can be used to add coastlines, borders and rivers to a raster image if the geographic projection of the image and the image extent in projection coordinates are known

.. image:: images/BMNG_clouds_201109181715_areaT2.png

Pycoast can add contours to either a PIL image object:

    >>> from PIL import Image
    >>> from pycoast import ContourWriter
    >>> img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriter('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, area_def, resolution='l', level=4)
    >>> cw.add_rivers(img, area_def, level=5, outline='blue')
    >>> cw.add_borders(img, area_def, outline=(255, 0, 0))
    >>> img.show()
    
or to an image file:

    >>> from pycoast import ContourWriter
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriter('/home/esn/data/gshhs')
    >>> cw.add_coastlines_to_file('BMNG_clouds_201109181715_areaT2.png', area_def, resolution='l', level=4)
    >>> cw.add_rivers_to_file('BMNG_clouds_201109181715_areaT2.png', area_def, level=5, outline='blue')
    >>> cw.add_borders_to_file('BMNG_clouds_201109181715_areaT2.png', area_def, outline=(255, 0, 0))
    
Where the :attr:`area_extent` is the extent of the image in projection coordinates as (x_ll, y_ll, x_ur, x_ur) measured at pixel edges.

The argument to :attr:`ContourWriter` must be replaced with your *GSHHS_DATA_ROOT*.

.. image:: images/euro_coast.png

The resulting (not so pretty) image shows the effect of the various arguments. The :attr:`resolution` keyword argument controls the resolution of the dataset used. It defaults to 'c' for coarse. Increasing the resolution also increases the processing time. The :attr:`level` keyword argument controls the detail level of the dataset used. It defaults to *1* for the lowest detail level.

Instead of a tuple for :attr:`area_def` a pyresample_ :attr:`AreaDefinition` object can be used.

See method docstrings for information about possible argument values see method docstrings.

Creating an image with coastlines only:

    >>> from PIL import Image
    >>> from pycoast import ContourWriter
    >>> img = Image.new('RGB', (425, 425))
    >>> proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
    >>> area_extent = (-5570248.4773392612, -5567248.074173444, 5567248.074173444, 5570248.4773392612)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriter('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, area_def, resolution='l')
    >>> img.show()    

.. image:: images/geos_coast.png

High quality contours using AGG
-------------------------------
The default plotting mode of pycoast uses PIL for rendering of contours. PIL does not support antialiasing and opacity. The AGG engine can be used for making high quality images using the aggdraw_ module.

First install the aggdraw_ module.

Tip: if the building of aggdraw files with::

    agg_array.h:523: error: cast from ‘agg::int8u*’ to ‘unsigned int’ loses precision
    
Try::

    export CFLAGS="-fpermissive"
    
before building.

Using pycaost with AGG for making antialiased drawing:

    >>> from PIL import Image
    >>> from pycoast import ContourWriterAGG
    >>> img = Image.new('RGB', (425, 425))
    >>> proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
    >>> area_extent = (-5570248.4773392612, -5567248.074173444, 5567248.074173444, 5570248.4773392612)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriterAGG('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, (proj4_string, area_extent), resolution='l', width=0.5)
    >>> img.show()
    
.. image:: images/geos_coast_agg.png

and making the not-so-nice image from the first example nice:

    >>> from PIL import Image
    >>> from pycoast import ContourWriterAGG
    >>> img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriterAGG('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, area_def, resolution='l', level=4)
    >>> cw.add_rivers(img, area_def, level=5, outline='blue', width=0.5, outline_opacity=127)
    >>> cw.add_borders(img, area_def, outline=(255, 0, 0), width=3, outline_opacity=32)
    >>> img.show()

.. image:: images/euro_coast_agg.png

See docstrings of :attr:`ContourWriterAGG` methods for argument descriptions.
    
.. _pyshp: http://code.google.com/p/pyshp/
.. _PIL: http://www.pythonware.com/products/pil/
.. _SOEST: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
.. _pyresample: http://code.google.com/p/pyresample/
.. _aggdraw: http://effbot.org/zone/aggdraw-index.htm
