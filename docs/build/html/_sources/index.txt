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

Tip: if the building of aggdraw fails with:

.. code-block:: bash

    agg_array.h:523: error: cast from ‘agg::int8u*’ to ‘unsigned int’ loses precision
    
Try:

.. code-block:: bash

    export CFLAGS="-fpermissive"
    
before building.

Using pycoast with AGG for making antialiased drawing:

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

Adding graticule to images
--------------------------
Pycoast can be used to add graticule to images. For PIL:

    >>> from PIL import Image, ImageFont
    >>> from pycoast import ContourWriter
    >>> import aggdraw
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriter('/home/esn/data/gshhs')
    >>> font = ImageFont.truetype("/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif.ttf",16)
    >>> img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    >>> cw.add_coastlines(img, area_def, resolution='l', level=4)
    >>> cw.add_grid(img, area_def, (10.0,10.0),(2.0,2.0), font,fill='blue',
    ...             outline='blue', minor_outline='blue')
    >>> img.show()

.. image:: images/euro_grid.png
The font argument is optional for PIL if it is not given a default font will be used.

and for AGG:

    >>> from PIL import Image, ImageFont
    >>> from pycoast import ContourWriterAGG
    >>> import aggdraw
    >>> proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    >>> area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    >>> area_def = (proj4_string, area_extent)
    >>> cw = ContourWriterAGG('/home/esn/data/gshhs')
    >>> font = aggdraw.Font("blue", "/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif.ttf", size=16)
    >>> img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    >>> cw.add_coastlines(img, area_def, resolution='l', level=4)
    >>> cw.add_grid(img, area_def, (10.0,10.0),(2.0,2.0),font,
    ...             outline='blue',outline_opacity=255,width=1.0,
    ...             minor_outline='lightblue',minor_outline_opacity=255,minor_width=0.5,
    ...             minor_is_tick=False)
    >>> img.show()

.. image:: images/euro_grid_agg.png

Note the difference in the optional font argument for PIL and AGG. With AGG the font argument is mandatory unless the keyword argument :attr:`write_text=False` is used.

Tip: If the adding graticule with AGG fails with something like:

.. code-block:: bash

    Traceback (most recent call last):
        File "grid_demo_AGG.py", line 13, in <module>
            font=aggdraw.Font("blue", "/usr/share/fonts/truetype/ttf-dejavu/DejaVuSerif.ttf", size=16)
    IOError: cannot load font (no text renderer)

make sure the FREETYPE_ROOT in setup.py of aggdraw points to the correct location e.g. set *FREETYPE_ROOT = "/usr"*

Testing
-------
In order to run the unittests the environment variable GSHHS_DATA_ROOT has to be defined pointing to your *GSHHS_DATA_ROOT*. E.g.:

.. code-block:: bash
    
    $ export GSHHS_DATA_ROOT=/home/esn/data/gshhs

or whatever matches your system and data location.

Subsequently the tests can be run using nosetest:

.. code-block:: bash

    $ cd <pycoast_dir>
    $ nosetests tests/  

    
.. _pyshp: http://code.google.com/p/pyshp/
.. _PIL: http://www.pythonware.com/products/pil/
.. _SOEST: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
.. _pyresample: http://code.google.com/p/pyresample/
.. _aggdraw: http://effbot.org/zone/aggdraw-index.htm
