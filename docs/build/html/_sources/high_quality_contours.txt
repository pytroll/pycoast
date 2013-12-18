High quality contours using AGG
-------------------------------

The default plotting mode of pycoast uses PIL_ for rendering of contours. PIL
does not support antialiasing and opacity. The AGG engine can be used for
making high quality images using the aggdraw_ module.

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

.. _PIL: http://www.pythonware.com/products/pil/
.. _aggdraw: http://effbot.org/zone/aggdraw-index.htm
