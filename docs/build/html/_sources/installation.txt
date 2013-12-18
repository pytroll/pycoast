Installation
------------
Pycoast depends on pyshp_ and PIL_.

Install pycoast and dependencies.

Download the zipped GSHHS and WDBII Shapefiles from SOEST_.
Unzip the files to a data directory (hereafter called *GSHHS_DATA_ROOT*).
The structure of *GSHHS_DATA_ROOT* should now be::

    .
    ├── GSHHS_shp
    │   ├── c
    │   ├── f
    │   ├── h
    │   ├── i
    │   └── l
    └── WDBII_shp
        ├── c
        ├── f
        ├── h
        ├── i
        └── l

Where each dir on the lowest level contains Shapefiles like *GSHHS_shp/c/GSHHS_c_L1.shp*

.. _SOEST: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
.. _PIL: http://www.pythonware.com/products/pil/
.. _pyshp: http://code.google.com/p/pyshp/
