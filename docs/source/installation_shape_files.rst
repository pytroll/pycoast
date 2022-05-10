Installation of shape files
---------------------------
Pycoast depends on pyshp_ and PIL_.

Install pycoast and dependencies.

Download the zipped GSHHS and WDBII Shapefiles from SOEST_.
Unzip the files to a data directory (hereafter *DB_DATA_ROOT*).
The absolute path/name of this directory is called *db_root_path*
in the code examples further down.
The structure of *DB_DATA_ROOT* should now be::

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

Where each dir on the lowest level contains Shapefiles like
*GSHHS_shp/c/GSHHS_c_L1.shp, WDBII_shp/WDBII_border_c_L1.shp*

.. _SOEST: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
.. _PIL: https://pillow.readthedocs.io/en/stable/
.. _pyshp: http://code.google.com/p/pyshp/
