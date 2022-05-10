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


Installation of city names
--------------------------

Add a subdirectory CITIES to your *DB_DATA_ROOT*.
Download cities2022.zip_ from github/pytroll/pycoast.
If this deep link does no work (anymore) you can try
do download it from github.com/pytroll/pycoast/ from
directory tests/test_data/ by pressing the *Download*
button. If all fails ask for help on the google
pytroll_ list or on MSG-1_. Put cities2022.zip into
CITIES. The structure of *DB_DATA_ROOT* should now be::

    .
    ├── GSHHS_shp
    │   ├── c
    │   ├── f
    │   ├── h
    │   ├── i
    │   └── l
    ├── WDBII_shp
    │   ├── c
    │   ├── f
    │   ├── h
    │   ├── i
    │   └── l
    └── CITIES

Where CITIES contains cities2022.zip. Enter this archive
and read file README_PyCoast.txt inside. Follow its instructions.
You will end up with an additional tab delimited text file *cities.txt* in CITIES.
The PyCoast API documentation explains in detail how to call the function *add_cities()*.

.. _SOEST: http://www.soest.hawaii.edu/pwessel/gshhs/index.html
.. _PIL: https://pillow.readthedocs.io/en/stable/
.. _pyshp: http://code.google.com/p/pyshp/
.. _cities2022.zip: https://raw.githubusercontent.com/lobsiger/pycoast/fix_cities/pycoast/tests/test_data/gshhs/CITIES/cities2022.zip
.. _pytroll: https://groups.google.com/g/pytroll/
.. _MSG-1: https://groups.io/g/MSG-1/