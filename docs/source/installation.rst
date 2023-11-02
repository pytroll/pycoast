Installation
============

The below sections describe how to install both the pycoast python library
and additional data files that maybe required to use some features of pycoast.

If you have any trouble with the installation of the package or the data files
described below, please file a bug report on GitHub:

https://github.com/pytroll/pycoast/

Package installation
--------------------

Pycoast can be installed in an existing Python environment via pip or in a
conda environment via ``conda`` using the conda-forge channel. To use pip:

.. code-block:: bash

    pip install pycoast

Alternatively, with conda:

.. code-block:: bash

    conda install -c conda-forge pycoast

Installation of shape files
---------------------------

To use the features of pycoast that draw country or other political borders,
rivers, and lakes, shapefiles from the
`SOEST GSHHG <https://www.soest.hawaii.edu/pwessel/gshhg/>`_ website must be
installed. Download the zipped GSHHS and WDBII shapefiles. At the time of
writing the current zip file can be found at:

https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip

Unzip the files to a data directory (hereafter *DB_DATA_ROOT*).
The absolute path/name of this directory is called *db_root_path*
in the code examples used elsewhere in the documentation.

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

To use the features of Pycoast that depend on city locations or names, one or
more files from `GeoNames <https://www.geonames.org/>`_ must be downloaded
and made available in the same *DB_DATA_ROOT* directory created in the above
GSHHG shapefile download. GeoNames releases multiple lists of city information
available from their file archive:

https://download.geonames.org/export/dump/

There are files that contain city information for cities with a population
larger than 500, 1000, 5000, and 15000. Only one of these files needs to be
downloaded depending on your needs. At the time of writing the URLs for
these files are:

* https://download.geonames.org/export/dump/cities500.zip
* https://download.geonames.org/export/dump/cities1000.zip
* https://download.geonames.org/export/dump/cities5000.zip
* https://download.geonames.org/export/dump/cities15000.zip

Once downloaded, extract the single cities .txt file inside and move it to
a new ``DB_DATA_ROOT/CITIES/`` directory. Currently, Pycoast requires that
the file be named "cities.txt". The structure of *DB_DATA_ROOT* should now be::

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
    └─── CITIES
        └── cities.txt


The PyCoast API documentation explains in detail how to use this city
information via the :meth:`~pycoast.cw_base.ContourWriterBase.add_cities` method.
