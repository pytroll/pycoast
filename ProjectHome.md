pycoast is a Python package to add coastlines, borders and rivers to raster images using data from the GSHHS and WDBII datasets

[Documentation](http://pycoast.googlecode.com/git/docs/build/html/index.html)
[Latest download](http://code.google.com/p/pycoast/downloads/detail?name=pycoast-0.5.2.tar.gz)

### News ###
  * **2013-02-19**: v0.5.2 released. Improved handling of grid labeling and addition of label placement control. Switch to use of FFT metric in unit tests to allow for small system dependent rendering differences.
  * **2013-01-24**: v0.5.1 released. Longitude markings now handles dateline crossing as well.
  * **2013-01-23**: v0.5.0 released. Implemented handling of areas crossing the dateline. Improved performance of grid line addition so it is now practical usable for globe projections. Unit tests now self contained.
  * **2012-05-15**: v0.4.0 released. Graticule can now be added to images.
  * **2011-12-05**: v0.3.1 released. Bugfix: Corrected add\_coastlines\_to\_file. Increased unittest coverage to include `*``_`to\_file methods.
  * **2011-12-02**: v0.3.0 released. Bugfixes to improve accuracy. Added offset arguments. Added MANIFEST.in and GPL text to all source files.
  * **2011-10-07**: v0.2.1 released. Generalized unittests to work with environment variable
  * **2011-09-27**: v0.2.0 released. Can now use AGG for high quality images.
  * **2011-09-26**: v0.1.2 released. Now displays sensible for polygons with non-valid projected parts (like the geos projection).
  * **2011-09-23**: First release.

![http://pycoast.googlecode.com/git/docs/source/images/euro_coast_grid_agg.png](http://pycoast.googlecode.com/git/docs/source/images/euro_coast_grid_agg.png)

TODO:
  * Re-factor the contourwriter into smaller code units - e.g. elementary drawing class needed.
  * Setup true unittesting (current testing slow and less portable)
  * Full ESRI shape support (in progress)
  * Cities support (in progress)