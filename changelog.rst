Changelog
=========

v1.1.0 (2017-12-08)
-------------------

- Update changelog. [Panu Lahtinen]

- Bump version: 1.0.0 → 1.1.0. [Panu Lahtinen]

- Merge pull request #11 from TAlonglong/develop. [Panu Lahtinen]

  make _iterate_db loop over a list of levels, L1..LX. GSHHS has 6 leve…

- Make _iterate_db loop over a list of levels, L1..LX. GSHHS has 6
  levels, Borders 3 levels, rivers 11 levels. [Trygve Aspenes]

- Merge pull request #10 from pytroll/feature-refactor. [David Hoese]

  Refactor shape file feature creation

- Add contributor map script. [davidh-ssec]

- Make GSHHS data root optional. [davidh-ssec]

- Remove special exception class. [davidh-ssec]

  It doesn't provide anything that a ValueError doesn't


- Reorganize methods for adding shapes to use the same code. [davidh-
  ssec]

- Change adding features to using a generator. [davidh-ssec]

v1.0.0 (2017-08-15)
-------------------

- Update changelog. [Martin Raspaud]

- Bump version: 0.6.1 → 1.0.0. [Martin Raspaud]

- Merge remote-tracking branch 'origin/master' into pre-master. [Martin
  Raspaud]

- Merge pull request #9 from pytroll/feature-python3. [Martin Raspaud]

  Feature python3

- Fix version import on python 3.4 (again) [davidh-ssec]

- Fix version import on python 3.4. [davidh-ssec]

- Remove unused import in setup.py. [davidh-ssec]

- Change version import to use importlib. [davidh-ssec]

- Use davidh-ssec/aggdraw for travis builds. [davidh-ssec]

- Fix version import in setup.py. [davidh-ssec]

- Fix version import in setup.py and add version to main package init.
  [davidh-ssec]

- Remove aggdraw from python 3 travis executions for now. [davidh-ssec]

- Add python 3 to travis tests. [davidh-ssec]

- Fix python 3 compatiblity in cw_X.py modules. [davidh-ssec]

- Fix tests and MANIFEST to be more modern. [davidh-ssec]

v0.6.1 (2017-05-18)
-------------------

- Update changelog. [Panu Lahtinen]

- Bump version: 0.6.0 → 0.6.1. [Panu Lahtinen]

- Add missing module name from Proj() call. [Panu Lahtinen]

- Create projection before using it. [Panu Lahtinen]

v0.6.0 (2017-05-09)
-------------------

- Update changelog. [Panu Lahtinen]

- Bump version: 0.5.5 → 0.6.0. [Panu Lahtinen]

- Merge branch 'master' into pre-master. [Panu Lahtinen]

- Add bumpversion. [Panu Lahtinen]

- Add bumpversion, gitchangelog and bring versio up-to-date. [Panu
  Lahtinen]

- Update changelog. [Panu Lahtinen]

- Merge branch 'master' into pre-master. [Panu Lahtinen]

- Merge pull request #5 from avalentino/noagg. [Martin Raspaud]

  Skip test on ContourWriterAGG if drawadd is not available

- Skip test on ContourWriterAGG if drawadd is not available. [Antonio
  Valentino]

- Merge pull request #8 from pytroll/feature_extent_degrees. [Panu
  Lahtinen]

  Feature extent degrees

- Drop Python 2.6 support. [Panu Lahtinen]

- Remove texts from the images. [Panu Lahtinen]

- Do more PEP8. [Panu Lahtinen]

- Make geocentric extents the default. [Panu Lahtinen]

- Fix line lengths (PEP8) [Panu Lahtinen]

- Use functions from cw_base instead local copies. [Panu Lahtinen]

- Make it possible to add coastlines to area defs with extents in
  degrees. [Panu Lahtinen]

- Minor PEP8. [Panu Lahtinen]

- Pep8 pretify. [Adam.Dybbroe]

- Pep8. [Adam.Dybbroe]

- Update fonts in two test images. [Martin Raspaud]

v0.5.4 (2016-02-21)
-------------------

- Update changelog. [Martin Raspaud]

- Bump version: 0.5.3 → 0.5.4. [Martin Raspaud]

- Merge branch 'pre-master' [Martin Raspaud]

- Merge branch 'pre-master' of github.com:pytroll/pycoast into pre-
  master. [Martin Raspaud]

- Add missing import of ConfigParser. [Panu Lahtinen]

- Remove skipping shapefiles crossing dateline, fixes global overlays.
  [Panu Lahtinen]

- Update test reference images with ones made in Ubuntu so automatic
  testing might work in Travis. [Panu Lahtinen]

v0.5.3 (2016-02-21)
-------------------

Fix
~~~

- Bugfix: The section is called "coasts", plural... [Martin Raspaud]

- Bugfix: the refactoring used only coastal style. [Martin Raspaud]

Other
~~~~~

- Update changelog. [Martin Raspaud]

- Bump version: 0.5.2 → 0.5.3. [Martin Raspaud]

- Merge branch 'master' into pre-master. [Martin Raspaud]

  Conflicts:
  	.gitignore
  	setup.py

- Add a gitignore file. [Martin Raspaud]

- Merge pull request #2 from mitkin/master. [Martin Raspaud]

  Update setup.py

- Update setup.py. [Mikhail Itkin]

- Fix pyproj missing dependency. [Martin Raspaud]

- Move the test directory to pycoast.tests (for consistency) [Martin
  Raspaud]

- Fix .travis.yml to install aggdraw first. [Martin Raspaud]

- Merge branch 'restructure' into pre-master. [Martin Raspaud]

- Split different readers out of __init__.py, adjust __init__.py so that
  everything still works in the same way as previously. [Panu Lahtinen]

- More generic ignores. [Panu Lahtinen]

- Add a test suite and .travis file for ci. [Martin Raspaud]

- Update the reference images to make the test pass with pillow/aggdraw.
  [Martin Raspaud]

- Add numpy to the list of dependencies. [Martin Raspaud]

- Added documentation about configuration files for pycoast. [Martin
  Raspaud]

- Add setup.cfg for easy rpm generation. [Martin Raspaud]

- Some refactoring and pep8. [Martin Raspaud]

- _add_shapefile_shape bug fixed. [s.cerino]

- Merge branch 'master' into pre-master. [Martin Raspaud]

- Add configuration file reading feature. [Martin Raspaud]

- Merge branch 'master' of https://code.google.com/p/pycoast. [Martin
  Raspaud]

- Fixed (sometimes fatal) ImageDraw import. [Hrobjartur Thorsteinsson]

  ImageDraw and other PIL modules should be imported
  directly from prevailing PIL package.


- Merge branch 'master' of https://code.google.com/p/pycoast. [Martin
  Raspaud]

- Merge branch 'master' of https://code.google.com/p/pycoast. [Martin
  Raspaud]

- Merge branch 'master' of https://code.google.com/p/pycoast. [Martin
  Raspaud]

  Conflicts:
  	pycoast/__init__.py


- Removing the rounding of the pixel indices. (Works with AGG and
  without). [Martin Raspaud]

- Docbuilds. [Hrobjartur Thorsteinsson]

  docbuilds


- Added documentation for polygons and shapefile methods. [Hrobjartur
  Thorsteinsson]

  Added documentation for polygons and shapefile methods.


- Add_polygon and add_shapefile_shape(s) integration testing.
  [Hrobjartur Thorsteinsson]

  add_polygon and add_shapefile_shape(s) integration testing.
  Also included preliminary test data.


- Work in progress setting up shape and cities support. [Hrobjartur
  Thorsteinsson]

  Work in progress setting up shape and cities support


- Removed print line from add_shape routine. [Hrobjartur Thorsteinsson]

  removed print line from add_shape routine


- Make pillow a dependency if PIL is not already there. [Martin Raspaud]

- Fixed fata ImageDraw import. [Hrobjartur Thorsteinsson]

  Fixed importing conflict, affecting some users
  seemingly with mixed installations of PIL/Pillow.

  all PIL imports should be from same package.
  made "from PIL import ImageDraw"


- Adding appertizer image at the front. [Adam Dybbroe]

- Rearranging documentation, and minor editorial stuff. [Adam Dybbroe]

- Bug fix: add_line / add_polygon. [Hrobjartur Thorsteinsson]

  Minor bug fix: add_line / add_polygon exception.


- Added custom shapefile and shape draw routines. [Hrobjartur
  Thorsteinsson]

  custom shapefile and shape draw routines.

  add_shapefile_shape(...)
  add_shapefile_shapes(...)
  add_line(...)
  add_polygon(...)


v0.5.2 (2013-02-19)
-------------------

- Built docs. [Esben S. Nielsen]

- Hrobs changes and FFT metric for unit test. [Esben S. Nielsen]

- Flexible grid labeling and placement implemented. [Esben S. Nielsen]

v0.5.1 (2013-01-24)
-------------------

- Lon markings now account for dateline too. [Esben S. Nielsen]

v0.5.0 (2013-01-23)
-------------------

- Updated doc image. [Esben S. Nielsen]

- Updated docs. [Esben S. Nielsen]

- Test updated. [Esben S. Nielsen]

- Implemented correct dateline handling and updated tests. [Esben S.
  Nielsen]

v0.4.0 (2012-09-20)
-------------------

- Added all of docs/build/html. [Esben S. Nielsen]

- Modified comment. [Esben S. Nielsen]

- Added graticule computation from Hrob. [Esben S. Nielsen]

v0.3.1 (2011-12-05)
-------------------

- Corrected bug in add_coastlines_to_file. [Esben S. Nielsen]

v0.3.0 (2011-12-02)
-------------------

- Bugfixing to improve accuracy. [Esben S. Nielsen]

- Added testing. [Esben S. Nielsen]

- Corrected docs. [Esben S. Nielsen]

- Corrected git doc mess. [Esben S. Nielsen]

- Updated docs. [Esben S. Nielsen]

- Added possiblility to use AGG. Changed API slightly. [Esben S.
  Nielsen]

- Docs messed up by git. Trying to clean. [Esben S. Nielsen]

- Added missing build doc files. [Esben S. Nielsen]

- Corrected invalid reprojection issue for projections like geos. [Esben
  S. Nielsen]

- Rebuild docs. [Esben S. Nielsen]

- Bumped up version. [Esben S. Nielsen]

- Corrected south pole filtering bug. [Esben S. Nielsen]

- Changed link to SOEST. [Esben S. Nielsen]

- Documented project. [Esben S. Nielsen]

- Added license and docs. [Esben S. Nielsen]

- Now handles poles. [Esben S. Nielsen]

- Added docstrings. [Esben S. Nielsen]

- Added test. [Esben S. Nielsen]

- Created package. [Esben S. Nielsen]

- Restructured pixel index calculation. [Esben S. Nielsen]

- Added borders and rivers. [Esben S. Nielsen]

- First version. [Esben S. Nielsen]

- First version. [Esben S. Nielsen]


