## Version 1.7.1 (2024/07/01)

### Issues Closed

* [Issue 57](https://github.com/pytroll/pycoast/issues/57) - Test failure on debian sid

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 126](https://github.com/pytroll/pycoast/pull/126) - Fix numpy 2 compatibility

#### Features added

* [PR 124](https://github.com/pytroll/pycoast/pull/124) - Switch to pytest-lazy-fixtures

In this release 2 pull requests were closed.


## Version 1.7.0 (2023/11/30)

### Issues Closed

* [Issue 107](https://github.com/pytroll/pycoast/issues/107) - Incompatibility with Pillow 10 ([PR 108](https://github.com/pytroll/pycoast/pull/108) by [@avalentino](https://github.com/avalentino))
* [Issue 95](https://github.com/pytroll/pycoast/issues/95) - Cached overlays are pale
* [Issue 82](https://github.com/pytroll/pycoast/issues/82) - Test failure with Pillow 9.4 ([PR 84](https://github.com/pytroll/pycoast/pull/84) by [@mraspaud](https://github.com/mraspaud))
* [Issue 49](https://github.com/pytroll/pycoast/issues/49) - add install instructions to docs

In this release 4 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 108](https://github.com/pytroll/pycoast/pull/108) - Fix compatibility with Pillow 10 (Draw.textsize versus Draw.textbbox) ([107](https://github.com/pytroll/pycoast/issues/107))
* [PR 98](https://github.com/pytroll/pycoast/pull/98) - Remove special handling of geographic (longlat) CRSes
* [PR 96](https://github.com/pytroll/pycoast/pull/96) - Fix cached images producing different results without caching

#### Features added

* [PR 105](https://github.com/pytroll/pycoast/pull/105) - [pre-commit.ci] pre-commit autoupdate
* [PR 86](https://github.com/pytroll/pycoast/pull/86) - Fix the pycoast tests further
* [PR 85](https://github.com/pytroll/pycoast/pull/85) - Cleanup tests
* [PR 84](https://github.com/pytroll/pycoast/pull/84) - Fix reference images for pillow 9.4 ([82](https://github.com/pytroll/pycoast/issues/82))
* [PR 83](https://github.com/pytroll/pycoast/pull/83) - Factorize font path computation

#### Documentation changes

* [PR 113](https://github.com/pytroll/pycoast/pull/113) - Update shapefile URL in docs to HTTPS
* [PR 96](https://github.com/pytroll/pycoast/pull/96) - Fix cached images producing different results without caching

In this release 10 pull requests were closed.


## Version 1.6.1 (2022/11/07)

### Pull Requests Merged

#### Bugs fixed

* [PR 70](https://github.com/pytroll/pycoast/pull/70) - Fix NameError when writing grid text in certain cases

In this release 1 pull request was closed.


## Version 1.6.0 (2022/10/17)

### Issues Closed

* [Issue 55](https://github.com/pytroll/pycoast/issues/55) - Additional text overlay options ([PR 56](https://github.com/pytroll/pycoast/pull/56) by [@howff](https://github.com/howff))
* [Issue 23](https://github.com/pytroll/pycoast/issues/23) - Pycoast v0.5.2  issue when drawing coastlines ([PR 59](https://github.com/pytroll/pycoast/pull/59) by [@lobsiger](https://github.com/lobsiger))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 65](https://github.com/pytroll/pycoast/pull/65) - Regenerate shapefile test image for pillow 9.2.0
* [PR 61](https://github.com/pytroll/pycoast/pull/61) - Fix 'add_cities' to use modern GeoNames data
* [PR 60](https://github.com/pytroll/pycoast/pull/60) - Fix major/minor grid line parameters not working in add_overlay_from_dict
* [PR 59](https://github.com/pytroll/pycoast/pull/59) - Fix for horizontal (merc) and vertical (tmerc) scratches in medium si… ([23](https://github.com/pytroll/pycoast/issues/23))

#### Features added

* [PR 63](https://github.com/pytroll/pycoast/pull/63) - Rewrite and refactor all of pycoast to meet minimum modern standards
* [PR 62](https://github.com/pytroll/pycoast/pull/62) - Add ability to overlay user specific shapefiles from dict
* [PR 56](https://github.com/pytroll/pycoast/pull/56) - Add 'coords_ref' option to 'add_points' to specify lon/lat or pixel coordinates ([55](https://github.com/pytroll/pycoast/issues/55))

#### Documentation changes

* [PR 63](https://github.com/pytroll/pycoast/pull/63) - Rewrite and refactor all of pycoast to meet minimum modern standards

In this release 8 pull requests were closed.


## Version 1.5.0 (2021/08/18)

### Issues Closed

* [Issue 48](https://github.com/pytroll/pycoast/issues/48) - Raster images without map projection
* [Issue 46](https://github.com/pytroll/pycoast/issues/46) - add_coastlines call is incorrect in documentation

In this release 2 issues were closed.

### Pull Requests Merged

#### Features added

* [PR 54](https://github.com/pytroll/pycoast/pull/54) - Add area and parameter hashing to overlay caching
* [PR 53](https://github.com/pytroll/pycoast/pull/53) - Switch to GitHub Actions for CI

In this release 2 pull requests were closed.


## Version 1.4.0 (2020/06/08)

### Issues Closed

* [Issue 36](https://github.com/pytroll/pycoast/issues/36) - Add points to an image with list of lat/long pairs ([PR 42](https://github.com/pytroll/pycoast/pull/42))

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 45](https://github.com/pytroll/pycoast/pull/45) - Fix cached overlays always being regenerated
* [PR 44](https://github.com/pytroll/pycoast/pull/44) - Fix adding textbox in the add_points method

#### Features added

* [PR 42](https://github.com/pytroll/pycoast/pull/42) - Add "add_points" method for point/symbol data ([36](https://github.com/pytroll/pycoast/issues/36))

In this release 3 pull requests were closed.


## Version 1.3.2 (2019/12/06)

### Issues Closed

* [Issue 39](https://github.com/pytroll/pycoast/issues/39) - distorted/strange coastlines with pyproj-2.4.2 ([PR 40](https://github.com/pytroll/pycoast/pull/40))

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 40](https://github.com/pytroll/pycoast/pull/40) - Fix compatibility with pyproj 2.4.2 and reduce generated warnings ([39](https://github.com/pytroll/pycoast/issues/39))
* [PR 38](https://github.com/pytroll/pycoast/pull/38) - Remove unnecessary casting in adding overlay from dictionary

In this release 2 pull requests were closed.


## Version 1.3.1 (2019/11/07)

### Pull Requests Merged

#### Bugs fixed

* [PR 37](https://github.com/pytroll/pycoast/pull/37) - Fix the add overlay function to accept `minor_is_tick` as a boolean

In this release 1 pull request was closed.


## Version 1.3.0 (2019/10/25)

### Issues Closed

* [Issue 29](https://github.com/pytroll/pycoast/issues/29) - pycoast compatability issue with pyproj 2+
* [Issue 26](https://github.com/pytroll/pycoast/issues/26) - inconsitency adding lat labelling on the right (on graticules)... ([PR 33](https://github.com/pytroll/pycoast/pull/33))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 33](https://github.com/pytroll/pycoast/pull/33) - Fix #26 ([26](https://github.com/pytroll/pycoast/issues/26))

#### Features added

* [PR 32](https://github.com/pytroll/pycoast/pull/32) - Add dict configuration
* [PR 30](https://github.com/pytroll/pycoast/pull/30) - Convert to RGBA mode when opening the image for adding coastlines or rivers to file.
* [PR 25](https://github.com/pytroll/pycoast/pull/25) - Add coordinate grid overlaying from configuration file

In this release 4 pull requests were closed.


## Version 1.2.3 (2019/06/06)

### Issues Closed

* [Issue 27](https://github.com/pytroll/pycoast/issues/27) - Not compatible with pyproj>=2.2.0 ([PR 28](https://github.com/pytroll/pycoast/pull/28))

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 28](https://github.com/pytroll/pycoast/pull/28) - Fix pyproj 2.2+ is_latlong compatibility ([27](https://github.com/pytroll/pycoast/issues/27))

In this release 1 pull request was closed.


## Version 1.2.2 (2019/01/16)

### Pull Requests Merged

#### Bugs fixed

* [PR 24](https://github.com/pytroll/pycoast/pull/24) - Fix default font loading for 'add_grid'

In this release 1 pull request was closed.


## Version 1.2.1 (2018/11/12)

### Pull Requests Merged

#### Documentation changes

* [PR 22](https://github.com/pytroll/pycoast/pull/22) - Add deprecation warning if using ContourWriter

In this release 1 pull request was closed.


## Version 1.2.0 (2018/11/12)

### Backwards Incompatibility

* `ContourWriter` is now named `ContourWriterPIL`
* `ContourWriterAGG` is now the preferred writer and is aliased to `ContourWriter`

### Issues Closed

* [Issue 17](https://github.com/pytroll/pycoast/issues/17) - Test failures with pyshp 2.0.0
* [Issue 16](https://github.com/pytroll/pycoast/issues/16) - Add/fix documentation and examples about using pyresample AreaDefs
* [Issue 14](https://github.com/pytroll/pycoast/issues/14) - Add documentation for `add_shapes` method
* [Issue 13](https://github.com/pytroll/pycoast/issues/13) - Wrongly positioned coastlines on near sided perspective projection

In this release 4 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 19](https://github.com/pytroll/pycoast/pull/19) - Fix base writer to work with shapely 2.0+ ([18](https://github.com/pytroll/pycoast/issues/18))

#### Documentation changes

* [PR 21](https://github.com/pytroll/pycoast/pull/21) - Add AUTHORS list
* [PR 20](https://github.com/pytroll/pycoast/pull/20) - Cleanup documentation and switch to versioneer for version number

In this release 3 pull requests were closed.


## v1.1.0 (2017-12-08)

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

## v1.0.0 (2017-08-15)

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

## v0.6.1 (2017-05-18)

- Update changelog. [Panu Lahtinen]

- Bump version: 0.6.0 → 0.6.1. [Panu Lahtinen]

- Add missing module name from Proj() call. [Panu Lahtinen]

- Create projection before using it. [Panu Lahtinen]

## v0.6.0 (2017-05-09)

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

## v0.5.4 (2016-02-21)

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

## v0.5.3 (2016-02-21)

- Bugfix: The section is called "coasts", plural... [Martin Raspaud]

- Bugfix: the refactoring used only coastal style. [Martin Raspaud]

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


## v0.5.2 (2013-02-19)

- Built docs. [Esben S. Nielsen]

- Hrobs changes and FFT metric for unit test. [Esben S. Nielsen]

- Flexible grid labeling and placement implemented. [Esben S. Nielsen]

## v0.5.1 (2013-01-24)

- Lon markings now account for dateline too. [Esben S. Nielsen]

## v0.5.0 (2013-01-23)

- Updated doc image. [Esben S. Nielsen]

- Updated docs. [Esben S. Nielsen]

- Test updated. [Esben S. Nielsen]

- Implemented correct dateline handling and updated tests. [Esben S.
  Nielsen]

## v0.4.0 (2012-09-20)

- Added all of docs/build/html. [Esben S. Nielsen]

- Modified comment. [Esben S. Nielsen]

- Added graticule computation from Hrob. [Esben S. Nielsen]

## v0.3.1 (2011-12-05)

- Corrected bug in add_coastlines_to_file. [Esben S. Nielsen]

## v0.3.0 (2011-12-02)

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
