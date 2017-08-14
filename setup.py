# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011, 2014, 2016  Esben S. Nielsen
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import sys
from setuptools import setup


def get_version():
    if sys.version_info >= (3, 5):
        import importlib
        spec = importlib.util.spec_from_file_location('version', 'pycoast/version.py')
        version = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(version)
    else:
        # python 2.7 doesn't have the `importlib.util` package
        # python 3.4 doesn't have the `module_from_spec` method
        import imp
        version = imp.load_source('pycoast.version', 'pycoast/version.py')
    return version.__version__

requires = ["pyshp", 'numpy', 'pyproj', 'pillow', 'six']

setup(name='pycoast',
      version=get_version(),
      description='Writing of coastlines, borders and rivers to images in Python',
      author='Esben S. Nielsen',
      author_email='esn@dmi.dk',
      packages=['pycoast', 'pycoast.tests'],
      install_requires=requires,
      test_suite='pycoast.tests.test_pycoast.suite',
      zip_safe=False,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Programming Language :: Python',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering'
      ]
      )
