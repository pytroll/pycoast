# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011-2018 PyCoast Developers
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from setuptools import setup
import versioneer

requires = ['aggdraw', 'pyshp', 'numpy', 'pyproj', 'pillow', 'six']

with open('README', 'r') as readme_file:
    long_description = readme_file.read()

setup(name='pycoast',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description='Writing of coastlines, borders and rivers to images in Python',
      long_description=long_description,
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
