#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2016 Adam.Dybbroe

# Author(s):

#   Adam.Dybbroe <a000680@c20671.ad.smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
"""

from PIL import Image
from pycoast import ContourWriterCairo
img = Image.open('/home/a000680/Downloads/ctt_1024.png')
proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
area_extent = (-5570248.4773392612, -5567248.074173444,
               5567248.074173444, 5570248.4773392612)
area_def = (proj4_string, area_extent)
cw = ContourWriterCairo('/home/a000680/data/shapes')
cw.add_coastlines(
    img, (proj4_string, area_extent), resolution='l', width=.5, outline='black')
img.show()
