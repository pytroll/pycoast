#pycoast, Writing of coastlines, borders and rivers to images in Python
# 
#Copyright (C) 2011  Esben S. Nielsen
#
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import unittest

from PIL import Image, ImageFont
import numpy as np

from pycoast import ContourWriter

def tmp(f):
    f.tmp = True
    return f

gshhs_root_dir = os.environ['GSHHS_DATA_ROOT']
test_file = 'test_image.png'
grid_file = 'test_grid.png'

class TestPycoast(unittest.TestCase):

    def setUp(self):
        img = Image.new('RGB', (640, 480))
        img.save(test_file)
        img.save(grid_file)
        
    def tearDown(self):    
        os.remove(test_file)
        os.remove(grid_file)
        
class TestPIL(TestPycoast):

    def test_europe(self):
        euro_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_europe.png'))
        euro_data = np.array(euro_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriter(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_rivers(img, area_def, level=5, outline='blue')
        cw.add_borders(img, area_def, outline=(255, 0, 0))
        res = np.array(img)
        self.failUnless(np.array_equal(euro_data, res), 'Writing of contours failed')
        
    def test_europe_file(self):
        euro_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_europe.png'))
        euro_data = np.array(euro_img)

        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriter(gshhs_root_dir)
        cw.add_coastlines_to_file(test_file, area_def, resolution='l', level=4)
        cw.add_rivers_to_file(test_file, area_def, level=5, outline='blue')
        cw.add_borders_to_file(test_file, area_def, outline=(255, 0, 0))
        
        img = Image.open(test_file)
        res = np.array(img)
        self.failUnless(np.array_equal(euro_data, res), 'Writing of contours failed')

    def test_geos(self):
        geos_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_geos.png'))
        geos_data = np.array(geos_img)
        
        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444, 5567248.074173444, 5570248.4773392612)
        area_def = (proj4_string, area_extent)
        cw = ContourWriter(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l')
        res = np.array(img)
        self.failUnless(np.array_equal(geos_data, res), 'Writing of geos contours failed')
        
    def test_grid(self):
        # Only test if functions execute
        img = Image.new('RGB', (640, 480))
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriter(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_grid(img, area_def, (10.0,10.0),(2.0,2.0), fill='blue',
                    outline='blue', minor_outline='blue')
                    
    def test_grid_file(self):
        # Only test if functions execute
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriter(gshhs_root_dir)

        cw.add_coastlines_to_file(grid_file, area_def, resolution='l', level=4)
        cw.add_grid_to_file(grid_file, area_def, (10.0,10.0),(2.0,2.0), fill='blue',
                    outline='blue', minor_outline='blue')

class TestPILAGG(TestPycoast):        

    def test_europe_agg(self):
        from pycoast import ContourWriterAGG
        euro_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_europe_agg.png'))
        euro_data = np.array(euro_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_rivers(img, area_def, level=5, outline='blue', width=0.5, outline_opacity=127)
        cw.add_borders(img, area_def, outline=(255, 0, 0), width=3, outline_opacity=32)
        res = np.array(img)
        self.failUnless(np.array_equal(euro_data, res), 'Writing of contours failed for AGG')
        
    def test_europe_agg_file(self):
        from pycoast import ContourWriterAGG
        euro_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_europe_agg.png'))
        euro_data = np.array(euro_img)

        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines_to_file(test_file, area_def, resolution='l', level=4)
        cw.add_rivers_to_file(test_file, area_def, level=5, outline='blue', 
                              width=0.5, outline_opacity=127)
        cw.add_borders_to_file(test_file, area_def, outline=(255, 0, 0), width=3, 
                               outline_opacity=32)
                               
        img = Image.open(test_file)
        res = np.array(img)
        self.failUnless(np.array_equal(euro_data, res), 'Writing of contours failed for AGG')


    def test_geos_agg(self):
        from pycoast import ContourWriterAGG
        geos_img = Image.open(os.path.join(os.path.dirname(__file__), 
                              'contours_geos_agg.png'))
        geos_data = np.array(geos_img)
        
        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444, 5567248.074173444, 5570248.4773392612)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, (proj4_string, area_extent), resolution='l', width=0.5)
        res = np.array(img)
        self.failUnless(np.array_equal(geos_data, res), 'Writing of geos contours failed for AGG')
    
    @tmp    
    def test_grid_agg(self):
        # Only test if functions execute
        from pycoast import ContourWriterAGG
        img = Image.new('RGB', (640, 480))
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)
        
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_grid(img, area_def, (10.0,10.0),(2.0,2.0), write_text=False,
                    outline='blue',outline_opacity=255,width=1.0,
                    minor_outline='lightblue',minor_outline_opacity=255,minor_width=0.5,
                    minor_is_tick=False)
    
    @tmp                
    def test_grid_agg_file(self):
        # Only test if functions execute
        from pycoast import ContourWriterAGG
        proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines_to_file(grid_file, area_def, resolution='l', level=4)
        cw.add_grid_to_file(grid_file, area_def, (10.0,10.0),(2.0,2.0), write_text=False,
                    outline='blue',outline_opacity=255,width=1.0,
                    minor_outline='lightblue',minor_outline_opacity=255,minor_width=0.5,
                    minor_is_tick=False)
    
