import os
import unittest

from PIL import Image
import numpy as np

from pycoast import ContourWriter

gshhs_root_dir = os.environ['GSHHS_DATA_ROOT']

class Test(unittest.TestCase):

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
