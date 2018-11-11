#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

import os
import unittest

import numpy as np
from PIL import Image, ImageFont


def tmp(f):
    f.tmp = True
    return f


def fft_proj_rms(a1, a2):
    """Compute the RMS of differences between FFT vectors of a1
    and projection of FFT vectors of a2.
    This metric is sensitive to large scale changes and image noise but
    insensitive to small rendering differences.
    """

    ms = 0

    for i in range(3):
        fr1 = np.fft.fft2(a1[:, :, i])
        fr2 = np.fft.fft2(a2[:, :, i])

        ps1 = np.log10(fr1 * fr1.conj()).real
        ps2 = np.log10(fr2 * fr2.conj()).real

        p1 = np.arctan2(fr1.imag, fr1.real)
        p2 = np.arctan2(fr2.imag, fr2.real)

        theta = p2 - p1
        l = ps2 * np.cos(theta)
        ms += ((l - ps1) ** 2).sum() / float(ps1.size)

    rms = np.sqrt(ms)

    return rms


def fft_metric(data1, data2, max_value=0.1):
    """Execute FFT metric
    """

    rms = fft_proj_rms(data1, data2)
    return rms <= max_value


gshhs_root_dir = os.path.join(os.path.dirname(__file__), 'test_data', 'gshhs')
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
        from pycoast import ContourWriterPIL
        euro_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_europe.png'))
        euro_data = np.array(euro_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterPIL(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_rivers(img, area_def, level=5, outline='blue')
        cw.add_borders(img, area_def, outline=(255, 0, 0))

        res = np.array(img)
        self.assertTrue(fft_metric(euro_data, res),
                        'Writing of contours failed')

    def test_europe_file(self):
        from pycoast import ContourWriterPIL
        euro_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_europe.png'))
        euro_data = np.array(euro_img)

        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterPIL(gshhs_root_dir)
        cw.add_coastlines_to_file(test_file, area_def, resolution='l', level=4)
        cw.add_rivers_to_file(test_file, area_def, level=5, outline='blue')
        cw.add_borders_to_file(test_file, area_def, outline=(255, 0, 0))

        img = Image.open(test_file)
        res = np.array(img)
        self.assertTrue(
            fft_metric(euro_data, res), 'Writing of contours failed')

    def test_geos(self):
        from pycoast import ContourWriterPIL
        geos_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_geos.png'))
        geos_data = np.array(geos_img)

        img = Image.new('RGB', (425, 425))
        proj4_string = \
            '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444,
                       5567248.074173444, 5570248.4773392612)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterPIL(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l')

        res = np.array(img)
        self.assertTrue(
            fft_metric(geos_data, res), 'Writing of geos contours failed')

    def test_grid(self):
        from pycoast import ContourWriterPIL
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_europe.png'))
        grid_data = np.array(grid_img)
        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),
                                               'test_data', 'DejaVuSerif.ttf'),
                                  16)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    font=font, fill='blue', write_text=False,
                    outline='blue', minor_outline='blue')

        res = np.array(img)
        self.assertTrue(fft_metric(grid_data, res), 'Writing of grid failed')

    def test_grid_geos(self):
        from pycoast import ContourWriterPIL
        geos_img = Image.open(
            os.path.join(os.path.dirname(__file__), 'grid_geos.png'))
        geos_data = np.array(geos_img)
        img = Image.new('RGB', (425, 425))
        proj4_string = \
            '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444,
                       5567248.074173444, 5570248.4773392612)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterPIL(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l')
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0), fill='blue',
                    outline='blue', minor_outline='blue',
                    write_text=False)

        res = np.array(img)
        self.assertTrue(
            fft_metric(geos_data, res), 'Writing of geos contours failed')

    def test_grid_file(self):
        from pycoast import ContourWriterPIL
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_europe.png'))
        grid_data = np.array(grid_img)
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines_to_file(grid_file, area_def, resolution='l', level=4)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),
                                               'test_data', 'DejaVuSerif.ttf'),
                                  16)
        cw.add_grid_to_file(grid_file, area_def, (10.0, 10.0), (2.0, 2.0),
                            font=font, fill='blue', write_text=False,
                            outline='blue', minor_outline='blue')

        img = Image.open(grid_file)
        res = np.array(img)
        self.assertTrue(fft_metric(grid_data, res), 'Writing of grid failed')

    def test_dateline_cross(self):
        from pycoast import ContourWriterPIL
        dl_img = Image.open(os.path.join(os.path.dirname(__file__),
                                         'dateline_cross.png'))
        dl_data = np.array(dl_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = '+proj=stere +lon_0=-170.00 +lat_0=60.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),
                                               'test_data',
                                               'DejaVuSerif.ttf'), 16)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    font=font, fill='blue', write_text=False,
                    outline='blue', minor_outline='blue',
                    lon_placement='b', lat_placement='lr')

        res = np.array(img)
        self.assertTrue(fft_metric(dl_data, res),
                        'Writing of dateline crossing data failed')

    def test_dateline_boundary_cross(self):
        from pycoast import ContourWriterPIL
        dl_img = Image.open(os.path.join(os.path.dirname(__file__),
                                         'dateline_boundary_cross.png'))
        dl_data = np.array(dl_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=140.00 +lat_0=60.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),
                                               'test_data',
                                               'DejaVuSerif.ttf'), 16)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    font=font, fill='blue',
                    outline='blue', minor_outline='blue', write_text=False,
                    lon_placement='b', lat_placement='lr')

        res = np.array(img)
        self.assertTrue(fft_metric(dl_data, res),
                        'Writing of dateline boundary crossing data failed')

    def test_grid_nh(self):
        from pycoast import ContourWriterPIL
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_nh.png'))
        grid_data = np.array(grid_img)
        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m'
        area_extent = (-5326849.0625, -5326849.0625,
                       5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = ImageFont.truetype(os.path.join(os.path.dirname(__file__),
                                               'test_data', 'DejaVuSerif.ttf'),
                                  10)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    font=font, fill='blue',
                    outline='blue', minor_outline='blue', write_text=False,
                    lon_placement='tblr', lat_placement='')

        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of nh grid failed')

    def test_add_polygon(self):
        from pycoast import ContourWriterPIL
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'nh_polygons.png'))
        grid_data = np.array(grid_img)
        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m'
        area_extent = (-5326849.0625, -5326849.0625,
                       5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        polygons = {
            'REYKJAVIK_ATC_A': ((-20.0, 73.0), (0.0, 73.0), (0.0, 61.0),
                                (-30.0, 61.0), (-39.0, 63.5), (-20, 70)),
            'REYKJAVIK_ATC_B': (
                (-39, 63.5), (-55 + 4 / 6.0, 63.5), (-57 + 45 / 60.0, 65),
                (-76, 76), (-75, 78), (-60, 82), (0, 90),
                (30, 82), (0, 82), (0, 73), (-20, 73), (-20, 70)),
            'REYKJAVIK_ATC':   (
                (0.0, 73.0), (0.0, 61.0), (-30.0, 61.0), (-39, 63.5),
                (-55 + 4 / 6.0, 63.5), (-57 + 45 / 60.0, 65),
                (-76, 76), (-75, 78), (-60, 82), (0, 90), (30, 82), (0, 82)),
            'ICELAND_BOX':     ((-25, 62.5), (-25, 67), (-13, 67), (-13, 62.5))
        }

        cw.add_polygon(img, area_def, polygons['REYKJAVIK_ATC'], outline='red')
        cw.add_polygon(img, area_def, polygons['ICELAND_BOX'], outline='green',
                       fill='gray')
        cw.add_coastlines(img, area_def, resolution='l', level=4)

        res = np.array(img)
        self.assertTrue(fft_metric(grid_data, res),
                        'Writing of nh polygons failed')

    def test_add_shapefile_shapes(self):
        from pycoast import ContourWriterPIL
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'brazil_shapefiles.png'))
        grid_data = np.array(grid_img)

        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m'
        area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterPIL(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_shapefile_shapes(img, area_def,
                                os.path.join(
                                    os.path.dirname(__file__),
                                    'test_data/shapes/Metareas.shp'),
                                outline='red')
        cw.add_shapefile_shape(img, area_def,
                               os.path.join(os.path.dirname(__file__),
                                            'test_data/shapes/divisao_politica/BR_Regioes.shp'), 3,
                               outline='blue')
        cw.add_shapefile_shape(img, area_def,
                               os.path.join(os.path.dirname(__file__),
                                            'test_data/shapes/divisao_politica/BR_Regioes.shp'), 4,
                               outline='blue', fill='green')

        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of Brazil shapefiles failed')


class TestPILAGG(TestPycoast):

    def test_europe_agg(self):
        from pycoast import ContourWriterAGG
        euro_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_europe_agg.png'))
        euro_data = np.array(euro_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_rivers(img, area_def, level=5, outline='blue', width=0.5,
                      outline_opacity=127)
        cw.add_borders(img, area_def, outline=(255, 0, 0),
                       width=3, outline_opacity=32)
        res = np.array(img)
        self.assertTrue(fft_metric(euro_data, res),
                        'Writing of contours failed for AGG')

    def test_europe_agg_file(self):
        from pycoast import ContourWriterAGG
        euro_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_europe_agg.png'))
        euro_data = np.array(euro_img)

        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines_to_file(test_file, area_def, resolution='l', level=4)
        cw.add_rivers_to_file(test_file, area_def, level=5, outline='blue',
                              width=0.5, outline_opacity=127)
        cw.add_borders_to_file(test_file, area_def, outline=(255, 0, 0),
                               width=3, outline_opacity=32)

        img = Image.open(test_file)
        res = np.array(img)
        self.assertTrue(fft_metric(euro_data, res),
                        'Writing of contours failed for AGG')

    def test_geos_agg(self):
        from pycoast import ContourWriterAGG
        geos_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'contours_geos_agg.png'))
        geos_data = np.array(geos_img)

        img = Image.new('RGB', (425, 425))
        proj4_string = \
            '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444,
                       5567248.074173444, 5570248.4773392612)
        # area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, (proj4_string, area_extent),
                          resolution='l', width=0.5)
        res = np.array(img)
        self.assertTrue(fft_metric(geos_data, res),
                        'Writing of geos contours failed for AGG')

    def test_grid_agg(self):
        from pycoast import ContourWriterAGG
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_europe_agg.png'))
        grid_data = np.array(grid_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0), write_text=False,
                    outline='blue', outline_opacity=255, width=1.0,
                    minor_outline='white', minor_outline_opacity=255,
                    minor_width=0.5, minor_is_tick=False)

        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of grid failed for AGG')

    def test_grid_agg_txt(self):
        from pycoast import ContourWriterAGG
        import aggdraw
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_europe_agg_txt.png'))
        grid_data = np.array(grid_img)

        img = Image.new('RGB', (640, 480))
        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = aggdraw.Font('blue', os.path.join(os.path.dirname(__file__),
                                                 'test_data',
                                                 'DejaVuSerif.ttf'), size=16,
                            opacity=200)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0), font=font,
                    outline='blue', outline_opacity=255, width=1.0,
                    minor_outline='white', minor_outline_opacity=255,
                    minor_width=0.5, minor_is_tick=False)

        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of grid failed for AGG')

    def test_grid_geos_agg(self):
        from pycoast import ContourWriterAGG
        geos_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_geos_agg.png'))
        geos_data = np.array(geos_img)
        img = Image.new('RGB', (425, 425))
        proj4_string = \
            '+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0'
        area_extent = (-5570248.4773392612, -5567248.074173444,
                       5567248.074173444, 5570248.4773392612)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution='l')
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    fill='blue', outline='blue', minor_outline='blue',
                    write_text=False)

        res = np.array(img)
        self.assertTrue(
            fft_metric(geos_data, res), 'Writing of geos contours failed')

    def test_grid_agg_file(self):
        from pycoast import ContourWriterAGG
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_europe_agg.png'))
        grid_data = np.array(grid_img)

        proj4_string = \
            '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines_to_file(grid_file, area_def, resolution='l', level=4)
        cw.add_grid_to_file(grid_file, area_def, (10.0, 10.0), (2.0, 2.0),
                            write_text=False, outline='blue',
                            outline_opacity=255, width=1.0,
                            minor_outline='white', minor_outline_opacity=255,
                            minor_width=0.5, minor_is_tick=False)
        img = Image.open(grid_file)
        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of grid failed for AGG')

    def test_grid_nh_agg(self):
        from pycoast import ContourWriterAGG
        import aggdraw
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'grid_nh_agg.png'))
        grid_data = np.array(grid_img)
        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m'
        area_extent = (-5326849.0625, -5326849.0625,
                       5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        font = aggdraw.Font('blue', os.path.join(os.path.dirname(__file__),
                                                 'test_data',
                                                 'DejaVuSerif.ttf'), size=10)
        cw.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0),
                    font=font, fill='blue',
                    outline='blue', minor_outline='blue',
                    lon_placement='tblr', lat_placement='')

        res = np.array(img)

        # NOTE: Experience inconsistency in ttf font writing between systems.
        # Still trying to figure out why this test sometimes fails to write
        # correct font markings.
        self.assertTrue(fft_metric(grid_data, res),
                        'Writing of nh grid failed for AGG')

    def test_add_polygon_agg(self):
        from pycoast import ContourWriterAGG
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'nh_polygons_agg.png'))
        grid_data = np.array(grid_img)

        img = Image.new('RGB', (425, 425))
        proj4_string = '+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m'
        area_extent = (-5326849.0625, -5326849.0625,
                       5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        polygons = {
            'REYKJAVIK_ATC_A': ((-20.0, 73.0), (0.0, 73.0), (0.0, 61.0),
                                (-30.0, 61.0), (-39.0, 63.5), (-20, 70)),
            'REYKJAVIK_ATC_B': ((-39, 63.5), (-55 + 4 / 6.0, 63.5),
                                (-57 + 45 / 60.0, 65), (-76, 76),
                                (-75, 78), (-60, 82), (0, 90), (30, 82),
                                (0, 82), (0, 73), (-20, 73), (-20, 70)),
            'REYKJAVIK_ATC': ((0.0, 73.0), (0.0, 61.0), (-30.0, 61.0),
                              (-39, 63.5), (-55 + 4 / 6.0, 63.5),
                              (-57 + 45 / 60.0, 65), (-76, 76), (-75, 78),
                              (-60, 82), (0, 90), (30, 82), (0, 82)),
            'ICELAND_BOX': ((-25, 62.5), (-25, 67), (-13, 67), (-13, 62.5))
        }

        cw.add_polygon(img, area_def, polygons['REYKJAVIK_ATC'],
                       outline='red', width=2)
        cw.add_polygon(img, area_def, polygons['ICELAND_BOX'],
                       outline='green', fill='gray', width=2)
        cw.add_coastlines(img, area_def, resolution='l', level=4)

        res = np.array(img)
        self.assertTrue(fft_metric(grid_data, res),
                        'Writing of nh polygons failed')

    def test_add_shapefile_shapes_agg(self):
        from pycoast import ContourWriterAGG
        grid_img = Image.open(os.path.join(os.path.dirname(__file__),
                                           'brazil_shapefiles_agg.png'))
        grid_data = np.array(grid_img)

        img = Image.new('RGB', (425, 425))
        proj4_string = \
            '+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m'
        area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution='l', level=4)
        cw.add_shapefile_shapes(img, area_def,
                                os.path.join(
                                    os.path.dirname(__file__),
                                    'test_data/shapes/Metareas.shp'),
                                outline='red', width=2)
        cw.add_shapefile_shape(img, area_def,
                               os.path.join(os.path.dirname(__file__),
                                            'test_data/shapes/divisao_politica/BR_Regioes.shp'), 3,
                               outline='blue')
        cw.add_shapefile_shape(img, area_def,
                               os.path.join(os.path.dirname(__file__),
                                            'test_data/shapes/divisao_politica/BR_Regioes.shp'), 4,
                               outline='blue', fill='green')

        res = np.array(img)
        self.assertTrue(
            fft_metric(grid_data, res), 'Writing of Brazil shapefiles failed')


def suite():
    loader = unittest.TestLoader()
    mysuite = unittest.TestSuite()
    mysuite.addTest(loader.loadTestsFromTestCase(TestPIL))
    mysuite.addTest(loader.loadTestsFromTestCase(TestPILAGG))

    return mysuite


if __name__ == "__main__":
    suite()
