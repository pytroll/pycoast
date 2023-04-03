#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011-2020 PyCoast Developers
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
"""Main unit tests for pycoast."""

import os
import pathlib
import time
from glob import glob

import aggdraw
import numpy as np
import pytest
import shapefile
from PIL import Image, ImageFont
from pyresample.geometry import AreaDefinition
from pytest_lazyfixture import lazy_fixture

LOCAL_DIR = os.path.dirname(__file__)

gshhs_root_dir = os.path.join(LOCAL_DIR, "test_data", "gshhs")
repos_root_dir = os.path.join(LOCAL_DIR, "..", "..")
test_filename = "test_image.png"
grid_filename = "test_grid.png"
p_coasts_filename = "test_coasts_p_mode.png"

font_path = os.path.join(LOCAL_DIR, "test_data", "DejaVuSerif.ttf")
agg_font_20_yellow = aggdraw.Font("yellow", font_path, opacity=255, size=20)
agg_font_20_orange = aggdraw.Font("orange", font_path, opacity=255, size=20)
pil_font_20 = ImageFont.truetype(font_path, 20)
pil_font_16 = ImageFont.truetype(font_path, 16)


@pytest.fixture
def cd_test_dir(monkeypatch):
    """Change directory to the pycoast/tests directory."""
    monkeypatch.chdir(LOCAL_DIR)


def fft_proj_rms(a1, a2):
    """Compute the RMS of differences between FFT vectors of a1 and projection of FFT vectors of a2.

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
        adjusted_ps2 = ps2 * np.cos(theta)
        ms += ((adjusted_ps2 - ps1) ** 2).sum() / float(ps1.size)

    rms = np.sqrt(ms)

    return rms


def fft_metric(data1, data2, max_value=0.1, plot_failure=False):
    """Execute FFT metric."""
    rms = fft_proj_rms(data1, data2)
    within_threshold = rms <= max_value
    if not within_threshold and plot_failure:
        import matplotlib.pyplot as plt

        fig, (ax1, ax2, ax3) = plt.subplots(3, 1)
        ax1.imshow(data1)
        ax2.imshow(data2)
        ax3.imshow(np.abs(data1.astype(np.float64) - data2.astype(np.float64)).astype(np.uint8))
        plt.show()
    return within_threshold


def europe():
    """Create a Europe area."""
    proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
    area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
    area_def = (proj4_string, area_extent)
    return area_def


EUROPE = europe()


def geos():
    """Create a geos area."""
    proj4_string = "+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0"
    area_extent = (
        -5570248.4773392612,
        -5567248.074173444,
        5567248.074173444,
        5570248.4773392612,
    )
    area_def = (proj4_string, area_extent)
    return area_def


GEOS = geos()


def germ():
    """Create an area covering Germany."""
    proj4_string = "+proj=stere +ellps=bessel +lat_0=90.0 +lon_0=5.0 +lat_ts=50.0 +a=6378144.0 +b=6356759.0"
    area_extent = [-155100.436345, -4441495.37946, 868899.563655, -3417495.37946]
    area_def = (proj4_string, area_extent)
    return area_def


GERM = germ()


def dateline_1():
    """Create an area covering the dateline."""
    proj4_string = "+proj=stere +lon_0=-170.00 +lat_0=60.00 +lat_ts=50.00 +ellps=WGS84"
    area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
    area_def = (proj4_string, area_extent)
    return area_def


DATELINE_1 = dateline_1()


def dateline_2():
    """Create another area covering the dateline."""
    proj4_string = "+proj=stere +lon_0=140.00 +lat_0=60.00 +lat_ts=50.00 +ellps=WGS84"
    area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
    area_def = (proj4_string, area_extent)
    return area_def


DATELINE_2 = dateline_2()


def nh():
    """Create the nh area."""
    proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
    area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
    area_def = (proj4_string, area_extent)
    return area_def


NH = nh()


def nh_1024():
    """Create the nh area in 1024x1024 pixels."""
    proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
    area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
    area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)
    return area_def


NH_1024 = nh_1024()


def nh_def(shape):
    """Create the nh area definition with custom shape."""
    proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
    area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 0.0)
    area_def = AreaDefinition("nh", "nh", "nh", proj4_string, shape[1], shape[0], area_extent)
    return area_def


def brazil():
    """Create a Brazil area."""
    proj4_string = "+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m"
    area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
    area_def = (proj4_string, area_extent)
    return area_def


BRAZIL = brazil()


def nh_425():
    """Create the nh area in 425 pixels size."""
    proj_dict = {
        "proj": "laea",
        "lat_0": 90.0,
        "lon_0": 0.0,
        "a": 6371228.0,
        "units": "m",
    }
    area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
    area_def = AreaDefinition("nh", "nh", "nh", proj_dict, 425, 425, area_extent)
    return area_def


NH_425 = nh_425()


def eurasia():
    """Create a Eurasia area."""
    proj4_string = "+proj=tmerc +ellps=WGS84 +lat_0=20.0 +lon_0=-50.0"
    area_extent = [-4865942.5, 1781111.9, 4865942.5, 7235767.2]
    area_def = (proj4_string, area_extent)
    return area_def


def north_atlantic():
    """Create a North Atlantic area."""
    proj4_string = "+proj=tmerc +ellps=WGS84 +lat_0=20.0 +lon_0=50.0"
    area_extent = [-4865942.5, 1781111.9, 4865942.5, 7235767.2]
    area_def = (proj4_string, area_extent)
    return area_def


def uk_and_ireland():
    """Create an area covering Ireland and the UK."""
    proj4_string = "+proj=stere +ellps=WGS84 +lon_0=-4.532 +lat_0=54.228"
    area_extent = (-600000.0, -600000.0, 600000.0, 600000.0)
    area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 800, 800, area_extent)
    return area_def


def south_america():
    """Create an area covering south America."""
    proj4_string = "+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m"
    area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
    area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 425, 425, area_extent)
    return area_def


def europe_1024():
    """Create a Europe area in 1024 pixels size."""
    proj4_string = "+proj=stere +ellps=WGS84 +lat_0=51.5 +lon_0=-0.1"
    area_extent = (-2000000.0, -2000000.0, 2000000.0, 2000000.0)
    area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)
    return area_def


def bering_straight():
    """Create an area covering the Bering straight."""
    proj4_string = "+proj=merc +ellps=WGS84 +lon_0=170.0"
    area_extent = [-3899875.0, 1795000.0, 5014125.0, 9600000.0]
    area_def = (proj4_string, area_extent)
    return area_def


def hawaii():
    """Create an area covering Hawai."""
    proj4_string = "+proj=tmerc +ellps=WGS84 +lon_0=-155.0"
    area_extent = [-3503550.0, -556597.5, 3503550.0, 5009377.3]
    area_def = (proj4_string, area_extent)
    return area_def


@pytest.fixture(scope="session")
def cw_pil():
    """Create a PIL ContourWriter."""
    from pycoast import ContourWriterPIL

    cw = ContourWriterPIL(gshhs_root_dir)
    return cw


@pytest.fixture
def cw_agg():
    """Create a PIL ContourWriter."""
    from pycoast import ContourWriterAGG

    cw = ContourWriterAGG(gshhs_root_dir)
    return cw


@pytest.fixture
def new_test_image(request, tmp_path):
    """Create a new test image, and save it to tmp_path if the test fails."""
    img_container = []

    def _new_test_image(mode, shape, filename, color=0):
        img = Image.new(mode, shape, color=color)
        img_container.append(filename)
        img_container.append(img)
        return img

    yield _new_test_image
    try:
        filename, img = img_container
    except ValueError:  # the fixture wasn't used
        return
    if request.node.rep_call.failed:
        img.save(os.path.join(tmp_path, filename))


def images_match(ref_image, test_image):
    """Check is images match."""
    return fft_metric(np.array(ref_image), np.array(test_image))


class _ContourWriterTestBase:
    """Base class for test classes that need example images."""

    def setup_method(self):
        img = Image.new("RGB", (640, 480))
        img.save(test_filename)
        img.save(grid_filename)
        img_p = Image.new("P", (640, 480))
        img_p.save(p_coasts_filename)

    def teardown_method(self):
        os.remove(test_filename)
        os.remove(grid_filename)
        os.remove(p_coasts_filename)


@pytest.fixture
def test_file_path(tmp_path):
    """Create a test image file on disk."""
    path = tmp_path / test_filename
    img = Image.new("RGB", (640, 480))
    img.save(path)
    yield path


@pytest.fixture
def grid_file_path(tmp_path):
    """Create a test grid image file on disk."""
    path = tmp_path / grid_filename
    img = Image.new("RGB", (640, 480))
    img.save(path)
    yield path


class TestContourWriterPIL:
    """Test PIL-based contour writer."""

    def test_europe_coastlines_rivers_and_borders(self, cw_pil, new_test_image):
        """Test coastlines, rivers and borders over Europe."""
        filename = "contours_europe.png"

        euro_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (640, 480), filename)

        area_def = EUROPE
        cw_pil.add_coastlines(img, area_def, resolution="l", level=4)
        cw_pil.add_rivers(img, area_def, level=5, outline="blue")
        cw_pil.add_borders(img, area_def, outline=(255, 0, 0), level=1)

        assert images_match(euro_img, img), "Writing of contours failed"

    def test_europe_coastlines_rivers_and_borders_on_file(self, cw_pil, test_file_path):
        """Test coastlines, rivers and borders over Europe on a file."""
        filename = "contours_europe.png"
        euro_img = Image.open(os.path.join(LOCAL_DIR, filename))

        area_def = EUROPE
        cw_pil.add_coastlines_to_file(test_file_path, area_def, resolution="l", level=4)
        cw_pil.add_rivers_to_file(test_file_path, area_def, level=5, outline="blue")
        cw_pil.add_borders_to_file(test_file_path, area_def, outline=(255, 0, 0))

        img = Image.open(test_file_path)
        assert images_match(euro_img, img), "Writing of contours to file failed"

    def test_geos(self, cw_pil, new_test_image):
        filename = "contours_geos.png"
        geos_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = GEOS
        cw_pil.add_coastlines(img, area_def, resolution="l")

        assert images_match(geos_img, img), "Writing of geos contours failed"

    @pytest.mark.parametrize(
        "filename, shape, area_def, level, grid_kwargs",
        [
            (
                "grid_europe.png",
                (640, 480),
                EUROPE,
                4,
                dict(fill="blue", write_text=False, outline="blue", minor_outline="blue"),
            ),
            (
                "grid_geos.png",
                (425, 425),
                GEOS,
                1,
                dict(fill="blue", write_text=False, outline="blue", minor_outline="blue"),
            ),
            (
                "grid_geos.png",
                (425, 425),
                GEOS,
                1,
                dict(fill="blue", write_text=True, outline="blue", minor_outline="blue"),
            ),
            (
                "grid_germ.png",
                (1024, 1024),
                GERM,
                4,
                dict(fill="yellow", write_text=True, outline="red", minor_outline="white"),
            ),
            (
                "dateline_cross.png",
                (640, 480),
                DATELINE_1,
                4,
                dict(
                    fill="blue",
                    write_text=False,
                    outline="blue",
                    minor_outline="blue",
                    lon_placement="b",
                    lat_placement="lr",
                ),
            ),
            (
                "dateline_boundary_cross.png",
                (640, 480),
                DATELINE_2,
                4,
                dict(
                    fill="blue",
                    write_text=False,
                    outline="blue",
                    minor_outline="blue",
                    lon_placement="b",
                    lat_placement="lr",
                ),
            ),
        ],
    )
    @pytest.mark.usefixtures("cd_test_dir")
    def test_grid(self, cw_pil, new_test_image, filename, shape, area_def, level, grid_kwargs):
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", shape, filename)

        cw_pil.add_coastlines(img, area_def, resolution="l", level=level)
        font = pil_font_16
        cw_pil.add_grid(img, area_def, (10.0, 10.0), (2.0, 2.0), font=font, **grid_kwargs)

        assert images_match(grid_img, img), "Writing of grid failed"

    def test_grid_nh(self, cw_pil, new_test_image):
        filename = "grid_nh.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = NH

        cw_pil.add_coastlines(img, area_def, resolution="l", level=4)
        font = ImageFont.truetype(font_path, 10)
        cw_pil.add_grid(
            img,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            font=font,
            fill="blue",
            outline="blue",
            minor_outline="blue",
            write_text=False,
            lon_placement="tblr",
            lat_placement="",
        )

        assert images_match(grid_img, img), "Writing of nh grid failed"

    def test_grid_file(self, cw_pil, grid_file_path):
        filename = "grid_europe.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        area_def = EUROPE

        cw_pil.add_coastlines_to_file(grid_file_path, area_def, resolution="l", level=4)
        font = pil_font_16
        cw_pil.add_grid_to_file(
            grid_file_path,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            font=font,
            fill="blue",
            write_text=False,
            outline="blue",
            minor_outline="blue",
        )

        img = Image.open(grid_file_path)
        assert images_match(grid_img, img), "Writing of grid failed"

    def test_add_polygon(self, cw_pil, new_test_image):
        filename = "nh_polygons.png"
        ref_image = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = NH

        polygons = {
            "REYKJAVIK_ATC_A": (
                (-20.0, 73.0),
                (0.0, 73.0),
                (0.0, 61.0),
                (-30.0, 61.0),
                (-39.0, 63.5),
                (-20, 70),
            ),
            "REYKJAVIK_ATC_B": (
                (-39, 63.5),
                (-55 + 4 / 6.0, 63.5),
                (-57 + 45 / 60.0, 65),
                (-76, 76),
                (-75, 78),
                (-60, 82),
                (0, 90),
                (30, 82),
                (0, 82),
                (0, 73),
                (-20, 73),
                (-20, 70),
            ),
            "REYKJAVIK_ATC": (
                (0.0, 73.0),
                (0.0, 61.0),
                (-30.0, 61.0),
                (-39, 63.5),
                (-55 + 4 / 6.0, 63.5),
                (-57 + 45 / 60.0, 65),
                (-76, 76),
                (-75, 78),
                (-60, 82),
                (0, 90),
                (30, 82),
                (0, 82),
            ),
            "ICELAND_BOX": ((-25, 62.5), (-25, 67), (-13, 67), (-13, 62.5)),
        }

        cw_pil.add_polygon(img, area_def, polygons["REYKJAVIK_ATC"], outline="red")
        cw_pil.add_polygon(img, area_def, polygons["ICELAND_BOX"], outline="green", fill="gray")
        cw_pil.add_coastlines(img, area_def, resolution="l", level=4)

        assert images_match(ref_image, img), "Writing of nh polygons failed"

    def test_add_points_pil(self, cw_pil, new_test_image):
        filename = "nh_points_pil.png"
        ref_image = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (1024, 1024), filename, color=(255, 255, 255))

        area_def = NH_1024

        cw_pil.add_coastlines(img, area_def, outline="black", resolution="l", level=4)
        cw_pil.add_borders(img, area_def, outline="black", level=1, resolution="c")

        points_list = [((13.4050, 52.5200), "Berlin")]
        cw_pil.add_points(
            img,
            area_def,
            points_list=points_list,
            font_file=font_path,
            symbol="asterisk",
            ptsize=6,
            outline="red",
            box_outline="black",
        )

        points_list = [((12.4964, 41.9028), "Rome")]
        cw_pil.add_points(
            img,
            area_def,
            points_list=points_list,
            font_file=font_path,
            symbol="square",
            ptsize=6,
            outline="blue",
            fill="yellow",
            box_outline="black",
        )

        assert images_match(ref_image, img), "Writing of nh points failed"

    def test_add_points_coordinate_conversion(self, cw_pil, new_test_image):
        """Check that a point added with lonlat coordinates matches the same point in pixel coordinates."""
        shape = (512, 1024)
        area_def = nh_def(shape)
        lonlat_coords = (13.4050, 52.5200)
        pixel_colrow = area_def.get_array_indices_from_lonlat(*lonlat_coords)
        negative_pixel_colrow = (pixel_colrow[0] - shape[1], pixel_colrow[1] - shape[0])

        img1 = Image.new("RGB", shape[::-1], (255, 255, 255))

        points_list = [(lonlat_coords, "Berlin")]
        cw_pil.add_points(
            img1,
            area_def,
            points_list=points_list,
            font_file=font_path,
            symbol="asterisk",
            ptsize=6,
            outline="red",
        )
        res1 = np.array(img1)
        assert (res1[..., 0] == 255).all()  # everything is either white or red
        assert (res1[..., 1] != 255).any()  # not a completely white/empty image
        assert (res1[..., 2] != 255).any()  # not a completely white/empty image

        for img_coords in (pixel_colrow, negative_pixel_colrow):
            img2 = Image.new("RGB", shape[::-1], (255, 255, 255))
            points_list = [(img_coords, "Berlin")]
            cw_pil.add_points(
                img2,
                area_def,
                points_list=points_list,
                font_file=font_path,
                symbol="asterisk",
                ptsize=6,
                outline="red",
                coord_ref="image",
            )
            res2 = np.array(img2)
            assert (res2 != 255).any()  # not a completely black/empty image
            np.testing.assert_allclose(res1, res2)

    def test_add_points_bad_image_coords(self, cw_pil):
        shape = (512, 1024)
        area_def = nh_def(shape)
        for pixel_colrow in (shape[::-1], (-10000, -10000)):
            img1 = Image.new("RGB", shape[::-1], (255, 255, 255))

            points_list = [(pixel_colrow, "Berlin")]
            cw_pil.add_points(
                img1,
                area_def,
                points_list=points_list,
                font_file=font_path,
                symbol="asterisk",
                ptsize=6,
                outline="red",
                coord_ref="image",
            )
            res1 = np.array(img1)
            np.testing.assert_allclose(res1, 255)  # no added points

    def test_add_points_bad_coord_ref(self, cw_pil):
        shape = (512, 1024)
        area_def = nh_def(shape)
        img1 = Image.new("RGB", shape[::-1], (255, 255, 255))
        points_list = [((0, 0), "Berlin")]
        with pytest.raises(ValueError):
            cw_pil.add_points(
                img1,
                area_def,
                points_list=points_list,
                font_file=font_path,
                symbol="asterisk",
                ptsize=6,
                outline="red",
                coord_ref="fake",
            )

    def test_add_shapefile_shapes(self, cw_pil, new_test_image):
        filename = "brazil_shapefiles.png"
        ref_image = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = BRAZIL

        cw_pil.add_coastlines(img, area_def, resolution="l", level=4)
        cw_pil.add_shapefile_shapes(
            img,
            area_def,
            os.path.join(LOCAL_DIR, "test_data/shapes/Metareas.shp"),
            outline="red",
        )
        cw_pil.add_shapefile_shape(
            img,
            area_def,
            os.path.join(
                LOCAL_DIR,
                "test_data/shapes/divisao_politica/BR_Regioes.shp",
            ),
            3,
            outline="blue",
        )
        cw_pil.add_shapefile_shape(
            img,
            area_def,
            os.path.join(
                LOCAL_DIR,
                "test_data/shapes/divisao_politica/BR_Regioes.shp",
            ),
            4,
            outline="blue",
            fill="green",
        )

        assert images_match(ref_image, img), "Writing of Brazil shapefiles failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_config_file_coasts_and_grid(self, cw_pil, new_test_image):
        overlay_config = os.path.join(LOCAL_DIR, "coasts_and_grid.ini")
        filename = "grid_nh.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        area_def = NH_425

        overlay = cw_pil.add_overlay_from_config(overlay_config, area_def)
        img = new_test_image("RGB", (425, 425), filename)
        img.paste(overlay, mask=overlay)
        assert images_match(grid_img, img), "Writing of nh grid failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_config_file_points_and_borders_pil(self, cw_pil, new_test_image):
        config_file = os.path.join(LOCAL_DIR, "nh_points_pil.ini")
        filename = "nh_points_cfg_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (1024, 1024), filename, color=(255, 255, 255))

        area_def = NH_1024

        cw_pil.add_overlay_from_config(config_file, area_def, img)

        assert images_match(grid_img, img), "Writing of nh points failed"

    def test_add_cities_pil(self, cw_pil, new_test_image):
        filename = "nh_cities_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (1024, 1024), filename, color=(255, 255, 255))

        area_def = NH_1024

        cw_pil.add_coastlines(img, area_def, outline="black", resolution="l", level=4)
        cw_pil.add_borders(img, area_def, outline="black", level=1, resolution="c")

        cities_list = ["Zurich", "Oslo", "Reykjavik", "Fairbanks", "Toronto"]
        cw_pil.add_cities(
            img,
            area_def,
            cities_list=cities_list,
            font_file=font_path,
            font_size=20,
            symbol="square",
            ptsize=16,
            outline="black",
            width=3,
            fill="blue",
            fill_opacity=128,
            box_outline="blue",
            box_linewidth=0.5,
            box_fill="yellow",
            box_opacity=200,
        )

        assert images_match(grid_img, img), "Writing of nh cities_pil failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_add_cities_cfg_pil(self, cw_pil, new_test_image):
        config_file = os.path.join(LOCAL_DIR, "nh_cities_pil.ini")
        filename = "nh_cities_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (1024, 1024), filename, color=(255, 255, 255))

        area_def = NH_1024

        cw_pil.add_overlay_from_config(config_file, area_def, img)

        assert images_match(grid_img, img), "Writing of nh cities_cfg_pil failed"

    def test_add_cities_from_dict_pil(self, cw_pil, new_test_image):
        filename = "nh_cities_from_dict_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (1024, 1024), filename, color=(255, 255, 255))

        area_def = europe_1024()

        overlays = {}
        overlays["coasts"] = {"level": 4, "resolution": "l", "outline": "black"}
        overlays["borders"] = {"level": 1, "outline": "black", "resolution": "c"}
        cities_list1 = [
            "Berlin",
            "Paris",
            "London",
            "Dublin",
            "Madrid",
            "Reykjavik",
            "Oslo",
            "Rome",
        ]
        cities_list2 = ["Freiburg", "Montelimar", "Huesca", "Marseille"]
        cities_list3 = ["Belp", "Bad Schwalbach", "Edinburgh", "Hilversum"]
        cities_type1 = {
            "cities_list": cities_list1,
            "font": font_path,
            "font_size": 26,
            "symbol": "circle",
            "ptsize": 24,
            "outline": "black",
            "fill": "red",
            "box_outline": "black",
        }
        cities_type2 = {
            "cities_list": cities_list2,
            "font": font_path,
            "font_size": 24,
            "symbol": "pentagon",
            "ptsize": 24,
            "outline": "red",
            "fill": "blue",
            "box_outline": "black",
        }
        cities_type3 = {
            "cities_list": cities_list3,
            "font": font_path,
            "font_size": 22,
            "symbol": "star5",
            "ptsize": 35,
            "outline": "green",
            "fill": "yellow",
            "box_outline": "black",
        }

        overlays["cities"] = [cities_type1, cities_type2, cities_type3]

        img = cw_pil.add_overlay_from_dict(overlays, area_def, background=img)

        assert images_match(grid_img, img), "Writing of nh_cities_from_dict_pil failed"

    def test_add_shapefiles_from_dict_pil(self, cw_pil, new_test_image):
        filename = "two_shapefiles_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = south_america()

        overlays = {}
        overlays["coasts"] = {"level": 4, "resolution": "l"}
        overlays["shapefiles"] = [
            {
                "filename": os.path.join(LOCAL_DIR, "test_data/shapes/Metareas.shp"),
                "outline": "magenta",
                "width": 2.5,
            },
            {
                "filename": os.path.join(
                    LOCAL_DIR,
                    "test_data/shapes/divisao_politica/BR_Regioes.shp",
                ),
                "outline": "red",
                "fill": "yellow",
                "fill_opacity": 50,
            },
        ]

        img = cw_pil.add_overlay_from_dict(overlays, area_def, background=img)

        assert images_match(grid_img, img), "Writing of two shapefiles from dict pil failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_add_one_shapefile_from_cfg_pil(self, cw_pil, new_test_image):
        config_file = os.path.join(LOCAL_DIR, "nh_one_shapefile.ini")
        filename = "one_shapefile_from_cfg_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (425, 425), filename)
        area_def = south_america()

        cw_pil.add_overlay_from_config(config_file, area_def, img)

        assert images_match(grid_img, img), "Writing one shapefile from cfg pil failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_add_grid_from_dict_pil(self, cw_pil, new_test_image):
        filename = "grid_from_dict_pil.png"
        grid_img = Image.open(os.path.join(LOCAL_DIR, filename))

        img = new_test_image("RGB", (800, 800), filename)
        area_def = uk_and_ireland()

        font = ImageFont.truetype(font_path, 40)

        overlays = {}
        overlays["coasts"] = {"width": 3.0, "level": 4, "resolution": "l"}
        overlays["grid"] = {
            "major_lonlat": (5, 5),
            "minor_lonlat": (1, 1),
            "outline": (255, 0, 0),
            "outline_opacity": 127,
            "minor_outline": (0, 0, 255),
            "minor_outline_opacity": 127,
            "width": 10.5,
            "minor_width": 5.5,
            "minor_is_tick": False,
            "write_text": True,
            "lat_placement": "lr",
            "lon_placement": "b",
            "font": font,
            "fill": "yellow",
        }
        # Fill is pil text color! Pil Font can be None, then a of default font is choosen

        img = cw_pil.add_overlay_from_dict(overlays, area_def, background=img)

        assert images_match(grid_img, img), "Writing grid from dict pil failed"


@pytest.mark.parametrize(
    "cw, filename, area_def, specific_kwargs",
    [
        (
            lazy_fixture("cw_pil"),
            "western_shapes_pil.png",
            north_atlantic(),
            dict(font=pil_font_16, fill="yellow", outline="red", minor_outline="red"),
        ),
        (
            lazy_fixture("cw_pil"),
            "eastern_shapes_pil.png",
            eurasia(),
            dict(font=pil_font_20, fill="yellow", outline="red", minor_outline="red"),
        ),
        (
            lazy_fixture("cw_agg"),
            "western_shapes_agg.png",
            north_atlantic(),
            dict(
                font=agg_font_20_orange,
                outline="blue",
                width=5.0,
                outline_opacity=100,
                minor_outline="blue",
                minor_width=5.0,
                minor_outline_opacity=200,
            ),
        ),
        (
            lazy_fixture("cw_agg"),
            "eastern_shapes_agg.png",
            eurasia(),
            dict(
                font=agg_font_20_orange,
                outline="blue",
                width=5.0,
                outline_opacity=100,
                minor_outline="blue",
                minor_width=5.0,
                minor_outline_opacity=200,
            ),
        ),
    ],
)
def test_shapes(new_test_image, cw, filename, area_def, specific_kwargs):
    """Test western/eastern shapes."""
    result_file = os.path.join(LOCAL_DIR, filename)
    grid_img = Image.open(result_file)

    img = new_test_image("RGB", (1000, 560), filename)

    cw.add_coastlines(img, area_def, resolution="l", level=2)

    cw.add_grid(
        img,
        area_def,
        (10.0, 10.0),
        (5.0, 5.0),
        write_text=True,
        lon_placement="lbr",
        lat_placement="",
        **specific_kwargs,
    )

    assert images_match(grid_img, img), "Writing of shapes failed"


@pytest.mark.parametrize(
    "cw, filename, shape, area_def, specific_kwargs",
    [
        (  # lon=175 +/-40, lat=16..65 | Avoid Eurasia scratch with asymmetric area_extent
            lazy_fixture("cw_pil"),
            "no_h_scratch_pil.png",
            (888, 781),
            bering_straight(),
            dict(font=pil_font_20, fill="yellow", outline="orange", minor_outline="orange"),
        ),
        (  # lon=155+/-30 lat=-5..45 | No Eurasia problem (Eurasia has always lat > 0.0)
            lazy_fixture("cw_pil"),
            "no_v_scratch_pil.png",
            (888, 705),
            hawaii(),
            dict(font=pil_font_20, fill="yellow", outline="orange", minor_outline="orange"),
        ),
        (  # lon=175 +/-40, lat=16..65 | Avoid Eurasia scratch with asymmetric area_extent
            lazy_fixture("cw_agg"),
            "no_h_scratch_agg.png",
            (888, 781),
            bering_straight(),
            dict(
                font=agg_font_20_yellow,
                outline="green",
                width=5.0,
                outline_opacity=100,
                minor_outline="green",
                minor_width=5.0,
                minor_outline_opacity=200,
            ),
        ),
        (  # lon=155+/-30 lat=-5..45 | No Eurasia problem (Eurasia has always lat > 0.0)
            lazy_fixture("cw_agg"),
            "no_v_scratch_agg.png",
            (888, 705),
            hawaii(),
            dict(
                font=agg_font_20_yellow,
                outline="green",
                width=5.0,
                outline_opacity=100,
                minor_outline="green",
                minor_width=5.0,
                minor_outline_opacity=200,
            ),
        ),
    ],
)
@pytest.mark.usefixtures("cd_test_dir")
def test_no_scratch(new_test_image, cw, filename, shape, area_def, specific_kwargs):
    """Test no scratches are visible."""
    result_file = os.path.join(LOCAL_DIR, filename)
    grid_img = Image.open(result_file)

    img = new_test_image("RGB", shape, filename)

    cw.add_coastlines(img, area_def, resolution="l", level=2)

    cw.add_grid(
        img,
        area_def,
        (10.0, 10.0),
        (5.0, 5.0),
        write_text=True,
        lon_placement="bt",
        lat_placement="lr",
        **specific_kwargs,
    )
    cw.add_rivers(img, area_def, level=5, outline="blue")
    cw.add_borders(img, area_def, outline="red")

    assert images_match(grid_img, img), "Writing of no_scratch failed"


class TestContourWriterAGG(_ContourWriterTestBase):
    """Test AGG contour writer."""

    def test_europe_agg(self):
        from pycoast import ContourWriterAGG

        euro_img = Image.open(os.path.join(LOCAL_DIR, "contours_europe_agg.png"))
        euro_data = np.array(euro_img)

        img = Image.new("RGB", (640, 480))
        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution="l", level=4)
        cw.add_rivers(img, area_def, level=5, outline="blue", width=0.5, outline_opacity=127)
        cw.add_borders(img, area_def, outline=(255, 0, 0), width=3, outline_opacity=32)
        res = np.array(img)
        assert fft_metric(euro_data, res), "Writing of contours failed for AGG"

    def test_europe_agg_file(self):
        from pycoast import ContourWriterAGG

        euro_img = Image.open(os.path.join(LOCAL_DIR, "contours_europe_agg.png"))
        euro_data = np.array(euro_img)

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines_to_file(test_filename, area_def, resolution="l", level=4)
        cw.add_rivers_to_file(test_filename, area_def, level=5, outline="blue", width=0.5, outline_opacity=127)
        cw.add_borders_to_file(test_filename, area_def, outline=(255, 0, 0), width=3, outline_opacity=32)

        img = Image.open(test_filename)
        res = np.array(img)
        assert fft_metric(euro_data, res), "Writing of contours failed for AGG"

    def test_geos_agg(self):
        from pycoast import ContourWriterAGG

        geos_img = Image.open(os.path.join(LOCAL_DIR, "contours_geos_agg.png"))
        geos_data = np.array(geos_img)

        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0"
        area_extent = (
            -5570248.4773392612,
            -5567248.074173444,
            5567248.074173444,
            5570248.4773392612,
        )
        # area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, (proj4_string, area_extent), resolution="l", width=0.5)
        res = np.array(img)
        assert fft_metric(geos_data, res), "Writing of geos contours failed for AGG"

    def test_grid_agg(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_europe_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (640, 480))
        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution="l", level=4)
        cw.add_grid(
            img,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            write_text=False,
            outline="blue",
            outline_opacity=255,
            width=1.0,
            minor_outline="white",
            minor_outline_opacity=255,
            minor_width=0.5,
            minor_is_tick=False,
        )

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of grid failed for AGG"

    def test_grid_agg_txt(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_europe_agg_txt.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (640, 480))
        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution="l", level=4)
        font = aggdraw.Font(
            "blue",
            font_path,
            size=16,
            opacity=200,
        )
        cw.add_grid(
            img,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            font=font,
            outline="blue",
            outline_opacity=255,
            width=1.0,
            minor_outline="white",
            minor_outline_opacity=255,
            minor_width=0.5,
            minor_is_tick=False,
        )

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of grid failed for AGG"

    def test_grid_geos_agg(self):
        from pycoast import ContourWriterAGG

        geos_img = Image.open(os.path.join(LOCAL_DIR, "grid_geos_agg.png"))
        geos_data = np.array(geos_img)
        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=geos +lon_0=0.0 +a=6378169.00 +b=6356583.80 +h=35785831.0"
        area_extent = (
            -5570248.4773392612,
            -5567248.074173444,
            5567248.074173444,
            5570248.4773392612,
        )
        area_def = (proj4_string, area_extent)
        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, resolution="l")
        cw.add_grid(
            img,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            fill="blue",
            outline="blue",
            minor_outline="blue",
            write_text=False,
        )

        res = np.array(img)
        assert fft_metric(geos_data, res), "Writing of geos contours failed"

    def test_grid_agg_file(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_europe_agg.png"))
        grid_data = np.array(grid_img)

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines_to_file(grid_filename, area_def, resolution="l", level=4)
        cw.add_grid_to_file(
            grid_filename,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            write_text=False,
            outline="blue",
            outline_opacity=255,
            width=1.0,
            minor_outline="white",
            minor_outline_opacity=255,
            minor_width=0.5,
            minor_is_tick=False,
        )
        img = Image.open(grid_filename)
        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of grid failed for AGG"

    def test_grid_nh_agg(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_nh_agg.png"))
        grid_data = np.array(grid_img)
        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution="l", level=4)
        font = aggdraw.Font(
            "blue",
            font_path,
            size=10,
        )
        cw.add_grid(
            img,
            area_def,
            (10.0, 10.0),
            (2.0, 2.0),
            font=font,
            fill="blue",
            outline="blue",
            minor_outline="blue",
            lon_placement="tblr",
            lat_placement="",
        )

        res = np.array(img)

        # NOTE: Experience inconsistency in ttf font writing between systems.
        # Still trying to figure out why this test sometimes fails to write
        # correct font markings.
        assert fft_metric(grid_data, res), "Writing of nh grid failed for AGG"

    def test_add_polygon_agg(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_polygons_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        polygons = {
            "REYKJAVIK_ATC_A": (
                (-20.0, 73.0),
                (0.0, 73.0),
                (0.0, 61.0),
                (-30.0, 61.0),
                (-39.0, 63.5),
                (-20, 70),
            ),
            "REYKJAVIK_ATC_B": (
                (-39, 63.5),
                (-55 + 4 / 6.0, 63.5),
                (-57 + 45 / 60.0, 65),
                (-76, 76),
                (-75, 78),
                (-60, 82),
                (0, 90),
                (30, 82),
                (0, 82),
                (0, 73),
                (-20, 73),
                (-20, 70),
            ),
            "REYKJAVIK_ATC": (
                (0.0, 73.0),
                (0.0, 61.0),
                (-30.0, 61.0),
                (-39, 63.5),
                (-55 + 4 / 6.0, 63.5),
                (-57 + 45 / 60.0, 65),
                (-76, 76),
                (-75, 78),
                (-60, 82),
                (0, 90),
                (30, 82),
                (0, 82),
            ),
            "ICELAND_BOX": ((-25, 62.5), (-25, 67), (-13, 67), (-13, 62.5)),
        }

        cw.add_polygon(img, area_def, polygons["REYKJAVIK_ATC"], outline="red", width=2)
        cw.add_polygon(
            img,
            area_def,
            polygons["ICELAND_BOX"],
            outline="green",
            fill="gray",
            width=2,
        )
        cw.add_coastlines(img, area_def, resolution="l", level=4)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh polygons failed"

    def test_add_points_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_points_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (1024, 1024), (255, 255, 255))
        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, outline="black", resolution="l", level=4)
        cw.add_borders(img, area_def, outline="black", width=3, level=1, resolution="c")

        points_list = [((2.3522, 48.8566), "Paris"), ((0.1278, 51.5074), "London")]
        cw.add_points(
            img,
            area_def,
            points_list=points_list,
            font_file=font_path,
            symbol="circle",
            ptsize=16,
            outline="black",
            width=3,
            fill="red",
            fill_opacity=128,
            box_outline="blue",
            box_linewidth=0.5,
            box_fill="yellow",
            box_opacity=200,
        )

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh points failed"

    def test_add_shapefile_shapes_agg(self):
        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "brazil_shapefiles_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m"
        area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_coastlines(img, area_def, resolution="l", level=4)
        cw.add_shapefile_shapes(
            img,
            area_def,
            os.path.join(LOCAL_DIR, "test_data/shapes/Metareas.shp"),
            outline="red",
            width=2,
        )
        cw.add_shapefile_shape(
            img,
            area_def,
            os.path.join(
                LOCAL_DIR,
                "test_data/shapes/divisao_politica/BR_Regioes.shp",
            ),
            3,
            outline="blue",
        )
        cw.add_shapefile_shape(
            img,
            area_def,
            os.path.join(
                LOCAL_DIR,
                "test_data/shapes/divisao_politica/BR_Regioes.shp",
            ),
            4,
            outline="blue",
            fill="green",
        )

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of Brazil shapefiles failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_config_file_coasts_and_grid(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        overlay_config = os.path.join(LOCAL_DIR, "coasts_and_grid_agg.ini")
        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_nh_cfg_agg.png"))
        grid_data = np.array(grid_img)
        proj_dict = {
            "proj": "laea",
            "lat_0": 90.0,
            "lon_0": 0.0,
            "a": 6371228.0,
            "units": "m",
        }
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
        area_def = AreaDefinition("nh", "nh", "nh", proj_dict, 850, 850, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)
        overlay = cw.add_overlay_from_config(overlay_config, area_def)
        img = Image.new("RGB", (850, 850), (255, 255, 255))
        img.paste(overlay, mask=overlay)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh grid failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_config_file_points_and_borders_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        config_file = os.path.join(LOCAL_DIR, "nh_points_agg.ini")

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_points_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (1024, 1024), (255, 255, 255))

        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_overlay_from_config(config_file, area_def, img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Add points with agg module from a config file failed"

    def test_coastlines_convert_to_rgba_agg(self):
        from pycoast import ContourWriterAGG

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = (proj4_string, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines_to_file(p_coasts_filename, area_def, resolution="l", level=4)

        img = Image.open(p_coasts_filename)
        image_mode = img.mode
        img.close()

        assert image_mode == "RGBA", "Conversion to RGBA failed."

    def test_add_cities_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_cities_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (1024, 1024), (255, 255, 255))
        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)
        cw.add_coastlines(img, area_def, outline="black", resolution="l", level=4)
        cw.add_borders(img, area_def, outline="black", width=3, level=1, resolution="c")

        cities_list = ["Zurich", "Oslo", "Reykjavik", "Fairbanks", "Toronto"]
        cw.add_cities(
            img,
            area_def,
            cities_list=cities_list,
            font_file=font_path,
            font_size=20,
            symbol="square",
            ptsize=16,
            outline="black",
            width=1,
            fill="blue",
            fill_opacity=128,
            box_outline="blue",
            box_linewidth=0.5,
            box_fill="yellow",
            box_opacity=200,
        )

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh cities_agg failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_add_cities_cfg_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        config_file = os.path.join(LOCAL_DIR, "nh_cities_agg.ini")

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_cities_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (1024, 1024), (255, 255, 255))

        proj4_string = "+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m"
        area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_overlay_from_config(config_file, area_def, img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh cities_cfg_agg failed"

    def test_add_cities_from_dict_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "nh_cities_from_dict_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (1024, 1024), (255, 255, 255))

        proj4_string = "+proj=stere +ellps=WGS84 +lat_0=51.5 +lon_0=-0.1"
        area_extent = (-2000000.0, -2000000.0, 2000000.0, 2000000.0)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 1024, 1024, area_extent)

        overlays = {}
        overlays["coasts"] = {"level": 4, "resolution": "l", "outline": "black"}
        overlays["borders"] = {"level": 1, "outline": "black", "resolution": "c"}
        cities_list1 = [
            "Berlin",
            "Paris",
            "London",
            "Dublin",
            "Madrid",
            "Reykjavik",
            "Oslo",
            "Rome",
        ]
        cities_list2 = ["Freiburg", "Montelimar", "Huesca", "Marseille"]
        cities_list3 = ["Belp", "Bad Schwalbach", "Edinburgh", "Hilversum"]
        cities_type1 = {
            "cities_list": cities_list1,
            "font": font_path,
            "font_size": 26,
            "symbol": "circle",
            "ptsize": 24,
            "outline": "black",
            "fill": "red",
            "box_outline": "black",
            "box_fill": "blue",
            "box_opacity": 50,
            "box_linewidth": 2.0,
        }
        cities_type2 = {
            "cities_list": cities_list2,
            "font": font_path,
            "font_size": 24,
            "symbol": "pentagon",
            "ptsize": 24,
            "outline": "red",
            "fill": "blue",
            "box_outline": "black",
            "box_fill": "yellow",
            "box_opacity": 50,
            "box_linewidth": 2.0,
        }
        cities_type3 = {
            "cities_list": cities_list3,
            "font": font_path,
            "font_size": 22,
            "symbol": "star5",
            "ptsize": 35,
            "outline": "green",
            "fill": "yellow",
            "box_outline": "black",
            "box_fill": "red",
            "box_opacity": 50,
            "box_linewidth": 2.0,
        }

        overlays["cities"] = [cities_type1, cities_type2, cities_type3]

        cw = ContourWriterAGG(gshhs_root_dir)

        img = cw.add_overlay_from_dict(overlays, area_def, background=img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing of nh_cities_from_dict_agg failed"

    def test_add_shapefiles_from_dict_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "two_shapefiles_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m"
        area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 425, 425, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        overlays = {}
        overlays["coasts"] = {"level": 4, "resolution": "l"}
        overlays["shapefiles"] = [
            {
                "filename": os.path.join(LOCAL_DIR, "test_data/shapes/Metareas.shp"),
                "outline": "magenta",
                "width": 2.5,
            },
            {
                "filename": os.path.join(
                    LOCAL_DIR,
                    "test_data/shapes/divisao_politica/BR_Regioes.shp",
                ),
                "outline": "red",
                "fill": "yellow",
                "fill_opacity": 50,
            },
        ]

        img = cw.add_overlay_from_dict(overlays, area_def, background=img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing two shapefiles from dict agg failed"

    @pytest.mark.usefixtures("cd_test_dir")
    def test_add_one_shapefile_from_cfg_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        config_file = os.path.join(LOCAL_DIR, "nh_one_shapefile.ini")

        grid_img = Image.open(os.path.join(LOCAL_DIR, "one_shapefile_from_cfg_agg.png"))
        grid_data = np.array(grid_img)

        img = Image.new("RGB", (425, 425))
        proj4_string = "+proj=merc +lon_0=-60 +lat_ts=-30.0 +a=6371228.0 +units=m"
        area_extent = (-2000000.0, -5000000.0, 5000000.0, 2000000.0)
        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 425, 425, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        cw.add_overlay_from_config(config_file, area_def, img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing one shapefile from cfg agg failed"

    def test_add_grid_from_dict_agg(self):
        from pyresample.geometry import AreaDefinition

        from pycoast import ContourWriterAGG

        grid_img = Image.open(os.path.join(LOCAL_DIR, "grid_from_dict_agg.png"))

        grid_data = np.array(grid_img)

        img = Image.new("RGB", (800, 800))
        proj4_string = "+proj=stere +ellps=WGS84 +lon_0=-4.532 +lat_0=54.228"
        area_extent = (-600000.0, -600000.0, 600000.0, 600000.0)

        area_def = AreaDefinition("nh", "nh", "nh", proj4_string, 800, 800, area_extent)

        cw = ContourWriterAGG(gshhs_root_dir)

        font = aggdraw.Font(
            "yellow",
            font_path,
            opacity=255,
            size=40,
        )

        overlays = {}
        overlays["coasts"] = {"width": 3.0, "level": 4, "resolution": "l"}
        overlays["grid"] = {
            "major_lonlat": (5, 5),
            "minor_lonlat": (1, 1),
            "outline": (255, 0, 0),
            "outline_opacity": 127,
            "minor_outline": (0, 0, 255),
            "minor_outline_opacity": 127,
            "width": 10.5,
            "minor_width": 5.5,
            "minor_is_tick": False,
            "write_text": True,
            "lat_placement": "lr",
            "lon_placement": "b",
            "font": font,
            "fill": "red",
        }
        # Fill has no agg effect! Agg Font can be None if and only if write_text is set to False

        img = cw.add_overlay_from_dict(overlays, area_def, background=img)

        res = np.array(img)
        assert fft_metric(grid_data, res), "Writing grid from dict agg failed"


class FakeAreaDef:
    """A fake area definition object."""

    def __init__(self, proj4_string, area_extent, x_size, y_size):
        self.proj_str = self.proj_dict = self.crs = proj4_string
        self.area_extent = area_extent
        self.width = x_size
        self.height = y_size
        self.area_id = "fakearea"


class TestFromConfig:
    """Test burning overlays from a config file."""

    def test_foreground(self):
        """Test generating a transparent foreground."""
        from pycoast import ContourWriterPIL

        euro_img = Image.open(os.path.join(LOCAL_DIR, "contours_europe_alpha.png"))
        euro_data = np.array(euro_img)

        # img = Image.new('RGB', (640, 480))
        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = FakeAreaDef(proj4_string, area_extent, 640, 480)
        cw = ContourWriterPIL(gshhs_root_dir)
        config_file = os.path.join(LOCAL_DIR, "test_data", "test_config.ini")
        img = cw.add_overlay_from_config(config_file, area_def)

        res = np.array(img)
        assert fft_metric(euro_data, res), "Writing of contours failed"

        overlays = {
            "coasts": {"level": [1, 2, 3, 4], "resolution": "l"},
            "borders": {"outline": (255, 0, 0), "resolution": "c"},
            "rivers": {"outline": "blue", "resolution": "c", "level": 5},
        }

        img = cw.add_overlay_from_dict(overlays, area_def)
        res = np.array(img)
        assert fft_metric(euro_data, res), "Writing of contours failed"

    def test_cache_generation_reuse(self, tmpdir):
        """Test generating a transparent foreground and cache it."""
        from pycoast import ContourWriterPIL

        euro_img = Image.open(os.path.join(LOCAL_DIR, "contours_europe_alpha.png"))
        euro_data = np.array(euro_img)

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = FakeAreaDef(proj4_string, area_extent, 640, 480)
        cw = ContourWriterPIL(gshhs_root_dir)

        overlays = {
            "cache": {"file": os.path.join(tmpdir, "pycoast_cache")},
            "coasts": {"level": 4, "resolution": "l"},
            "borders": {"outline": (255, 0, 0), "resolution": "c"},
            "rivers": {"outline": "blue", "resolution": "c", "level": 5},
        }

        # Create the original cache file
        img = cw.add_overlay_from_dict(overlays, area_def)
        res = np.array(img)
        cache_glob = glob(os.path.join(tmpdir, "pycoast_cache_*.png"))
        assert len(cache_glob) == 1
        cache_filename = cache_glob[0]
        assert fft_metric(euro_data, res), "Writing of contours failed"
        assert os.path.isfile(cache_filename)
        mtime = os.path.getmtime(cache_filename)

        # Reuse the generated cache file
        img = cw.add_overlay_from_dict(overlays, area_def)
        res = np.array(img)
        assert fft_metric(euro_data, res), "Writing of contours failed"
        assert os.path.isfile(cache_filename)
        assert os.path.getmtime(cache_filename) == mtime

        # Regenerate cache file
        current_time = time.time()
        cw.add_overlay_from_dict(overlays, area_def, current_time)
        mtime = os.path.getmtime(cache_filename)
        assert mtime > current_time
        assert fft_metric(euro_data, res), "Writing of contours failed"

        cw.add_overlay_from_dict(overlays, area_def, current_time)
        assert os.path.getmtime(cache_filename) == mtime
        assert fft_metric(euro_data, res), "Writing of contours failed"
        overlays["cache"]["regenerate"] = True
        cw.add_overlay_from_dict(overlays, area_def)

        assert os.path.getmtime(cache_filename) != mtime
        assert fft_metric(euro_data, res), "Writing of contours failed"

        overlays.pop("cache")
        overlays["grid"] = {
            "outline": (255, 255, 255),
            "outline_opacity": 175,
            "minor_outline": (200, 200, 200),
            "minor_outline_opacity": 127,
            "width": 1.0,
            "minor_width": 0.5,
            "minor_is_tick": True,
            "write_text": True,
            "lat_placement": "lr",
            "lon_placement": "b",
        }
        cw.add_overlay_from_dict(overlays, area_def)
        os.remove(cache_filename)

    def test_caching_with_param_changes(self, tmpdir):
        """Testing caching when changing parameters."""
        from pycoast import ContourWriterPIL

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = FakeAreaDef(proj4_string, area_extent, 640, 480)
        cw = ContourWriterPIL(gshhs_root_dir)

        font = ImageFont.truetype(font_path)
        overlays = {
            "cache": {"file": os.path.join(tmpdir, "pycoast_cache")},
            "grid": {"font": font},
        }

        # Create the original cache file
        cw.add_overlay_from_dict(overlays, area_def)
        cache_glob = glob(os.path.join(tmpdir, "pycoast_cache_*.png"))
        assert len(cache_glob) == 1
        cache_filename = cache_glob[0]
        assert os.path.isfile(cache_filename)
        mtime = os.path.getmtime(cache_filename)

        # Reuse the generated cache file
        cw.add_overlay_from_dict(overlays, area_def)
        cache_glob = glob(os.path.join(tmpdir, "pycoast_cache_*.png"))
        assert len(cache_glob) == 1
        assert os.path.isfile(cache_filename)
        assert os.path.getmtime(cache_filename) == mtime

        # Remove the font option, should produce the same result
        # font is not considered when caching
        del overlays["grid"]["font"]
        cw.add_overlay_from_dict(overlays, area_def)
        cache_glob = glob(os.path.join(tmpdir, "pycoast_cache_*.png"))
        assert len(cache_glob) == 1
        assert os.path.isfile(cache_filename)
        assert os.path.getmtime(cache_filename) == mtime

        # Changing a parameter should create a new cache file
        overlays = {
            "cache": {"file": os.path.join(tmpdir, "pycoast_cache")},
            "grid": {"width": 2.0},
        }
        cw.add_overlay_from_dict(overlays, area_def)
        cache_glob = glob(os.path.join(tmpdir, "pycoast_cache_*.png"))
        assert len(cache_glob) == 2
        assert os.path.isfile(cache_filename)
        new_cache_filename = cache_glob[0] if cache_glob[0] != cache_filename else cache_glob[1]
        # original cache file should be unchanged
        assert os.path.getmtime(cache_filename) == mtime
        # new cache file should be...new
        assert os.path.getmtime(new_cache_filename) != mtime

    @pytest.mark.parametrize("background_mode", ["RGB", "RGBA"])
    @pytest.mark.parametrize("include_background_pattern", [False, True])
    @pytest.mark.parametrize("upper_right_opacity", [32, 127, 255])
    def test_cache_nocache_consistency(
        self, tmp_path, include_background_pattern, background_mode, upper_right_opacity
    ):
        """Test that an image generated with an image looks the same when using a cached foreground."""
        from pycoast import ContourWriterAGG

        proj4_string = "+proj=longlat +ellps=WGS84"
        area_extent = (-10.0, -10.0, 10.0, 10.0)
        area_def = FakeAreaDef(proj4_string, area_extent, 200, 200)
        cw = ContourWriterAGG(gshhs_root_dir)

        # create test shapefiles
        test_shape_filename1 = tmp_path / "test_shapes1"
        _create_polygon_shapefile(test_shape_filename1, [[[-10.0, 10.0], [-5.0, 10.0], [-5.0, 5.0], [-10.0, 5.0]]])
        test_shape_filename2 = tmp_path / "test_shapes2"
        _create_polygon_shapefile(test_shape_filename2, [[[5.0, 10.0], [10.0, 10.0], [10.0, 5.0], [5.0, 5.0]]])

        overlays = {
            "cache": {"file": os.path.join(tmp_path, "pycoast_cache")},
            "shapefiles": [
                {
                    "filename": str(test_shape_filename1),
                    "fill": (0, 255, 0),
                    "outline": (255, 255, 255),
                    "fill_opacity": 255,
                },
                {
                    "filename": str(test_shape_filename2),
                    "fill": (0, 255, 0),
                    "outline": (255, 255, 0),
                    "fill_opacity": upper_right_opacity,
                },
            ],
        }

        # Create the original cache file
        background_img1 = _create_background_image(include_background_pattern, background_mode)
        cached_image1 = cw.add_overlay_from_dict(overlays, area_def, background=background_img1)

        # Reuse the generated cache file
        background_img2 = _create_background_image(include_background_pattern, background_mode)
        cached_image2 = cw.add_overlay_from_dict(overlays, area_def, background=background_img2)

        # Create without cache
        overlays.pop("cache")
        background_img3 = _create_background_image(include_background_pattern, background_mode)
        cw.add_overlay_from_dict(overlays, area_def, background=background_img3)

        # Manually (no dict, no cache)
        background_img4 = _create_background_image(include_background_pattern, background_mode)
        for shape_params in overlays["shapefiles"]:
            cw.add_shapefile_shapes(background_img4, area_def, **shape_params)

        # two transparent overlay images
        np.testing.assert_allclose(np.array(cached_image1), np.array(cached_image2), atol=0)
        # cached overlay applied to background image should always be equal
        np.testing.assert_allclose(np.array(background_img1), np.array(background_img2), atol=0)
        # cached overlay applied to background image and regenerated overlay on background image should be equal
        # but due to floating point differences they are off by a little bit
        np.testing.assert_allclose(
            np.array(background_img1, dtype=np.float32), np.array(background_img3, dtype=np.float32), atol=1
        )
        # no cache version and manual version should be the same
        np.testing.assert_allclose(np.array(background_img3), np.array(background_img4), atol=0)

    def test_get_resolution(self):
        """Get the automagical resolution computation."""
        from pycoast import get_resolution_from_area

        proj4_string = "+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84"
        area_extent = (-3363403.31, -2291879.85, 2630596.69, 2203620.1)
        area_def = FakeAreaDef(proj4_string, area_extent, 640, 480)
        assert get_resolution_from_area(area_def) == "l"
        area_def = FakeAreaDef(proj4_string, area_extent, 6400, 4800)
        assert get_resolution_from_area(area_def) == "h"


def _create_polygon_shapefile(fn: pathlib.Path, polygon_coords: list) -> None:
    with shapefile.Writer(fn) as test_shapefile:
        test_shapefile.field("name", "C")
        test_shapefile.poly(polygon_coords)
        test_shapefile.record("test")


def _create_background_image(add_pattern: bool, background_mode: str) -> Image:
    num_bands = len(background_mode)
    img_data = np.zeros((200, 200, num_bands), dtype=np.uint8)
    if background_mode[-1] == "A":
        img_data[..., -1] = 255
    if add_pattern:
        img_data[6:30, 6:30, 0] = 127
        img_data[6:30, -30:-6, 0] = 127
    return Image.fromarray(img_data, mode=background_mode)
