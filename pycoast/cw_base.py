#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011-2022 PyCoast Developers
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
"""Base class for contour writers."""
from __future__ import annotations

import ast
import configparser
import hashlib
import json
import logging
import math
import os
from typing import Generator

import numpy as np
import pyproj
import shapefile
from PIL import Image

try:
    from pyresample import AreaDefinition
except ImportError:
    AreaDefinition = None

logger = logging.getLogger(__name__)


def get_resolution_from_area(area_def):
    """Get the best resolution for an area definition."""
    x_size = area_def.width
    y_size = area_def.height
    prj = Proj(area_def.crs if hasattr(area_def, "crs") else area_def.proj_str)
    if prj.is_latlong():
        x_ll, y_ll = prj(area_def.area_extent[0], area_def.area_extent[1])
        x_ur, y_ur = prj(area_def.area_extent[2], area_def.area_extent[3])
        x_resolution = (x_ur - x_ll) / x_size
        y_resolution = (y_ur - y_ll) / y_size
    else:
        x_resolution = (area_def.area_extent[2] - area_def.area_extent[0]) / x_size
        y_resolution = (area_def.area_extent[3] - area_def.area_extent[1]) / y_size
    res = min(x_resolution, y_resolution)

    if res > 25000:
        return "c"
    elif res > 5000:
        return "l"
    elif res > 1000:
        return "i"
    elif res > 200:
        return "h"
    else:
        return "f"


class _CoordConverter:
    """Convert coordinates from one space to in-bound image pixel column and row.

    Convert the coordinate (x,y) in the coordinates
    reference system ('lonlat' or 'image') into an image
    x,y coordinate.
    Uses the area_def methods if coord_ref is 'lonlat'.
    Raises ValueError if pixel coordinates are outside the image bounds
    defined by area_def.width and area_def.height.

    """

    def __init__(self, coord_ref: str, area_def: AreaDefinition):
        self._area_def = self._check_area_def(area_def)
        convert_methods = {
            "lonlat": self._lonlat_to_pixels,
            "image": self._image_to_pixels,
        }
        if coord_ref not in convert_methods:
            pretty_coord_refs = [
                f"'{cr_name}'" for cr_name in sorted(convert_methods.keys())
            ]
            raise ValueError(f"'coord_ref' must be one of {pretty_coord_refs}.")
        self._convert_method = convert_methods[coord_ref]

    def _check_area_def(self, area_def):
        if AreaDefinition is None:
            raise ImportError(
                "Missing required 'pyresample' module, please "
                "install it with 'pip install pyresample' or "
                "'conda install pyresample'."
            )
        if not isinstance(area_def, AreaDefinition):
            raise ValueError("'area_def' must be an instance of AreaDefinition")
        return area_def

    def __call__(self, x, y):
        return self._convert_method(x, y)

    def _lonlat_to_pixels(self, x, y):
        return self._area_def.get_array_indices_from_lonlat(x, y)

    def _image_to_pixels(self, x, y):
        area_def = self._area_def
        x, y = (int(x), int(y))
        if x < 0:
            x += area_def.width
        if y < 0:
            y += area_def.height
        if x < 0 or y < 0 or x >= area_def.width or y >= area_def.height:
            raise ValueError(
                "Image pixel coords out of image bounds "
                f"(width={area_def.width}, height={area_def.height})."
            )
        return x, y


def hash_dict(dict_to_hash: dict) -> str:
    """Hash dict object by serializing with json."""
    dhash = hashlib.sha256()
    encoded = json.dumps(dict_to_hash, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


class Proj(pyproj.Proj):
    """Wrapper around pyproj to add in 'is_latlong'."""

    def is_latlong(self):
        if hasattr(self, "crs"):
            return self.crs.is_geographic
        # pyproj<2.0
        return super(Proj, self).is_latlong()


class ContourWriterBase(object):
    """Base class for contourwriters. Do not instantiate.

    :Parameters:
        db_root_path : str
            Path to root dir of GSHHS and WDBII shapefiles

    """

    _draw_module = "FIXME"
    # This is a flag to make _add_grid aware of which draw.text
    # subroutine, from PIL, aggdraw or cairo is being used
    # (unfortunately they are not fully compatible).

    def __init__(self, db_root_path=None):
        if db_root_path is None:
            self.db_root_path = os.environ.get("GSHHS_DATA_ROOT")
        else:
            self.db_root_path = db_root_path

    @property
    def is_agg(self) -> bool:
        """Get if we are using the 'agg' backend."""
        return self._draw_module == "AGG"

    def _draw_text(self, draw, position, txt, font, align="cc", **kwargs):
        """Draw text with agg module."""
        txt_width, txt_height = draw.textsize(txt, font)
        x_pos, y_pos = position
        ax, ay = align.lower()
        if ax == "r":
            x_pos = x_pos - txt_width
        elif ax == "c":
            x_pos = x_pos - txt_width / 2

        if ay == "b":
            y_pos = y_pos - txt_height
        elif ay == "c":
            y_pos = y_pos - txt_width / 2

        self._engine_text_draw(draw, x_pos, y_pos, txt, font, **kwargs)

    def _engine_text_draw(self, draw, pos, txt, font, **kwargs):
        raise NotImplementedError("Text drawing undefined for render engine")

    def _draw_grid_labels(self, draw, xys, linetype, txt, font, **kwargs):
        """Draw text with default PIL module."""
        if font is None:
            # NOTE: Default font does not use font size in PIL writer
            font = self._get_font(kwargs.get("fill", "black"), font, 12)
        placement_def = kwargs[linetype].lower()
        for xy in xys:
            # note xy[0] is xy coordinate pair,
            # xy[1] is required alignment e.g. 'tl','lr','lc','cc'...
            ax, ay = xy[1].lower()
            if ax in placement_def or ay in placement_def:
                self._draw_text(draw, xy[0], txt, font, align=xy[1], **kwargs)

    def _find_line_intercepts(self, xys, size, margins):
        """Find intercepts of poly-line xys with image boundaries offset by margins.

        Returns an array of coordinates.

        """
        x_size, y_size = size

        def is_in_box(x_y, extents):
            x, y = x_y
            xmin, xmax, ymin, ymax = extents
            return xmin < x < xmax and ymin < y < ymax

        def crossing(x1, x2, lim):
            return (x1 < lim) != (x2 < lim)

        # set box limits
        xlim1 = margins[0]
        ylim1 = margins[1]
        xlim2 = x_size - margins[0]
        ylim2 = y_size - margins[1]

        # only consider crossing within a box a little bigger than grid boundary
        search_box = (-10, x_size + 10, -10, y_size + 10)

        # loop through line steps and detect crossings
        intercepts = []
        align_left = "LC"
        align_right = "RC"
        align_top = "CT"
        align_bottom = "CB"
        prev_xy = xys[0]
        for xy in xys[1:]:
            if not is_in_box(xy, search_box):
                prev_xy = xy
                continue
            # crossing LHS
            if crossing(prev_xy[0], xy[0], xlim1):
                x = xlim1
                y = xy[1]
                intercepts.append(((x, y), align_left))
            # crossing RHS
            elif crossing(prev_xy[0], xy[0], xlim2):
                x = xlim2
                y = xy[1]
                intercepts.append(((x, y), align_right))
            # crossing Top
            elif crossing(prev_xy[1], xy[1], ylim1):
                x = xy[0]
                y = ylim1
                intercepts.append(((x, y), align_top))
            # crossing Bottom
            elif crossing(prev_xy[1], xy[1], ylim2):
                x = xy[0]  # - txt_width/2
                y = ylim2  # - txt_height
                intercepts.append(((x, y), align_bottom))
            prev_xy = xy

        return intercepts

    def _add_grid(
        self,
        image,
        area_def,
        Dlon,
        Dlat,
        dlon,
        dlat,
        font=None,
        write_text=True,
        **kwargs,
    ):
        """Add a lat lon grid to image."""
        try:
            proj_def = area_def.crs if hasattr(area_def, "crs") else area_def.proj_dict
            area_extent = area_def.area_extent
        except AttributeError:
            proj_def = area_def[0]
            area_extent = area_def[1]

        draw = self._get_canvas(image)

        # use kwargs for major lines ... but reform for minor lines:
        minor_line_kwargs = kwargs.copy()
        minor_line_kwargs["outline"] = kwargs["minor_outline"]
        if self.is_agg:
            minor_line_kwargs["outline_opacity"] = kwargs["minor_outline_opacity"]
            minor_line_kwargs["width"] = kwargs["minor_width"]

        # text margins (at sides of image frame)
        y_text_margin = 4
        x_text_margin = 4

        # Area and projection info
        x_size, y_size = image.size
        prj = Proj(proj_def)

        x_offset = 0
        y_offset = 0

        # Calculate min and max lons and lats of interest
        lon_min, lon_max, lat_min, lat_max = _get_lon_lat_bounding_box(
            area_extent, x_size, y_size, prj
        )

        # Handle dateline crossing
        if lon_max < lon_min:
            lon_max = 360 + lon_max

        # Draw lonlat grid lines ...
        # create adjustment of line lengths to avoid cluttered pole lines
        if lat_max == 90.0:
            shorten_max_lat = Dlat
        else:
            shorten_max_lat = 0.0

        if lat_min == -90.0:
            increase_min_lat = Dlat
        else:
            increase_min_lat = 0.0

        # major lon lines
        round_lon_min = lon_min - (lon_min % Dlon)
        maj_lons = np.arange(round_lon_min, lon_max, Dlon)
        maj_lons[maj_lons > 180] = maj_lons[maj_lons > 180] - 360

        # minor lon lines (ticks)
        min_lons = np.arange(round_lon_min, lon_max, dlon)
        min_lons[min_lons > 180] = min_lons[min_lons > 180] - 360

        # Get min_lons not in maj_lons
        min_lons = np.lib.arraysetops.setdiff1d(min_lons, maj_lons)

        # lats along major lon lines
        lin_lats = np.arange(
            lat_min + increase_min_lat,
            lat_max - shorten_max_lat,
            float(lat_max - lat_min) / y_size,
        )
        # lin_lats in rather high definition so that it can be used to
        # position text labels near edges of image...

        # perhaps better to find the actual length of line in pixels...
        round_lat_min = lat_min - (lat_min % Dlat)

        # major lat lines
        maj_lats = np.arange(round_lat_min + increase_min_lat, lat_max, Dlat)

        # minor lon lines (ticks)
        min_lats = np.arange(
            round_lat_min + increase_min_lat, lat_max - shorten_max_lat, dlat
        )

        # Get min_lats not in maj_lats
        min_lats = np.lib.arraysetops.setdiff1d(min_lats, maj_lats)

        # lons along major lat lines (extended slightly to avoid missing the end)
        lin_lons = np.linspace(lon_min, lon_max + Dlon / 5.0, max(x_size, y_size) // 5)

        # MINOR LINES ######
        if not kwargs["minor_is_tick"]:
            minor_lat_lines = [[(x, lat) for x in lin_lons] for lat in min_lats]
            minor_lon_lines = [[(lon, x) for x in lin_lats] for lon in min_lons]
            minor_lines = minor_lat_lines + minor_lon_lines
            self._draw_minor_grid_lines(
                minor_lines,
                draw,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset,
                y_offset,
                minor_line_kwargs,
            )

        # MAJOR LINES AND MINOR TICKS ######
        # major lon lines and tick marks:
        for lon in maj_lons:
            # Draw 'minor' tick lines dlat separation along the lon
            if kwargs["minor_is_tick"]:
                tick_lons = np.linspace(lon - Dlon / 20.0, lon + Dlon / 20.0, 5)
                minor_tick_lines = [[(x, lat) for x in tick_lons] for lat in min_lats]
                self._draw_minor_grid_lines(
                    minor_tick_lines,
                    draw,
                    area_extent,
                    x_size,
                    y_size,
                    prj,
                    x_offset,
                    y_offset,
                    minor_line_kwargs,
                )

            # Draw 'major' lines
            lonlats = [(lon, x) for x in lin_lats]
            index_arrays = self._grid_line_index_array_generator(
                [lonlats], area_extent, x_size, y_size, prj, x_offset, y_offset
            )
            for index_array in index_arrays:
                self._draw_line(draw, index_array.flatten().tolist(), **kwargs)

            # add lon text markings at each end of longitude line
            if write_text:
                txt = self._grid_lon_label(lon)
                xys = self._find_line_intercepts(
                    index_array, image.size, (x_text_margin, y_text_margin)
                )

                self._draw_grid_labels(draw, xys, "lon_placement", txt, font, **kwargs)

        # major lat lines and tick marks:
        for lat in maj_lats:
            # Draw 'minor' tick dlon separation along the lat
            if kwargs["minor_is_tick"]:
                tick_lats = np.linspace(lat - Dlat / 20.0, lat + Dlat / 20.0, 5)
                minor_tick_lines = [[(lon, x) for x in tick_lats] for lon in min_lons]
                self._draw_minor_grid_lines(
                    minor_tick_lines,
                    draw,
                    area_extent,
                    x_size,
                    y_size,
                    prj,
                    x_offset,
                    y_offset,
                    minor_line_kwargs,
                )

            # Draw 'major' lines
            lonlats = [(x, lat) for x in lin_lons]
            index_arrays = self._grid_line_index_array_generator(
                [lonlats], area_extent, x_size, y_size, prj, x_offset, y_offset
            )
            for index_array in index_arrays:
                self._draw_line(draw, index_array.flatten().tolist(), **kwargs)

            # add lat text markings at each end of parallels ...
            if write_text:
                txt = self._grid_lat_label(lat)
                xys = self._find_line_intercepts(
                    index_array, image.size, (x_text_margin, y_text_margin)
                )
                self._draw_grid_labels(draw, xys, "lat_placement", txt, font, **kwargs)

        # Draw cross on poles ...
        if lat_max == 90.0:
            crosslats = np.arange(
                90.0 - Dlat / 2.0, 90.0, float(lat_max - lat_min) / y_size
            )
            cross_lines = [
                [(lon, x) for x in crosslats] for lon in (0.0, 90.0, 180.0, -90.0)
            ]
            self._draw_minor_grid_lines(
                cross_lines,
                draw,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset,
                y_offset,
                kwargs,
            )

        if lat_min == -90.0:
            crosslats = np.arange(
                -90.0, -90.0 + Dlat / 2.0, float(lat_max - lat_min) / y_size
            )
            cross_lines = [
                [(lon, x) for x in crosslats] for lon in (0.0, 90.0, 180.0, -90.0)
            ]
            self._draw_minor_grid_lines(
                cross_lines,
                draw,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset,
                y_offset,
                kwargs,
            )
        self._finalize(draw)

    def _draw_minor_grid_lines(
        self,
        minor_lines,
        draw,
        area_extent,
        x_size,
        y_size,
        prj,
        x_offset,
        y_offset,
        kwargs,
    ):
        for minor_line_lonlats in minor_lines:
            index_arrays, _ = _get_pixel_index(
                minor_line_lonlats,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset=x_offset,
                y_offset=y_offset,
            )
            if not index_arrays:
                continue
            for index_array in index_arrays:
                self._draw_line(draw, index_array.flatten().tolist(), **kwargs)

    def _grid_line_index_array_generator(
        self, grid_lines, area_extent, x_size, y_size, prj, x_offset, y_offset
    ):
        for grid_line_lonlats in grid_lines:
            index_arrays, is_reduced = _get_pixel_index(
                grid_line_lonlats,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset=x_offset,
                y_offset=y_offset,
            )
            # Skip empty datasets
            if not index_arrays:
                continue
            yield from index_arrays

    def _grid_lon_label(self, lon):
        # FIXME: Use f-strings or just pass the direction
        if lon > 0.0:
            txt = "%.2dE" % (lon)
        else:
            txt = "%.2dW" % (-lon)
        return txt

    def _grid_lat_label(self, lat):
        if lat >= 0.0:
            txt = "%.2dN" % (lat)
        else:
            txt = "%.2dS" % (-lat)
        return txt

    def _find_bounding_box(self, xys):
        lons = [x for (x, y) in xys]
        lats = [y for (x, y) in xys]
        return [min(lons), min(lats), max(lons), max(lats)]

    def _add_shapefile_shapes(
        self, image, area_def, filename, feature_type=None, **kwargs
    ):
        """Draw all shapes (polygon/poly-lines) from a shape file onto a PIL Image."""
        sf = shapefile.Reader(filename)
        return self.add_shapes(
            image, area_def, sf.shapes(), feature_type=feature_type, **kwargs
        )

    def _add_shapefile_shape(
        self, image, area_def, filename, shape_id, feature_type=None, **kwargs
    ):
        """Draw a single shape (polygon/poly-line) definition.

        Accesses single shape using shape_id from a custom shape file.

        """
        sf = shapefile.Reader(filename)
        shape = sf.shape(shape_id)
        return self.add_shapes(
            image, area_def, [shape], feature_type=feature_type, **kwargs
        )

    def _add_line(self, image, area_def, lonlats, **kwargs):
        """Draw a custom polyline.

        Lon and lat coordinates given by the list lonlat.

        """
        # create dummpy shapelike object
        shape = type("", (), {})()
        shape.points = lonlats
        shape.parts = [0]
        shape.bbox = self._find_bounding_box(lonlats)
        self.add_shapes(image, area_def, [shape], feature_type="line", **kwargs)

    def _add_polygon(self, image, area_def, lonlats, **kwargs):
        """Draw a custom polygon.

        Lon and lat coordinates given by the list lonlat.

        """
        # create dummpy shapelike object
        shape = type("", (), {})()
        shape.points = lonlats
        shape.parts = [0]
        shape.bbox = self._find_bounding_box(lonlats)
        self.add_shapes(image, area_def, [shape], feature_type="polygon", **kwargs)

    def add_shapes(
        self,
        image,
        area_def,
        shapes,
        feature_type=None,
        x_offset=0,
        y_offset=0,
        **kwargs,
    ):
        """Draw shape objects to PIL image.

        :Parameters:
            image : Image
                PIL Image to draw shapes on
            area_def : (proj_str, area_extent) or AreaDefinition
                Geolocation information for the provided image
            shapes: iterable
                Series of shape objects from pyshp. Can also be a series
                of 2-element tuples where the first element is the shape
                object and the second is a dictionary of additional drawing
                parameters for this shape.
            feature_type: str
                'polygon' or 'line' or None for what to draw shapes as.
                Default is to draw the shape with the type in the shapefile.
            kwargs:
                Extra drawing keyword arguments for all shapes

        .. versionchanged: 1.2.0

            Interface changed to have `shapes` before `feature_type` to allow
            `feature_type` to be optional and default to `None`.

        """
        try:
            proj_def = area_def.crs if hasattr(area_def, "crs") else area_def.proj_dict
            area_extent = area_def.area_extent
        except AttributeError:
            proj_def = area_def[0]
            area_extent = area_def[1]

        draw = self._get_canvas(image)

        # Area and projection info
        x_size, y_size = image.size
        prj = Proj(proj_def)

        # Calculate min and max lons and lats of interest
        lon_min, lon_max, lat_min, lat_max = _get_lon_lat_bounding_box(
            area_extent, x_size, y_size, prj
        )

        # Iterate through shapes
        for shape in shapes:
            if isinstance(shape, (list, tuple)):
                new_kwargs = kwargs.copy()
                if shape[1]:
                    new_kwargs.update(shape[1])
                shape = shape[0]
            else:
                new_kwargs = kwargs

            if self._polygon_is_irrelevant(lon_min, lon_max, lat_min, lat_max, shape):
                continue

            self._add_shape(
                shape,
                feature_type,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset,
                y_offset,
                draw,
                new_kwargs,
            )
        self._finalize(draw)

    @staticmethod
    def _polygon_is_irrelevant(lon_min, lon_max, lat_min, lat_max, shape):
        # Check if polygon is possibly relevant
        s_lon_ll, s_lat_ll, s_lon_ur, s_lat_ur = shape.bbox
        if lon_min <= lon_max:
            # Area_extent west or east of dateline
            shape_is_outside_lon = lon_max < s_lon_ll or lon_min > s_lon_ur
        else:
            # Area_extent spans over dateline
            shape_is_outside_lon = lon_max < s_lon_ll and lon_min > s_lon_ur
        shape_is_outside_lat = lat_max < s_lat_ll or lat_min > s_lat_ur
        return shape_is_outside_lon or shape_is_outside_lat

    def _add_shape(
        self,
        shape,
        feature_type,
        area_extent,
        x_size,
        y_size,
        prj,
        x_offset,
        y_offset,
        draw,
        new_kwargs,
    ):
        ftype = self._feature_type_for_shape(shape, feature_type)

        # iterate over shape parts (some shapes split into parts)
        # dummy shape part object
        parts = list(shape.parts) + [len(shape.points)]
        for i in range(len(parts) - 1):
            # Get pixel index coordinates of shape
            points = shape.points[parts[i] : parts[i + 1]]
            index_arrays, is_reduced = _get_pixel_index(
                points,
                area_extent,
                x_size,
                y_size,
                prj,
                x_offset=x_offset,
                y_offset=y_offset,
            )

            # Skip empty datasets
            if len(index_arrays) == 0:
                return

            # Make PIL draw the polygon or line
            for index_array in index_arrays:
                if ftype == "polygon" and not is_reduced:
                    # Draw polygon if dataset has not been reduced
                    self._draw_polygon(
                        draw, index_array.flatten().tolist(), **new_kwargs
                    )
                elif ftype == "line" or is_reduced:
                    # Draw line
                    self._draw_line(draw, index_array.flatten().tolist(), **new_kwargs)
                else:
                    raise ValueError("Unknown contour type: %s" % ftype)

    @staticmethod
    def _feature_type_for_shape(shape, feature_type):
        if feature_type is not None:
            return feature_type.lower()
        if shape.shapeType == shapefile.POLYLINE:
            ftype = "line"
        elif shape.shapeType == shapefile.POLYGON:
            ftype = "polygon"
        else:
            raise ValueError("Unsupported shape type: " + str(shape.shapeType))
        return ftype

    def _add_feature(
        self,
        image,
        area_def,
        feature_type,
        db_name,
        tag=None,
        zero_pad=False,
        resolution="c",
        level=1,
        x_offset=0,
        y_offset=0,
        db_root_path=None,
        **kwargs,
    ):
        """Add a contour feature to image."""
        shape_generator = self._iterate_db(
            db_name, tag, resolution, level, zero_pad, db_root_path=db_root_path
        )

        return self.add_shapes(
            image,
            area_def,
            shape_generator,
            feature_type=feature_type,
            x_offset=x_offset,
            y_offset=y_offset,
            **kwargs,
        )

    def _iterate_db(self, db_name, tag, resolution, level, zero_pad, db_root_path=None):
        """Iterate through datasets."""
        if db_root_path is None:
            db_root_path = self.db_root_path
        if db_root_path is None:
            raise ValueError("'db_root_path' must be specified to use this method")
        levels = range(1, level + 1) if not isinstance(level, list) else level
        format_string, format_params = self._get_db_shapefile_format_and_params(
            db_name, resolution, tag, zero_pad
        )
        shapefile_root_dir = os.path.join(db_root_path, f"{db_name}_shp", resolution)
        for i in levels:
            level_format_params = format_params + (i,)
            shapefilename = os.path.join(
                shapefile_root_dir, format_string % level_format_params
            )
            try:
                s = shapefile.Reader(shapefilename)
                shapes = s.shapes()
            except AttributeError:
                raise ValueError("Could not find shapefile %s" % shapefilename)

            yield from shapes

    @staticmethod
    def _get_db_shapefile_format_and_params(db_name, resolution, tag, zero_pad):
        format_string = "%s_%s_"
        format_params = (db_name, resolution)
        if tag is not None:
            format_string += "%s_"
            format_params = (db_name, tag, resolution)

        if zero_pad:
            format_string += "L%02i.shp"
        else:
            format_string += "L%s.shp"
        return format_string, format_params

    def _finalize(self, draw):
        """Do any need finalization of the drawing."""
        pass

    def _config_to_dict(self, config_file):
        """Convert a config file to a dict."""
        config = configparser.ConfigParser()
        try:
            with open(config_file, "r"):
                logger.info("Overlays config file %s found", str(config_file))
            config.read(config_file)
        except IOError:
            logger.error("Overlays config file %s does not exist!", str(config_file))
            raise
        except configparser.NoSectionError:
            logger.error("Error in %s", str(config_file))
            raise

        SECTIONS = [
            "cache",
            "coasts",
            "rivers",
            "borders",
            "shapefiles",
            "grid",
            "cities",
            "points",
        ]
        overlays = {}
        for section in config.sections():
            if section in SECTIONS:
                overlays[section] = {}
                for option in config.options(section):
                    val = config.get(section, option)
                    try:
                        overlays[section][option] = ast.literal_eval(val)
                    except ValueError:
                        overlays[section][option] = val
        return overlays

    def add_overlay_from_dict(
        self, overlays, area_def, cache_epoch=None, background=None
    ):
        """Create and return a transparent image adding all the overlays contained in the `overlays` dict.

        Optionally caches overlay results for faster rendering of images with
        the same provided AreaDefinition and parameters. Cached results are
        identified by hashing the AreaDefinition and the overlays dictionary.

        .. warning::

            Font objects are ignored in parameter hashing as they can't be easily hashed.
            Therefore font changes will not trigger a new rendering for the cache.

        :Parameters:
            overlays : dict
                overlays configuration
            area_def : object
                Area Definition of the creating image
            cache_epoch: seconds since epoch
                The latest time allowed for cache the cache file. If the cache
                file is older than this (mtime), the cache should be
                regenerated. Defaults to 0 meaning always reuse the cached
                file if it exists. Requires "cache" to be configured in the
                provided dictionary (see below).
            background: pillow image instance
                The image on which to write the overlays on. If it's None (default),
                a new image is created, otherwise the provide background is
                used and changed *in place*.


            The keys in `overlays` that will be taken into account are:
            cache, coasts, rivers, borders, shapefiles, grid, cities, points

            For all of them except `cache`, the items are the same as the
            corresponding functions in pycoast, so refer to the docstrings of
            these functions (add_coastlines, add_rivers, add_borders,
            add_shapefile_shapes, add_grid, add_cities, add_points).
            For cache, two parameters are configurable:

            - `file`: specify the directory and the prefix
                  of the file to save the caches decoration to (for example
                  /var/run/black_coasts_red_borders)
            - `regenerate`: True or False (default) to force the overwriting
                  of an already cached file.

        """
        overlay_helper = _OverlaysFromDict(
            self, overlays, area_def, cache_epoch, background
        )
        return overlay_helper.apply_overlays()

    def add_overlay_from_config(self, config_file, area_def, background=None):
        """Create and return a transparent image adding all the overlays contained in a configuration file.

        :Parameters:
            config_file : str
                Configuration file name
            area_def : object
                Area Definition of the creating image

        """
        overlays = self._config_to_dict(config_file)
        return self.add_overlay_from_dict(
            overlays, area_def, os.path.getmtime(config_file), background
        )

    def draw_star(self, draw, symbol, x, y, ptsize, **kwargs):
        # 5 <= n <= 8, symbol = string in ['star8', 'star7', 'star6', 'star5']
        n = int(symbol[4])
        alpha2 = math.pi / n
        # r1 = outer radius (defaults to 0.5 * ptsize), r1 > r2 = inner radius
        r1 = 0.5 * ptsize
        r2 = r1 / (math.cos(alpha2) + math.sin(alpha2) * math.tan(2.0 * alpha2))
        xy = []
        alpha = 0.0
        # Walk from star top ray CW around the symbol
        for i in range(2 * n):
            if (i % 2) == 0:
                xy.append(x + r1 * math.sin(alpha))
                xy.append(y - r1 * math.cos(alpha))
            else:
                xy.append(x + r2 * math.sin(alpha))
                xy.append(y - r2 * math.cos(alpha))
            alpha += alpha2
        self._draw_polygon(draw, xy, **kwargs)

    def draw_hexagon(self, draw, x, y, ptsize, **kwargs):
        xy = [
            x + 0.25 * ptsize,
            y - 0.4330127 * ptsize,
            x + 0.50 * ptsize,
            y,
            x + 0.25 * ptsize,
            y + 0.4330127 * ptsize,
            x - 0.25 * ptsize,
            y + 0.4330127 * ptsize,
            x - 0.50 * ptsize,
            y,
            x - 0.25 * ptsize,
            y - 0.4330127 * ptsize,
        ]
        self._draw_polygon(draw, xy, **kwargs)

    def draw_pentagon(self, draw, x, y, ptsize, **kwargs):
        xy = [
            x,
            y - 0.5 * ptsize,
            x + 0.4755283 * ptsize,
            y - 0.1545085 * ptsize,
            x + 0.2938926 * ptsize,
            y + 0.4045085 * ptsize,
            x - 0.2938926 * ptsize,
            y + 0.4045085 * ptsize,
            x - 0.4755283 * ptsize,
            y - 0.1545085 * ptsize,
        ]
        self._draw_polygon(draw, xy, **kwargs)

    def draw_triangle(self, draw, x, y, ptsize, **kwargs):
        xy = [
            x,
            y - 0.5 * ptsize,
            x + 0.4330127 * ptsize,
            y + 0.25 * ptsize,
            x - 0.4330127 * ptsize,
            y + 0.25 * ptsize,
        ]
        self._draw_polygon(draw, xy, **kwargs)

    def add_cities(
        self,
        image,
        area_def,
        cities_list,
        font_file,
        font_size=12,
        symbol="circle",
        ptsize=6,
        outline="black",
        fill="white",
        db_root_path=None,
        **kwargs,
    ):
        """Add cities (symbol and UTF-8 names as description) to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            area_def : object
                Area Definition of the provided image
            cities_list : list of city names ['City1', 'City2', City3, ..., 'CityN']
              | a list of UTF-8 or ASCII strings. If either of these strings is found
              | in file db_root_path/CITIES/cities.red, longitude and latitude is read
              | and the city is added like a point with its UTF-8 name as description
              | e.g. cities_list = ['Zurich', 'Oslo'] will add cities 'Zürich', 'Oslo'.
              | Check the README_PyCoast.txt in archive cities2022.zip for more details.
            font_file : str
                Path to font file
            font_size : int
                Size of font
            symbol : string
                type of symbol, one of the elelments from the list
                ['circle', 'hexagon', 'pentagon', 'square', 'triangle',
                'star8', 'star7', 'star6', 'star5', 'asterisk']
            ptsize : int
                Size of the point.
            outline : str or (R, G, B), optional
                Line color of the symbol
            fill : str or (R, G, B), optional
                Filling color of the symbol

        :Optional keyword arguments:
            width : float
                Line width of the symbol
            outline_opacity : int, optional {0; 255}
                Opacity of the line of the symbol.
            fill_opacity : int, optional {0; 255}
                Opacity of the filling of the symbol
            box_outline : str or (R, G, B), optional
                Line color of the textbox borders.
            box_linewidth : float
                Line width of the the borders of the textbox
            box_fill : str or (R, G, B), optional
                Filling color of the background of the textbox.
            box_opacity : int, optional {0; 255}
                Opacity of the background filling of the textbox.
        """
        if db_root_path is None:
            db_root_path = self.db_root_path
        if db_root_path is None:
            raise ValueError("'db_root_path' must be specified to use this method")

        try:
            from pyresample.geometry import AreaDefinition
        except ImportError:
            raise ImportError(
                "Missing required 'pyresample' module, please install it."
            )

        if not isinstance(area_def, AreaDefinition):
            raise ValueError(
                "Expected 'area_def' is an instance of AreaDefinition object"
            )

        draw = self._get_canvas(image)

        # cities.red is a reduced version of the files avalable at http://download.geonames.org
        # Fields: 0=name (UTF-8), 1=asciiname, 2=longitude [°E], 3=latitude [°N], 4=countrycode
        cities_filename = os.path.join(
            db_root_path, os.path.join("CITIES", "cities.txt")
        )
        cities_parser = GeoNamesCitiesParser(cities_filename)
        for city_name, lon, lat in cities_parser.iter_cities_names_lon_lat(cities_list):
            try:
                x, y = area_def.get_array_indices_from_lonlat(lon, lat)
            except ValueError:
                logger.info(
                    "City %s is out of the area, it will not be added to the image.",
                    city_name + " " + str((lon, lat)),
                )
                continue
            # add symbol
            if ptsize != 0:
                half_ptsize = int(round(ptsize / 2.0))
                dot_box = [
                    x - half_ptsize,
                    y - half_ptsize,
                    x + half_ptsize,
                    y + half_ptsize,
                ]

                width = kwargs.get("width", 1.0)
                outline_opacity = kwargs.get("outline_opacity", 255)
                fill_opacity = kwargs.get("fill_opacity", 255)
                self._draw_point_element(
                    draw,
                    symbol,
                    dot_box,
                    x,
                    y,
                    width,
                    ptsize,
                    outline,
                    outline_opacity,
                    fill,
                    fill_opacity,
                )
                text_position = [x + ptsize, y]
            else:
                text_position = [x, y]

            font = self._get_font(outline, font_file, font_size)
            new_kwargs = kwargs.copy()
            box_outline = new_kwargs.pop("box_outline", "white")
            box_opacity = new_kwargs.pop("box_opacity", 0)

            # add text_box
            self._draw_text_box(
                draw,
                text_position,
                city_name,
                font,
                outline,
                box_outline,
                box_opacity,
                **new_kwargs,
            )
            logger.info("%s added", city_name + " " + str((lon, lat)))
        self._finalize(draw)

    def add_points(
        self,
        image,
        area_def,
        points_list,
        font_file,
        font_size=12,
        symbol="circle",
        ptsize=6,
        outline="black",
        fill="white",
        coord_ref="lonlat",
        **kwargs,
    ):
        """Add a symbol and/or text at the point(s) of interest to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            area_def : object
                Area Definition of the provided image
            points_list : list [((x, y), desc)]
              | a list of points defined with (x, y) in float and a desc in string
              | [((x1, y1), desc1), ((x2, y2), desc2)]
              | See coord_ref (below) for the meaning of x, y.
              | x : float
              |    longitude or pixel x of a point
              | y : float
              |    latitude or pixel y of a point
              | desc : str
              |    description of a point
            font_file : str
                Path to font file
            font_size : int
                Size of font
            symbol : string
                type of symbol, one of the elelments from the list
                ['circle', 'hexagon', 'pentagon', 'square', 'triangle',
                'star8', 'star7', 'star6', 'star5, 'asterisk']
            ptsize : int
                Size of the point (should be zero if symbol:None).
            outline : str or (R, G, B), optional
                Line color of the symbol
            fill : str or (R, G, B), optional
                Filling color of the symbol

        :Optional keyword arguments:
            coord_ref : string
                The interpretation of x,y in points_list:
                'lonlat' (the default: x is degrees E, y is degrees N),
                or 'image' (x is pixels right, y is pixels down).
                If image coordinates are negative they are interpreted
                relative to the end of the dimension like standard Python
                indexing.
            width : float
                Line width of the symbol
            outline_opacity : int, optional {0; 255}
                Opacity of the line of the symbol.
            fill_opacity : int, optional {0; 255}
                Opacity of the filling of the symbol
            box_outline : str or (R, G, B), optional
                Line color of the textbox borders.
            box_linewidth : float
                Line width of the the borders of the textbox
            box_fill : str or (R, G, B), optional
                Filling color of the background of the textbox.
            box_opacity : int, optional {0; 255}
                Opacity of the background filling of the textbox.
        """
        coord_converter = _CoordConverter(coord_ref, area_def)
        draw = self._get_canvas(image)

        # Iterate through points list
        for point in points_list:
            (x, y), desc = point
            try:
                x, y = coord_converter(x, y)
            except ValueError:
                logger.info(
                    f"Point ({x}, {y}) is out of the image area, it will not be added to the image."
                )
                continue
            if ptsize != 0:
                half_ptsize = int(round(ptsize / 2.0))
                dot_box = [
                    x - half_ptsize,
                    y - half_ptsize,
                    x + half_ptsize,
                    y + half_ptsize,
                ]

                width = kwargs.get("width", 1.0)
                outline_opacity = kwargs.get("outline_opacity", 255)
                fill_opacity = kwargs.get("fill_opacity", 255)
                self._draw_point_element(
                    draw,
                    symbol,
                    dot_box,
                    x,
                    y,
                    width,
                    ptsize,
                    outline,
                    outline_opacity,
                    fill,
                    fill_opacity,
                )
            elif desc is None:
                logger.error(
                    "'ptsize' is 0 and 'desc' is None, nothing will be added to the image."
                )

            if desc is not None:
                text_position = [
                    x + ptsize,
                    y,
                ]  # draw the text box next to the point
                font = self._get_font(outline, font_file, font_size)

                new_kwargs = kwargs.copy()

                box_outline = new_kwargs.pop("box_outline", "white")
                box_opacity = new_kwargs.pop("box_opacity", 0)

                # add text_box
                self._draw_text_box(
                    draw,
                    text_position,
                    desc,
                    font,
                    outline,
                    box_outline,
                    box_opacity,
                    **new_kwargs,
                )

            logger.debug("Point %s has been added to the image", str((x, y)))

        self._finalize(draw)

    def _draw_point_element(
        self,
        draw,
        symbol,
        dot_box,
        x,
        y,
        width,
        ptsize,
        outline,
        outline_opacity,
        fill,
        fill_opacity,
    ):
        if symbol == "circle":
            # a 'circle' or a 'dot' i.e. circle with fill
            self._draw_ellipse(
                draw,
                dot_box,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        # All regular polygons are drawn horizontally based
        elif symbol == "hexagon":
            self.draw_hexagon(
                draw,
                x,
                y,
                ptsize,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        elif symbol == "pentagon":
            self.draw_pentagon(
                draw,
                x,
                y,
                ptsize,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        elif symbol == "square":
            self._draw_rectangle(
                draw,
                dot_box,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        elif symbol == "triangle":
            self.draw_triangle(
                draw,
                x,
                y,
                ptsize,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        # All stars are drawn with one vertical ray on top
        elif symbol in ["star8", "star7", "star6", "star5"]:
            self.draw_star(
                draw,
                symbol,
                x,
                y,
                ptsize,
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
                fill=fill,
                fill_opacity=fill_opacity,
            )
        elif symbol == "asterisk":  # an '*' sign
            self._draw_asterisk(
                draw,
                ptsize,
                (x, y),
                outline=outline,
                width=width,
                outline_opacity=outline_opacity,
            )
        elif symbol:
            raise ValueError("Unsupported symbol type: " + str(symbol))


def _get_lon_lat_bounding_box(area_extent, x_size, y_size, prj):
    """Get extreme lon and lat values."""
    bbox_lons, bbox_lats = _get_bounding_box_lonlat_sides(
        area_extent, x_size, y_size, prj
    )
    lons_s1, lons_s2, lons_s3, lons_s4 = bbox_lons
    lats_s1, lats_s2, lats_s3, lats_s4 = bbox_lats
    angle_sum = _get_angle_sum(lons_s1, lons_s2, lons_s3, lons_s4)

    if round(angle_sum) == -360:
        # Covers NP
        lat_min = min(lats_s1.min(), lats_s2.min(), lats_s3.min(), lats_s4.min())
        lat_max = 90
        lon_min = -180
        lon_max = 180
    elif round(angle_sum) == 360:
        # Covers SP
        lat_min = -90
        lat_max = max(lats_s1.max(), lats_s2.max(), lats_s3.max(), lats_s4.max())
        lon_min = -180
        lon_max = 180
    elif round(angle_sum) == 0:
        # Covers no poles
        if (
            np.sign(lons_s1[0]) * np.sign(lons_s1[-1]) == -1
            and lons_s1.min() * lons_s1.max() < -25000
        ):
            # End points of left side on different side of dateline
            lon_min = lons_s1[lons_s1 > 0].min()
        else:
            lon_min = lons_s1.min()

        if (
            np.sign(lons_s3[0]) * np.sign(lons_s3[-1]) == -1
            and lons_s3.min() * lons_s3.max() < -25000
        ):
            # End points of right side on different side of dateline
            lon_max = lons_s3[lons_s3 < 0].max()
        else:
            lon_max = lons_s3.max()

        lat_min = lats_s4.min()
        lat_max = lats_s2.max()
    else:
        # Pretty strange
        lat_min = -90
        lat_max = 90
        lon_min = -180
        lon_max = 180

    # Catch inf/1e30 or other invalid values
    if not (-180 <= lon_min <= 180):
        lon_min = -180
    if not (-180 <= lon_max <= 180):
        lon_max = 180
    if not (-90 <= lat_min <= 90):
        lat_min = -90
    if not (-90 <= lat_max <= 90):
        lat_max = 90

    return lon_min, lon_max, lat_min, lat_max


def _get_bounding_box_lonlat_sides(area_extent, x_size, y_size, prj):
    x_ll, y_ll, x_ur, y_ur = area_extent
    x_range = np.linspace(x_ll, x_ur, num=x_size)
    y_range = np.linspace(y_ll, y_ur, num=y_size)

    if prj.is_latlong():
        lons_s1, lats_s1 = x_ll * np.ones(y_range.size), y_range
        lons_s2, lats_s2 = x_range, y_ur * np.ones(x_range.size)
        lons_s3, lats_s3 = x_ur * np.ones(y_range.size), y_range
        lons_s4, lats_s4 = x_range, y_ll * np.ones(x_range.size)
    else:
        lons_s1, lats_s1 = prj(np.ones(y_range.size) * x_ll, y_range, inverse=True)
        lons_s2, lats_s2 = prj(x_range, np.ones(x_range.size) * y_ur, inverse=True)
        lons_s3, lats_s3 = prj(np.ones(y_range.size) * x_ur, y_range, inverse=True)
        lons_s4, lats_s4 = prj(x_range, np.ones(x_range.size) * y_ll, inverse=True)
    return (lons_s1, lons_s2, lons_s3, lons_s4), (lats_s1, lats_s2, lats_s3, lats_s4)


def _get_angle_sum(lons_s1, lons_s2, lons_s3, lons_s4):
    angle_sum = 0
    prev = None
    for lon in np.concatenate((lons_s1, lons_s2, lons_s3[::-1], lons_s4[::-1])):
        if not np.isfinite(lon):
            continue
        if prev is not None:
            delta = lon - prev
            if abs(delta) > 180:
                delta = (abs(delta) - 360) * np.sign(delta)
            angle_sum += delta
        prev = lon
    return angle_sum


def _get_pixel_index(shape, area_extent, x_size, y_size, prj, x_offset=0, y_offset=0):
    """Map coordinates of shape to image coordinates."""
    # Get shape data as array and reproject
    shape_data = np.array(shape.points if hasattr(shape, "points") else shape)
    lons = shape_data[:, 0]
    lats = shape_data[:, 1]

    if prj.is_latlong():
        x_ll, y_ll = prj(area_extent[0], area_extent[1])
        x_ur, y_ur = prj(area_extent[2], area_extent[3])
    else:
        x_ll, y_ll, x_ur, y_ur = area_extent

    x, y = prj(lons, lats)

    # Handle out of bounds
    i = 0
    segments = []
    if (x >= 1e30).any() or (y >= 1e30).any():
        # Split polygon in line segments within projection
        is_reduced = True
        if x[0] >= 1e30 or y[0] >= 1e30:
            in_segment = False
        else:
            in_segment = True

        for j in range(x.size):
            if x[j] >= 1e30 or y[j] >= 1e30:
                if in_segment:
                    segments.append((x[i:j], y[i:j]))
                    in_segment = False
            elif not in_segment:
                in_segment = True
                i = j
        if in_segment:
            segments.append((x[i:], y[i:]))

    else:
        is_reduced = False
        segments = [(x, y)]

    # Convert to pixel index coordinates
    l_x = (x_ur - x_ll) / x_size
    l_y = (y_ur - y_ll) / y_size

    index_arrays = []
    for x, y in segments:
        n_x = ((-x_ll + x) / l_x) + 0.5 + x_offset
        n_y = ((y_ur - y) / l_y) + 0.5 + y_offset

        index_array = np.vstack((n_x, n_y)).T
        index_arrays.append(index_array)

    return index_arrays, is_reduced


class GeoNamesCitiesParser:
    """Helper for parsing citiesN.txt files from GeoNames.org."""

    def __init__(self, cities_filename: str):
        self._cities_file = open(cities_filename, mode="r", encoding="utf-8")

    def iter_cities_names_lon_lat(
        self, cities_list: list[str]
    ) -> Generator[tuple[str, float, float], None, None]:
        for city_row in self._cities_file:
            city_info = city_row.split("\t")
            if not city_info or not (
                city_info[1] in cities_list or city_info[2] in cities_list
            ):
                continue
            city_name, lon, lat = city_info[1], float(city_info[5]), float(city_info[4])
            yield city_name, lon, lat


class _OverlaysFromDict:
    """Helper for drawing overlays from a dictionary of parameters."""

    def __init__(self, contour_writer, overlays, area_def, cache_epoch, background):
        self._cw = contour_writer
        self._overlays = overlays
        self._is_cached = False
        self._cache_filename = None
        self._background = background
        self._area_def = area_def

        foreground = None
        if "cache" in overlays:
            cache_filename, foreground = self._get_cached_filename_and_foreground(
                cache_epoch
            )
            self._cache_filename = cache_filename
            self._is_cached = foreground is not None

        if foreground is None:
            if self._cache_filename is None and self._background is not None:
                foreground = self._background
            else:
                x_size = area_def.width
                y_size = area_def.height
                foreground = Image.new("RGBA", (x_size, y_size), (0, 0, 0, 0))

        self._foreground = foreground

    def _get_cached_filename_and_foreground(self, cache_epoch):
        cache_file = self._generate_cache_filename(
            self._overlays["cache"]["file"],
            self._area_def,
            self._overlays,
        )
        regenerate = self._overlays["cache"].get("regenerate", False)
        foreground = self._apply_cached_image(
            cache_file, cache_epoch, self._background, regenerate=regenerate
        )
        return cache_file, foreground

    @staticmethod
    def _apply_cached_image(cache_file, cache_epoch, background, regenerate=False):
        try:
            config_time = cache_epoch or 0
            cache_time = os.path.getmtime(cache_file)
            # Cache file will be used only if it's newer than config file
            if config_time is not None and config_time < cache_time and not regenerate:
                foreground = Image.open(cache_file)
                logger.info("Using image in cache %s", cache_file)
                if background is not None:
                    background.paste(foreground, mask=foreground.split()[-1])
                return foreground
            logger.info("Regenerating cache file.")
        except OSError:
            logger.info(
                "No overlay image found, new overlay image will be saved in cache."
            )
        return None

    def _write_and_apply_new_cached_image(self):
        try:
            self._foreground.save(self._cache_filename)
        except IOError as e:
            logger.error("Can't save cache: %s", str(e))
        if self._background is not None:
            self._background.paste(self._foreground, mask=self._foreground.split()[-1])

    def _generate_cache_filename(self, cache_prefix, area_def, overlays_dict):
        area_hash = hash(area_def)
        base_dir, file_prefix = os.path.split(cache_prefix)
        params_to_hash = self._prepare_hashable_dict(overlays_dict)
        param_hash = hash_dict(params_to_hash)
        return os.path.join(base_dir, f"{file_prefix}_{area_hash}_{param_hash}.png")

    @staticmethod
    def _prepare_hashable_dict(nonhashable_dict):
        params_to_hash = {}
        # avoid wasteful deep copy by only doing two levels of copying
        for overlay_name, overlay_dict in nonhashable_dict.items():
            if overlay_name == "cache":
                continue
            params_to_hash[overlay_name] = overlay_dict.copy()
        # font objects are not hashable
        for font_cat in ("cities", "points", "grid"):
            if font_cat in params_to_hash:
                params_to_hash[font_cat].pop("font", None)
        return params_to_hash

    def apply_overlays(self):
        if self._is_cached:
            return self._foreground

        overlays = self._overlays
        self._add_coasts_rivers_borders_from_dict(overlays)
        if "shapefiles" in overlays:
            self._add_shapefiles_from_dict(overlays["shapefiles"])
        if "grid" in overlays:
            self._add_grid_from_dict(overlays["grid"])
        if "cities" in overlays:
            self._add_cities_from_dict(overlays["cities"])
        for param_key in ["points", "text"]:
            if param_key not in overlays:
                continue
            self._add_points_from_dict(overlays[param_key])

        if self._cache_filename is not None:
            self._write_and_apply_new_cached_image()
        return self._foreground

    def _add_coasts_rivers_borders_from_dict(self, overlays):
        default_resolution = get_resolution_from_area(self._area_def)
        DEFAULT = {
            "level": 1,
            "outline": "white",
            "width": 1,
            "fill": None,
            "fill_opacity": 255,
            "outline_opacity": 255,
            "x_offset": 0,
            "y_offset": 0,
            "resolution": default_resolution,
        }

        for section, fun in zip(
            ["coasts", "rivers", "borders"],
            [self._cw.add_coastlines, self._cw.add_rivers, self._cw.add_borders],
        ):
            if section not in overlays:
                continue
            params = DEFAULT.copy()
            params.update(overlays[section])

            if section != "coasts":
                params.pop("fill_opacity", None)
                params.pop("fill", None)

            if not self._cw.is_agg:
                for key in ["width", "outline_opacity", "fill_opacity"]:
                    params.pop(key, None)

            fun(self._foreground, self._area_def, **params)
            logger.info("%s added", section.capitalize())

    def _add_shapefiles_from_dict(self, shapefiles):
        # Backward compatibility and config.ini
        if isinstance(shapefiles, dict):
            shapefiles = [shapefiles]

        DEFAULT_FILENAME = None
        DEFAULT_OUTLINE = "white"
        DEFAULT_FILL = None
        for params in shapefiles:
            params = params.copy()  # don't modify the user's dictionary
            params.setdefault("filename", DEFAULT_FILENAME)
            params.setdefault("outline", DEFAULT_OUTLINE)
            params.setdefault("fill", DEFAULT_FILL)
            if not self._cw.is_agg:
                for key in ["width", "outline_opacity", "fill_opacity"]:
                    params.pop(key, None)
            self._cw.add_shapefile_shapes(
                self._foreground,
                self._area_def,
                feature_type=None,
                x_offset=0,
                y_offset=0,
                **params,
            )

    def _add_grid_from_dict(self, grid_dict):
        if "major_lonlat" in grid_dict or "minor_lonlat" in grid_dict:
            Dlonlat = grid_dict.get("major_lonlat", (10.0, 10.0))
            dlonlat = grid_dict.get("minor_lonlat", (2.0, 2.0))
        else:
            Dlonlat = (
                grid_dict.get("lon_major", 10.0),
                grid_dict.get("lat_major", 10.0),
            )
            dlonlat = (
                grid_dict.get("lon_minor", 2.0),
                grid_dict.get("lat_minor", 2.0),
            )
        outline = grid_dict.get("outline", "white")
        write_text = grid_dict.get("write_text", True)
        if isinstance(write_text, str):
            write_text = write_text.lower() in ["true", "yes", "1", "on"]
        font = grid_dict.get("font", None)
        font_size = int(grid_dict.get("font_size", 10))
        fill = grid_dict.get("fill", outline)
        fill_opacity = grid_dict.get("fill_opacity", 255)
        if isinstance(font, str):
            if self._cw.is_agg:
                from aggdraw import Font

                font = Font(fill, font, opacity=fill_opacity, size=font_size)
            else:
                from PIL.ImageFont import truetype

                font = truetype(font, font_size)
        minor_outline = grid_dict.get("minor_outline", "white")
        minor_is_tick = grid_dict.get("minor_is_tick", True)
        if isinstance(minor_is_tick, str):
            minor_is_tick = minor_is_tick.lower() in ["true", "yes", "1"]
        lon_placement = grid_dict.get("lon_placement", "tb")
        lat_placement = grid_dict.get("lat_placement", "lr")

        grid_kwargs = {}
        if self._cw.is_agg:
            width = float(grid_dict.get("width", 1.0))
            minor_width = float(grid_dict.get("minor_width", 0.5))
            outline_opacity = grid_dict.get("outline_opacity", 255)
            minor_outline_opacity = grid_dict.get("minor_outline_opacity", 255)
            grid_kwargs["width"] = width
            grid_kwargs["minor_width"] = minor_width
            grid_kwargs["outline_opacity"] = outline_opacity
            grid_kwargs["minor_outline_opacity"] = minor_outline_opacity

        self._cw.add_grid(
            self._foreground,
            self._area_def,
            Dlonlat,
            dlonlat,
            font=font,
            write_text=write_text,
            fill=fill,
            outline=outline,
            minor_outline=minor_outline,
            minor_is_tick=minor_is_tick,
            lon_placement=lon_placement,
            lat_placement=lat_placement,
            **grid_kwargs,
        )

    def _add_cities_from_dict(self, cities_dict):
        # Backward compatibility and config.ini
        if isinstance(cities_dict, dict):
            cities_dict = [cities_dict]

        DEFAULT_FONTSIZE = 12
        DEFAULT_SYMBOL = "circle"
        DEFAULT_PTSIZE = 6
        DEFAULT_OUTLINE = "black"
        DEFAULT_FILL = "white"

        for params in cities_dict:
            params = params.copy()
            cities_list = params.pop("cities_list")
            font_file = params.pop("font")
            font_size = int(params.pop("font_size", DEFAULT_FONTSIZE))
            symbol = params.pop("symbol", DEFAULT_SYMBOL)
            ptsize = int(params.pop("ptsize", DEFAULT_PTSIZE))
            outline = params.pop("outline", DEFAULT_OUTLINE)
            fill = params.pop("fill", DEFAULT_FILL)

            self._cw.add_cities(
                self._foreground,
                self._area_def,
                cities_list,
                font_file,
                font_size,
                symbol,
                ptsize,
                outline,
                fill,
                **params,
            )

    def _add_points_from_dict(self, points_dict):
        DEFAULT_FONTSIZE = 12
        DEFAULT_SYMBOL = "circle"
        DEFAULT_PTSIZE = 6
        DEFAULT_OUTLINE = "black"
        DEFAULT_FILL = "white"

        params = points_dict.copy()
        points_list = list(params.pop("points_list"))
        font_file = params.pop("font")
        font_size = int(params.pop("font_size", DEFAULT_FONTSIZE))
        symbol = params.pop("symbol", DEFAULT_SYMBOL)
        ptsize = int(params.pop("ptsize", DEFAULT_PTSIZE))
        outline = params.pop("outline", DEFAULT_OUTLINE)
        fill = params.pop("fill", DEFAULT_FILL)

        self._cw.add_points(
            self._foreground,
            self._area_def,
            points_list,
            font_file,
            font_size,
            symbol,
            ptsize,
            outline,
            fill,
            **params,
        )
