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
"""Base class for contour writers."""

import os
import hashlib
import json
import shapefile
import numpy as np
from PIL import Image
import pyproj
import logging
import ast

import configparser

logger = logging.getLogger(__name__)


def get_resolution_from_area(area_def):
    """Get the best resolution for an area definition."""
    x_size = area_def.width
    y_size = area_def.height
    prj = Proj(area_def.crs if hasattr(area_def, 'crs') else area_def.proj_str)
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


def hash_dict(dict_to_hash: dict) -> str:
    """Hash dict object by serializing with json."""
    dhash = hashlib.sha256()
    encoded = json.dumps(dict_to_hash, sort_keys=True).encode()
    dhash.update(encoded)
    return dhash.hexdigest()


class Proj(pyproj.Proj):
    """Wrapper around pyproj to add in 'is_latlong'."""

    def is_latlong(self):
        if hasattr(self, 'crs'):
            return self.crs.is_geographic
        # pyproj<2.0
        return super(Proj, self).is_latlong()


class ContourWriterBase(object):
    """Base class for contourwriters. Do not instantiate.

    :Parameters:
        db_root_path : str
            Path to root dir of GSHHS and WDBII shapefiles

    """

    _draw_module = None
    # This is a flag to make _add_grid aware of which draw.text
    # subroutine, from PIL, aggdraw or cairo is being used
    # (unfortunately they are not fully compatible).

    def __init__(self, db_root_path=None):
        if db_root_path is None:
            self.db_root_path = os.environ.get('GSHHS_DATA_ROOT')
        else:
            self.db_root_path = db_root_path

    def _draw_text(self, draw, position, txt, font, align='cc', **kwargs):
        """Draw text with agg module."""
        txt_width, txt_height = draw.textsize(txt, font)
        x_pos, y_pos = position
        ax, ay = align.lower()
        if ax == 'r':
            x_pos = x_pos - txt_width
        elif ax == 'c':
            x_pos = x_pos - txt_width / 2

        if ay == 'b':
            y_pos = y_pos - txt_height
        elif ay == 'c':
            y_pos = y_pos - txt_width / 2

        self._engine_text_draw(draw, x_pos, y_pos, txt, font, **kwargs)

    def _engine_text_draw(self, draw, pos, txt, font, **kwargs):
        raise NotImplementedError('Text drawing undefined for render engine')

    def _draw_grid_labels(self, draw, xys, linetype, txt, font, **kwargs):
        """Draw text with default PIL module."""
        if font is None:
            # NOTE: Default font does not use font size in PIL writer
            font = self._get_font(kwargs.get('fill', 'black'), font, 12)
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
            if xmin < x < xmax and ymin < y < ymax:
                return True
            else:
                return False

        def crossing(x1, x2, lim):
            if (x1 < lim) != (x2 < lim):
                return True
            else:
                return False

        # set box limits
        xlim1 = margins[0]
        ylim1 = margins[1]
        xlim2 = x_size - margins[0]
        ylim2 = y_size - margins[1]

        # only consider crossing within a box a little bigger than grid
        # boundary
        search_box = (-10, x_size + 10, -10, y_size + 10)

        # loop through line steps and detect crossings
        intercepts = []
        align_left = 'LC'
        align_right = 'RC'
        align_top = 'CT'
        align_bottom = 'CB'
        prev_xy = xys[0]
        for i in range(1, len(xys) - 1):
            xy = xys[i]
            if is_in_box(xy, search_box):
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

    def _add_grid(self, image, area_def,
                  Dlon, Dlat,
                  dlon, dlat,
                  font=None, write_text=True, **kwargs):
        """Add a lat lon grid to image."""
        try:
            proj_def = area_def.crs if hasattr(area_def, 'crs') else area_def.proj_dict
            area_extent = area_def.area_extent
        except AttributeError:
            proj_def = area_def[0]
            area_extent = area_def[1]

        draw = self._get_canvas(image)

        is_agg = self._draw_module == "AGG"

        # use kwargs for major lines ... but reform for minor lines:
        minor_line_kwargs = kwargs.copy()
        minor_line_kwargs['outline'] = kwargs['minor_outline']
        if is_agg:
            minor_line_kwargs['outline_opacity'] = \
                kwargs['minor_outline_opacity']
            minor_line_kwargs['width'] = kwargs['minor_width']

        # text margins (at sides of image frame)
        y_text_margin = 4
        x_text_margin = 4

        # Area and projection info
        x_size, y_size = image.size
        prj = Proj(proj_def)

        x_offset = 0
        y_offset = 0

        # Calculate min and max lons and lats of interest
        lon_min, lon_max, lat_min, lat_max = \
            _get_lon_lat_bounding_box(area_extent, x_size, y_size, prj)

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
        round_lon_min = (lon_min - (lon_min % Dlon))
        maj_lons = np.arange(round_lon_min, lon_max, Dlon)
        maj_lons[maj_lons > 180] = maj_lons[maj_lons > 180] - 360

        # minor lon lines (ticks)
        min_lons = np.arange(round_lon_min, lon_max, dlon)
        min_lons[min_lons > 180] = min_lons[min_lons > 180] - 360

        # Get min_lons not in maj_lons
        min_lons = np.lib.arraysetops.setdiff1d(min_lons, maj_lons)

        # lats along major lon lines
        lin_lats = np.arange(lat_min + increase_min_lat,
                             lat_max - shorten_max_lat,
                             float(lat_max - lat_min) / y_size)
        # lin_lats in rather high definition so that it can be used to
        # position text labels near edges of image...

        # perhaps better to find the actual length of line in pixels...

        round_lat_min = (lat_min - (lat_min % Dlat))

        # major lat lines
        maj_lats = np.arange(round_lat_min + increase_min_lat, lat_max, Dlat)

        # minor lon lines (ticks)
        min_lats = np.arange(round_lat_min + increase_min_lat,
                             lat_max - shorten_max_lat,
                             dlat)

        # Get min_lats not in maj_lats
        min_lats = np.lib.arraysetops.setdiff1d(min_lats, maj_lats)

        # lons along major lat lines (extended slightly to avoid missing the
        # end)
        lin_lons = np.linspace(lon_min, lon_max + Dlon / 5.0, max(x_size, y_size) // 5)

        # MINOR LINES ######
        if not kwargs['minor_is_tick']:
            # minor lat lines
            for lat in min_lats:
                lonlats = [(x, lat) for x in lin_lons]
                index_arrays, is_reduced = _get_pixel_index(lonlats,
                                                            area_extent,
                                                            x_size, y_size,
                                                            prj,
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)
                del is_reduced
                # Skip empty datasets
                if len(index_arrays) == 0:
                    continue
                # make PIL draw the tick line...
                for index_array in index_arrays:
                    self._draw_line(draw,
                                    index_array.flatten().tolist(),
                                    **minor_line_kwargs)
            # minor lon lines
            for lon in min_lons:
                lonlats = [(lon, x) for x in lin_lats]
                index_arrays, is_reduced = _get_pixel_index(lonlats,
                                                            area_extent,
                                                            x_size, y_size,
                                                            prj,
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)
                # Skip empty datasets
                if len(index_arrays) == 0:
                    continue
                # make PIL draw the tick line...
                for index_array in index_arrays:
                    self._draw_line(draw,
                                    index_array.flatten().tolist(),
                                    **minor_line_kwargs)

        # MAJOR LINES AND MINOR TICKS ######
        # major lon lines and tick marks:
        for lon in maj_lons:
            # Draw 'minor' tick lines dlat separation along the lon
            if kwargs['minor_is_tick']:
                tick_lons = np.linspace(lon - Dlon / 20.0,
                                        lon + Dlon / 20.0,
                                        5)

                for lat in min_lats:
                    lonlats = [(x, lat) for x in tick_lons]
                    index_arrays, is_reduced = \
                        _get_pixel_index(lonlats,
                                         area_extent,
                                         x_size, y_size,
                                         prj,
                                         x_offset=x_offset,
                                         y_offset=y_offset)
                    # Skip empty datasets
                    if len(index_arrays) == 0:
                        continue
                    # make PIL draw the tick line...
                    for index_array in index_arrays:
                        self._draw_line(draw,
                                        index_array.flatten().tolist(),
                                        **minor_line_kwargs)

            # Draw 'major' lines
            lonlats = [(lon, x) for x in lin_lats]
            index_arrays, is_reduced = _get_pixel_index(lonlats, area_extent,
                                                        x_size, y_size,
                                                        prj,
                                                        x_offset=x_offset,
                                                        y_offset=y_offset)
            # Skip empty datasets
            if len(index_arrays) == 0:
                continue

            # make PIL draw the lines...
            for index_array in index_arrays:
                self._draw_line(draw,
                                index_array.flatten().tolist(),
                                **kwargs)

            # add lon text markings at each end of longitude line
            if write_text:
                if lon > 0.0:
                    txt = "%.2dE" % (lon)
                else:
                    txt = "%.2dW" % (-lon)
                xys = self._find_line_intercepts(index_array, image.size,
                                                 (x_text_margin,
                                                  y_text_margin))

                self._draw_grid_labels(draw, xys, 'lon_placement',
                                       txt, font, **kwargs)

        # major lat lines and tick marks:
        for lat in maj_lats:
            # Draw 'minor' tick dlon separation along the lat
            if kwargs['minor_is_tick']:
                tick_lats = np.linspace(lat - Dlat / 20.0,
                                        lat + Dlat / 20.0,
                                        5)
                for lon in min_lons:
                    lonlats = [(lon, x) for x in tick_lats]
                    index_arrays, is_reduced = \
                        _get_pixel_index(lonlats, area_extent,
                                         x_size, y_size,
                                         prj,
                                         x_offset=x_offset,
                                         y_offset=y_offset)
                    # Skip empty datasets
                    if len(index_arrays) == 0:
                        continue
                    # make PIL draw the tick line...
                    for index_array in index_arrays:
                        self._draw_line(draw,
                                        index_array.flatten().tolist(),
                                        **minor_line_kwargs)

            # Draw 'major' lines
            lonlats = [(x, lat) for x in lin_lons]
            index_arrays, is_reduced = _get_pixel_index(lonlats, area_extent,
                                                        x_size, y_size,
                                                        prj,
                                                        x_offset=x_offset,
                                                        y_offset=y_offset)
            # Skip empty datasets
            if len(index_arrays) == 0:
                continue

            # make PIL draw the lines...
            for index_array in index_arrays:
                self._draw_line(draw, index_array.flatten().tolist(), **kwargs)

            # add lat text markings at each end of parallels ...
            if write_text:
                if lat >= 0.0:
                    txt = "%.2dN" % (lat)
                else:
                    txt = "%.2dS" % (-lat)
                xys = self._find_line_intercepts(index_array, image.size,
                                                 (x_text_margin,
                                                  y_text_margin))
                self._draw_grid_labels(draw, xys, 'lat_placement',
                                       txt, font, **kwargs)

        # Draw cross on poles ...
        if lat_max == 90.0:
            crosslats = np.arange(90.0 - Dlat / 2.0, 90.0,
                                  float(lat_max - lat_min) / y_size)
            for lon in (0.0, 90.0, 180.0, -90.0):
                lonlats = [(lon, x) for x in crosslats]
                index_arrays, is_reduced = _get_pixel_index(lonlats,
                                                            area_extent,
                                                            x_size, y_size,
                                                            prj,
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)
                # Skip empty datasets
                if len(index_arrays) == 0:
                    continue

                # make PIL draw the lines...
                for index_array in index_arrays:
                    self._draw_line(draw,
                                    index_array.flatten().tolist(),
                                    **kwargs)
        if lat_min == -90.0:
            crosslats = np.arange(-90.0, -90.0 + Dlat / 2.0,
                                  float(lat_max - lat_min) / y_size)
            for lon in (0.0, 90.0, 180.0, -90.0):
                lonlats = [(lon, x) for x in crosslats]
                index_arrays, is_reduced = _get_pixel_index(lonlats,
                                                            area_extent,
                                                            x_size, y_size,
                                                            prj,
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)
                # Skip empty datasets
                if len(index_arrays) == 0:
                    continue

                # make PIL draw the lines...
                for index_array in index_arrays:
                    self._draw_line(draw,
                                    index_array.flatten().tolist(),
                                    **kwargs)
        self._finalize(draw)

    def _find_bounding_box(self, xys):
        lons = [x for (x, y) in xys]
        lats = [y for (x, y) in xys]
        return [min(lons), min(lats), max(lons), max(lats)]

    def _add_shapefile_shapes(self, image, area_def, filename,
                              feature_type=None, **kwargs):
        """Draw all shapes (polygon/poly-lines) from a shape file onto a PIL Image."""
        sf = shapefile.Reader(filename)
        return self.add_shapes(image, area_def, sf.shapes(), feature_type=feature_type, **kwargs)

    def _add_shapefile_shape(self, image, area_def, filename, shape_id,
                             feature_type=None, **kwargs):
        """Draw a single shape (polygon/poly-line) definition.

        Accesses single shape using shape_id from a custom shape file.

        """
        sf = shapefile.Reader(filename)
        shape = sf.shape(shape_id)
        return self.add_shapes(image, area_def, [shape], feature_type=feature_type, **kwargs)

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

    def add_shapes(self, image, area_def, shapes, feature_type=None, x_offset=0, y_offset=0, **kwargs):
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
            proj_def = area_def.crs if hasattr(area_def, 'crs') else area_def.proj_dict
            area_extent = area_def.area_extent
        except AttributeError:
            proj_def = area_def[0]
            area_extent = area_def[1]

        draw = self._get_canvas(image)

        # Area and projection info
        x_size, y_size = image.size
        prj = Proj(proj_def)

        # Calculate min and max lons and lats of interest
        lon_min, lon_max, lat_min, lat_max = _get_lon_lat_bounding_box(area_extent,
                                                                       x_size,
                                                                       y_size,
                                                                       prj)

        # Iterate through shapes
        for shape in shapes:
            if isinstance(shape, (list, tuple)):
                new_kwargs = kwargs.copy()
                if shape[1]:
                    new_kwargs.update(shape[1])
                shape = shape[0]
            else:
                new_kwargs = kwargs

            if feature_type is None:
                if shape.shapeType == shapefile.POLYLINE:
                    ftype = "line"
                elif shape.shapeType == shapefile.POLYGON:
                    ftype = "polygon"
                else:
                    raise ValueError("Unsupported shape type: " + str(shape.shapeType))
            else:
                ftype = feature_type.lower()

            # Check if polygon is possibly relevant
            s_lon_ll, s_lat_ll, s_lon_ur, s_lat_ur = shape.bbox
            if lon_min <= lon_max:
                # Area_extent west or east of dateline
                shape_is_outside_lon = lon_max < s_lon_ll or lon_min > s_lon_ur
            else:
                # Area_extent spans over dateline
                shape_is_outside_lon = lon_max < s_lon_ll and lon_min > s_lon_ur
            shape_is_outside_lat = lat_max < s_lat_ll or lat_min > s_lat_ur
            if shape_is_outside_lon or shape_is_outside_lat:
                # Polygon is irrelevant
                continue

            # iterate over shape parts (some shapes split into parts)
            # dummy shape part object
            parts = list(shape.parts) + [len(shape.points)]
            for i in range(len(parts) - 1):
                # Get pixel index coordinates of shape
                points = shape.points[parts[i]:parts[i + 1]]
                index_arrays, is_reduced = _get_pixel_index(points,
                                                            area_extent,
                                                            x_size, y_size,
                                                            prj,
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)

                # Skip empty datasets
                if len(index_arrays) == 0:
                    continue

                # Make PIL draw the polygon or line
                for index_array in index_arrays:
                    if ftype == 'polygon' and not is_reduced:
                        # Draw polygon if dataset has not been reduced
                        self._draw_polygon(draw, index_array.flatten().tolist(), **new_kwargs)
                    elif ftype == 'line' or is_reduced:
                        # Draw line
                        self._draw_line(draw, index_array.flatten().tolist(), **new_kwargs)
                    else:
                        raise ValueError('Unknown contour type: %s' % ftype)

        self._finalize(draw)

    def _add_feature(self, image, area_def, feature_type,
                     db_name, tag=None, zero_pad=False, resolution='c',
                     level=1, x_offset=0, y_offset=0, db_root_path=None,
                     **kwargs):
        """Add a contour feature to image."""
        shape_generator = self._iterate_db(
            db_name, tag, resolution, level, zero_pad,
            db_root_path=db_root_path
        )

        return self.add_shapes(image, area_def, shape_generator, feature_type=feature_type,
                               x_offset=x_offset, y_offset=y_offset, **kwargs)

    def _iterate_db(self, db_name, tag, resolution, level, zero_pad, db_root_path=None):
        """Iterate through datasets."""
        if db_root_path is None:
            db_root_path = self.db_root_path
        if db_root_path is None:
            raise ValueError("'db_root_path' must be specified to use this method")

        format_string = '%s_%s_'
        if tag is not None:
            format_string += '%s_'

        if zero_pad:
            format_string += 'L%02i.shp'
        else:
            format_string += 'L%s.shp'

        if not isinstance(level, list):
            level = range(1, level + 1)

        for i in level:

            # One shapefile per level
            if tag is None:
                shapefilename = \
                    os.path.join(db_root_path, '%s_shp' % db_name,
                                 resolution, format_string %
                                 (db_name, resolution, i))
            else:
                shapefilename = \
                    os.path.join(db_root_path, '%s_shp' % db_name,
                                 resolution, format_string %
                                 (db_name, tag, resolution, i))
            try:
                s = shapefile.Reader(shapefilename)
                shapes = s.shapes()
            except AttributeError:
                raise ValueError('Could not find shapefile %s'
                                 % shapefilename)

            for shape in shapes:
                yield shape

    def _finalize(self, draw):
        """Do any need finalization of the drawing."""
        pass

    def _config_to_dict(self, config_file):
        """Convert a config file to a dict."""
        config = configparser.ConfigParser()
        try:
            with open(config_file, 'r'):
                logger.info("Overlays config file %s found", str(config_file))
            config.read(config_file)
        except IOError:
            logger.error("Overlays config file %s does not exist!",
                         str(config_file))
            raise
        except configparser.NoSectionError:
            logger.error("Error in %s", str(config_file))
            raise

        SECTIONS = ['cache', 'coasts', 'rivers', 'borders', 'cities', 'points', 'grid']
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

    def add_overlay_from_dict(self, overlays, area_def, cache_epoch=None, background=None):
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
            cache, coasts, rivers, borders, cities, points, grid

            For all of them except `cache`, the items are the same as the
            corresponding functions in pycoast, so refer to the docstrings of
            these functions (add_coastlines, add_rivers, add_borders,
            add_grid, add_cities, add_points). For cache, two parameters are
            configurable:

            - `file`: specify the directory and the prefix
                  of the file to save the caches decoration to (for example
                  /var/run/black_coasts_red_borders)
            - `regenerate`: True or False (default) to force the overwriting
                  of an already cached file.

        """
        # Cache management
        cache_file = None
        if 'cache' in overlays:
            cache_file = self._generate_cache_filename(
                overlays['cache']['file'],
                area_def,
                overlays,
            )
            regenerate = overlays['cache'].get('regenerate', False)
            foreground = self._apply_cached_image(cache_file, cache_epoch, background, regenerate=regenerate)
            if foreground is not None:
                return foreground

        x_size = area_def.width
        y_size = area_def.height
        if cache_file is None and background is not None:
            foreground = background
        else:
            foreground = Image.new('RGBA', (x_size, y_size), (0, 0, 0, 0))

        is_agg = self._draw_module == "AGG"

        # Coasts, rivers, borders
        default_resolution = get_resolution_from_area(area_def)

        DEFAULT = {'level': 1,
                   'outline': 'white',
                   'width': 1,
                   'fill': None,
                   'fill_opacity': 255,
                   'outline_opacity': 255,
                   'x_offset': 0,
                   'y_offset': 0,
                   'resolution': default_resolution}

        for section, fun in zip(['coasts', 'rivers', 'borders'],
                                [self.add_coastlines, self.add_rivers, self.add_borders]):
            if section not in overlays:
                continue
            params = DEFAULT.copy()
            params.update(overlays[section])

            if section != "coasts":
                params.pop('fill_opacity', None)
                params.pop('fill', None)

            if not is_agg:
                for key in ['width', 'outline_opacity', 'fill_opacity']:
                    params.pop(key, None)

            fun(foreground, area_def, **params)
            logger.info("%s added", section.capitalize())

        # Cities management
        if 'cities' in overlays:
            DEFAULT_FONT_SIZE = 12
            DEFAULT_OUTLINE = "yellow"

            citylist = [s.lstrip()
                        for s in overlays['cities']['list'].split(',')]
            font_file = overlays['cities']['font']
            font_size = int(overlays['cities'].get('font_size',
                                                   DEFAULT_FONT_SIZE))
            outline = overlays['cities'].get('outline', DEFAULT_OUTLINE)
            pt_size = int(overlays['cities'].get('pt_size', None))
            box_outline = overlays['cities'].get('box_outline', None)
            box_opacity = int(overlays['cities'].get('box_opacity', 255))

            self.add_cities(foreground, area_def, citylist, font_file,
                            font_size, pt_size, outline, box_outline,
                            box_opacity)
        # Points management
        if 'points' in overlays:
            DEFAULT_FONTSIZE = 12
            DEFAULT_SYMBOL = 'circle'
            DEFAULT_PTSIZE = 6
            DEFAULT_OUTLINE = 'white'
            DEFAULT_FILL = None

            params = overlays['points'].copy()

            points_list = list(params.pop('points_list'))
            font_file = params.pop('font')
            font_size = int(params.pop('font_size', DEFAULT_FONTSIZE))

            symbol = params.pop('symbol', DEFAULT_SYMBOL)
            ptsize = int(params.pop('ptsize', DEFAULT_PTSIZE))

            outline = params.pop('outline', DEFAULT_OUTLINE)
            fill = params.pop('fill', DEFAULT_FILL)

            self.add_points(foreground, area_def, points_list, font_file, font_size,
                            symbol, ptsize, outline, fill, **params)

        # Grids overlay
        if 'grid' in overlays:
            if 'major_lonlat' in overlays['grid'] or 'minor_lonlat' in overlays['grid']:
                Dlonlat = overlays['grid'].get('major_lonlat', (10.0, 10.0))
                dlonlat = overlays['grid'].get('minor_lonlat', (2.0, 2.0))
            else:
                Dlonlat = (overlays['grid'].get('lon_major', 10.0), overlays['grid'].get('lat_major', 10.0))
                dlonlat = (overlays['grid'].get('lon_minor', 2.0), overlays['grid'].get('lat_minor', 2.0))
            outline = overlays['grid'].get('outline', 'white')
            write_text = overlays['grid'].get('write_text', True)
            if isinstance(write_text, str):
                write_text = write_text.lower() in ['true', 'yes', '1', 'on']
            font = overlays['grid'].get('font', None)
            font_size = int(overlays['grid'].get('font_size', 10))
            fill = overlays['grid'].get('fill', outline)
            fill_opacity = overlays['grid'].get('fill_opacity', 255)
            if isinstance(font, str):
                if is_agg:
                    from aggdraw import Font
                    font = Font(fill, font, opacity=fill_opacity, size=font_size)
                else:
                    from PIL.ImageFont import truetype
                    font = truetype(font, font_size)
            minor_outline = overlays['grid'].get('minor_outline', 'white')
            minor_is_tick = overlays['grid'].get('minor_is_tick', True)
            if isinstance(minor_is_tick, str):
                minor_is_tick = minor_is_tick.lower() in ['true', 'yes', '1']
            lon_placement = overlays['grid'].get('lon_placement', 'tb')
            lat_placement = overlays['grid'].get('lat_placement', 'lr')

            grid_kwargs = {}
            if is_agg:
                width = float(overlays['grid'].get('width', 1.0))
                minor_width = float(overlays['grid'].get('minor_width', 0.5))
                outline_opacity = overlays['grid'].get('outline_opacity', 255)
                minor_outline_opacity = overlays['grid'].get('minor_outline_opacity', 255)
                grid_kwargs['width'] = width
                grid_kwargs['minor_width'] = minor_width
                grid_kwargs['outline_opacity'] = outline_opacity
                grid_kwargs['minor_outline_opacity'] = minor_outline_opacity

            self.add_grid(foreground, area_def, Dlonlat, dlonlat,
                          font=font, write_text=write_text, fill=fill,
                          outline=outline, minor_outline=minor_outline,
                          minor_is_tick=minor_is_tick,
                          lon_placement=lon_placement,
                          lat_placement=lat_placement,
                          **grid_kwargs)

        if cache_file is not None:
            self._write_and_apply_new_cached_image(cache_file, foreground, background)
        return foreground

    @staticmethod
    def _apply_cached_image(cache_file, cache_epoch, background, regenerate=False):
        try:
            config_time = cache_epoch or 0
            cache_time = os.path.getmtime(cache_file)
            # Cache file will be used only if it's newer than config file
            if config_time is not None and config_time < cache_time and not regenerate:
                foreground = Image.open(cache_file)
                logger.info('Using image in cache %s', cache_file)
                if background is not None:
                    background.paste(foreground, mask=foreground.split()[-1])
                return foreground
            logger.info("Regenerating cache file.")
        except OSError:
            logger.info("No overlay image found, new overlay image will be saved in cache.")
        return None

    @staticmethod
    def _write_and_apply_new_cached_image(cache_file, foreground, background):
        try:
            foreground.save(cache_file)
        except IOError as e:
            logger.error("Can't save cache: %s", str(e))
        if background is not None:
            background.paste(foreground, mask=foreground.split()[-1])

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
            if overlay_name == 'cache':
                continue
            params_to_hash[overlay_name] = overlay_dict.copy()
        # font objects are not hashable
        for font_cat in ('cities', 'points', 'grid'):
            if font_cat in params_to_hash:
                params_to_hash[font_cat].pop('font', None)
        return params_to_hash

    def add_overlay_from_config(self, config_file, area_def, background=None):
        """Create and return a transparent image adding all the overlays contained in a configuration file.

        :Parameters:
            config_file : str
                Configuration file name
            area_def : object
                Area Definition of the creating image

        """
        overlays = self._config_to_dict(config_file)
        return self.add_overlay_from_dict(overlays, area_def,
                                          os.path.getmtime(config_file), background)

    def add_cities(self, image, area_def, citylist, font_file, font_size,
                   ptsize, outline, box_outline, box_opacity, db_root_path=None):
        """Add cities (point and name) to a PIL image object."""
        if db_root_path is None:
            db_root_path = self.db_root_path
        if db_root_path is None:
            raise ValueError("'db_root_path' must be specified to use this method")

        draw = self._get_canvas(image)

        # read shape file with points
        # Sc-Kh shapefilename = os.path.join(self.db_root_path,
        # "cities_15000_alternativ.shp")
        shapefilename = os.path.join(
            db_root_path, os.path.join("CITIES",
                                       "cities_15000_alternativ.shp"))
        try:
            s = shapefile.Reader(shapefilename)
            shapes = s.shapes()
        except AttributeError:
            raise ValueError('Could not find shapefile %s'
                             % shapefilename)

        font = self._get_font(outline, font_file, font_size)

        # Iterate through shapes
        for i, shape in enumerate(shapes):
            # Select cities with name
            record = s.record(i)
            if record[3] in citylist:

                city_name = record[3]

                # use only parts of _get_pixel_index
                # Get shape data as array and reproject
                shape_data = np.array(shape.points)
                lons = shape_data[:, 0][0]
                lats = shape_data[:, 1][0]

                try:
                    (x, y) = area_def.get_xy_from_lonlat(lons, lats)
                except ValueError as exc:
                    logger.debug("Point not added (%s)", str(exc))
                else:
                    # add_dot
                    if ptsize is not None:
                        dot_box = [x - ptsize, y - ptsize,
                                   x + ptsize, y + ptsize]
                        self._draw_ellipse(
                            draw, dot_box, fill=outline, outline=outline)
                        text_position = [x + 9, y - 5]  # FIXME
                    else:
                        text_position = [x, y]

                # add text_box
                    self._draw_text_box(draw, text_position, city_name, font, outline,
                                        box_outline, box_opacity)
                    logger.info("%s added", str(city_name))

        self._finalize(draw)

    def add_points(self, image, area_def, points_list, font_file, font_size=12,
                   symbol='circle', ptsize=6, outline='black', fill='white', **kwargs):
        """Add a symbol and/or text at the point(s) of interest to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            area_def : object
                Area Definition of the provided image
            points_list : list [((lon, lat), desc)]
              | a list of points defined with (lon, lat) in float and a desc in string
              | [((lon1, lat1), desc1), ((lon2, lat2), desc2)]
              | lon : float
              |    longitude of a point
              | lat : float
              |    latitude of a point
              | desc : str
              |    description of a point
            font_file : str
                Path to font file
            font_size : int
                Size of font
            symbol : string
                type of symbol, one of the elelment from the list
                ['circle', 'square', 'asterisk']
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
        try:
            from pyresample.geometry import AreaDefinition
        except ImportError:
            raise ImportError("Missing required 'pyresample' module, please install it.")

        if not isinstance(area_def, AreaDefinition):
            raise ValueError("Expected 'area_def' is an instance of AreaDefinition object")

        draw = self._get_canvas(image)

        # Iterate through points list
        for point in points_list:
            (lon, lat), desc = point
            try:
                x, y = area_def.get_xy_from_lonlat(lon, lat)
            except ValueError:
                logger.info("Point %s is out of the area, it will not be added to the image.",
                            str((lon, lat)))
            else:
                if ptsize != 0:
                    half_ptsize = int(round(ptsize / 2.))

                    dot_box = [x - half_ptsize, y - half_ptsize,
                               x + half_ptsize, y + half_ptsize]

                    width = kwargs.get('width', 1.)
                    outline_opacity = kwargs.get('outline_opacity', 255)
                    fill_opacity = kwargs.get('fill_opacity', 0)

                    # draw the symbol at the (x, y) position
                    if symbol == 'circle':  # a 'circle' or a 'dot' i.e circle with fill
                        self._draw_ellipse(draw, dot_box,
                                           outline=outline, width=width,
                                           outline_opacity=outline_opacity,
                                           fill=fill, fill_opacity=fill_opacity)
                    elif symbol == 'square':
                        self._draw_rectangle(draw, dot_box,
                                             outline=outline, width=width,
                                             outline_opacity=outline_opacity,
                                             fill=fill, fill_opacity=fill_opacity)
                    elif symbol == 'asterisk':  # an '*' sign
                        self._draw_asterisk(draw, ptsize, (x, y),
                                            outline=outline, width=width,
                                            outline_opacity=outline_opacity)
                    else:
                        raise ValueError("Unsupported symbol type: " + str(symbol))

                elif desc is None:
                    logger.error("'ptsize' is 0 and 'desc' is None, nothing will be added to the image.")

                if desc is not None:
                    text_position = [x + ptsize, y]  # draw the text box next to the point
                    font = self._get_font(outline, font_file, font_size)

                    new_kwargs = kwargs.copy()

                    box_outline = new_kwargs.pop('box_outline', 'white')
                    box_opacity = new_kwargs.pop('box_opacity', 0)

                    # add text_box
                    self._draw_text_box(draw, text_position, desc, font, outline,
                                        box_outline, box_opacity, **new_kwargs)

            logger.debug("Point %s has been added to the image", str((lon, lat)))

        self._finalize(draw)


def _get_lon_lat_bounding_box(area_extent, x_size, y_size, prj):
    """Get extreme lon and lat values."""
    x_ll, y_ll, x_ur, y_ur = area_extent
    x_range = np.linspace(x_ll, x_ur, num=x_size)
    y_range = np.linspace(y_ll, y_ur, num=y_size)

    if prj.is_latlong():
        lons_s1, lats_s1 = x_ll * np.ones(y_range.size), y_range
        lons_s2, lats_s2 = x_range, y_ur * np.ones(x_range.size)
        lons_s3, lats_s3 = x_ur * np.ones(y_range.size), y_range
        lons_s4, lats_s4 = x_range, y_ll * np.ones(x_range.size)
    else:
        lons_s1, lats_s1 = prj(np.ones(y_range.size) * x_ll, y_range,
                               inverse=True)
        lons_s2, lats_s2 = prj(x_range, np.ones(x_range.size) * y_ur,
                               inverse=True)
        lons_s3, lats_s3 = prj(np.ones(y_range.size) * x_ur, y_range,
                               inverse=True)
        lons_s4, lats_s4 = prj(x_range, np.ones(x_range.size) * y_ll,
                               inverse=True)

    angle_sum = 0
    prev = None
    for lon in np.concatenate((lons_s1, lons_s2,
                               lons_s3[::-1], lons_s4[::-1])):
        if not np.isfinite(lon):
            continue
        if prev is not None:
            delta = lon - prev
            if abs(delta) > 180:
                delta = (abs(delta) - 360) * np.sign(delta)
            angle_sum += delta
        prev = lon

    if round(angle_sum) == -360:
        # Covers NP
        lat_min = min(lats_s1.min(), lats_s2.min(),
                      lats_s3.min(), lats_s4.min())
        lat_max = 90
        lon_min = -180
        lon_max = 180
    elif round(angle_sum) == 360:
        # Covers SP
        lat_min = -90
        lat_max = max(lats_s1.max(), lats_s2.max(),
                      lats_s3.max(), lats_s4.max())
        lon_min = -180
        lon_max = 180
    elif round(angle_sum) == 0:
        # Covers no poles
        if np.sign(lons_s1[0]) * np.sign(lons_s1[-1]) == -1 and \
                lons_s1.min() * lons_s1.max() < -25000:
            # End points of left side on different side of dateline
            lon_min = lons_s1[lons_s1 > 0].min()
        else:
            lon_min = lons_s1.min()

        if np.sign(lons_s3[0]) * np.sign(lons_s3[-1]) == -1 and \
                lons_s3.min() * lons_s3.max() < -25000:
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


def _get_pixel_index(shape, area_extent, x_size, y_size, prj,
                     x_offset=0, y_offset=0):
    """Map coordinates of shape to image coordinates."""
    # Get shape data as array and reproject
    shape_data = np.array(shape.points if hasattr(shape, 'points') else shape)
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
