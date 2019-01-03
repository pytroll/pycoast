#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pycoast, Writing of coastlines, borders, and rivers to images in Python
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

from PIL import Image, ImageFont
from PIL import ImageDraw
import logging

from pycoast.cw_base import ContourWriterBase

logger = logging.getLogger(__name__)


class ContourWriterPIL(ContourWriterBase):
    """Adds countours from GSHHS and WDBII to images.

    :Parameters:
        db_root_path : str
            Path to root dir of GSHHS and WDBII shapefiles

    """

    _draw_module = "PIL"
    # This  is a  flag  to make  _add_grid aware  of  which text  draw
    # routine from PIL  or from aggdraw should  be used (unfortunately
    # they are not fully compatible)

    def _get_canvas(self, image):
        """Returns PIL image object."""

        return ImageDraw.Draw(image)

    def _engine_text_draw(self, draw, x_pos, y_pos, txt, font, **kwargs):
        draw.text((x_pos, y_pos), txt, font=font, fill=kwargs['fill'])

    def _draw_polygon(self, draw, coordinates, **kwargs):
        """Draw polygon."""
        draw.polygon(coordinates, fill=kwargs['fill'], outline=kwargs['outline'])

    def _draw_ellipse(self, draw, coordinates, **kwargs):
        """Draw ellipse."""
        draw.ellipse(coordinates, fill=kwargs['fill'], outline=kwargs['outline'])

    def _draw_rectangle(self, draw, coordinates, **kwargs):
        """Draw rectangle."""
        draw.rectangle(coordinates, fill=kwargs['fill'], outline=kwargs['outline'])

    def _draw_text_box(self, draw, text_position, text, font, outline,
                       box_outline, box_opacity):
        """Add text box in xy."""
        if box_outline is not None:
            logger.warning(
                "Box background will not added; please install aggdraw lib")

        self._draw_text(
            draw, text_position, text, font, align="no", fill=outline)

    def _draw_line(self, draw, coordinates, **kwargs):
        """Draw line."""
        draw.line(coordinates, fill=kwargs['outline'])

    def add_shapefile_shapes(self, image, area_def, filename, feature_type=None,
                             fill=None, outline='white',
                             x_offset=0, y_offset=0):
        """Add shape file shapes from an ESRI shapefile.

        Note: Currently only supports lon-lat formatted coordinates.

        :Parameters:
            image : object
                PIL image object
            area_def : list [proj4_string, area_extent]
              | proj4_string : str
              |     Projection of area as Proj.4 string
              | area_extent : list
              |     Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            filename : str
                Path to ESRI shape file
            feature_type : 'polygon' or 'line',
                only to override the shape type defined in shapefile, optional
            fill : str or (R, G, B), optional
                Polygon fill color
            fill_opacity : int, optional {0; 255}
                Opacity of polygon fill
            outline : str or (R, G, B), optional
                line color
            outline_opacity : int, optional {0; 255}
                Opacity of lines
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """
        self._add_shapefile_shapes(image=image, area_def=area_def,
                                   filename=filename, feature_type=feature_type,
                                   x_offset=x_offset, y_offset=y_offset,
                                   fill=fill,
                                   outline=outline)

    def add_shapefile_shape(self, image, area_def, filename, shape_id,
                            feature_type=None,
                            fill=None, outline='white',
                            x_offset=0, y_offset=0):
        """Add a single shape file shape from an ESRI shapefile.

        Note: To add all shapes in file use the 'add_shape_file_shapes' routine.
        Note: Currently only supports lon-lat formatted coordinates.

        :Parameters:
            image : object
                PIL image object
            area_def : list [proj4_string, area_extent]
              | proj4_string : str
              |     Projection of area as Proj.4 string
              | area_extent : list
              |     Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            filename : str
                Path to ESRI shape file
            shape_id : int
                integer id of shape in shape file {0; ... }
            feature_type : 'polygon' or 'line',
                only to override the shape type defined in shapefile, optional
            fill : str or (R, G, B), optional
                Polygon fill color
            outline : str or (R, G, B), optional
                line color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """
        self._add_shapefile_shape(image=image,
                                  area_def=area_def, filename=filename,
                                  shape_id=shape_id,
                                  feature_type=feature_type,
                                  x_offset=x_offset, y_offset=y_offset,
                                  fill=fill,
                                  outline=outline)

    def add_line(self, image, area_def, lonlats,
                 fill=None, outline='white', x_offset=0, y_offset=0):
        """Add a user defined poly-line from a list of (lon,lat) coordinates.

        :Parameters:
            image : object
                PIL image object
            area_def : list [proj4_string, area_extent]
              | proj4_string : str
              |     Projection of area as Proj.4 string
              | area_extent : list
              |     Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            lonlats : list of lon lat pairs
                e.g. [(10,20),(20,30),...,(20,20)]
            outline : str or (R, G, B), optional
                line color
            width : float, optional
                line width
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """
        self._add_line(image, area_def, lonlats, x_offset=x_offset, y_offset=y_offset,
                       fill=fill, outline=outline)

    def add_polygon(self, image, area_def, lonlats,
                    fill=None, outline='white', x_offset=0, y_offset=0):
        """Add a user defined polygon from a list of (lon,lat) coordinates.

        :Parameters:
            image : object
                PIL image object
            area_def : list [proj4_string, area_extent]
              | proj4_string : str
              |     Projection of area as Proj.4 string
              | area_extent : list
              |     Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            lonlats : list of lon lat pairs
                e.g. [(10,20),(20,30),...,(20,20)]
            fill : str or (R, G, B), optional
                Polygon fill color
            outline : str or (R, G, B), optional
                line color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """
        self._add_polygon(image, area_def, lonlats, x_offset=x_offset,
                          y_offset=y_offset, fill=fill, outline=outline)

    def add_grid(self, image, area_def, Dlonlat, dlonlat,
                 font=None, write_text=True, fill=None, outline='white',
                 minor_outline='white', minor_is_tick=True,
                 lon_placement='tb', lat_placement='lr'):
        """Add a lon-lat grid to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            proj4_string : str
                Projection of area as Proj.4 string
            Dlonlat: (float, float)
                Major grid line separation
            dlonlat: (float, float)
                Minor grid line separation
            font: PIL ImageFont object, optional
                Font for major line markings
            write_text : boolean, optional
                Deterine if line markings are enabled
            fill : str or (R, G, B), optional
                Text color
            outline : str or (R, G, B), optional
                Major line color
            minor_outline : str or (R, G, B), optional
                Minor line/tick color
            minor_is_tick : boolean, optional
                Use tick minor line style (True) or full minor line style (False)

        """
        Dlon, Dlat = Dlonlat
        dlon, dlat = dlonlat
        self._add_grid(image, area_def, Dlon, Dlat, dlon, dlat, font=font,
                       write_text=write_text, fill=fill, outline=outline,
                       minor_outline=minor_outline, minor_is_tick=minor_is_tick,
                       lon_placement=lon_placement, lat_placement=lat_placement)

    def add_grid_to_file(self, filename, area_def, Dlonlat, dlonlat,
                         font=None, write_text=True, fill=None, outline='white',
                         minor_outline='white', minor_is_tick=True,
                         lon_placement='tb', lat_placement='lr'):
        """Add a lon-lat grid to an image file.

        :Parameters:
            image : object
                PIL image object
            proj4_string : str
                Projection of area as Proj.4 string
            Dlonlat: (float, float)
                Major grid line separation
            dlonlat: (float, float)
                Minor grid line separation
            font: PIL ImageFont object, optional
                Font for major line markings
            write_text : boolean, optional
                Deterine if line markings are enabled
            fill : str or (R, G, B), optional
                Text color
            outline : str or (R, G, B), optional
                Major line color
            minor_outline : str or (R, G, B), optional
                Minor line/tick color
            minor_is_tick : boolean, optional
                Use tick minor line style (True) or full minor line style (False)

        """
        image = Image.open(filename)
        self.add_grid(image, area_def, Dlonlat, dlonlat, font=font,
                      write_text=write_text, fill=fill, outline=outline,
                      minor_outline=minor_outline,
                      minor_is_tick=minor_is_tick,
                      lon_placement=lon_placement, lat_placement=lat_placement)
        image.save(filename)

    def add_coastlines(self, image, area_def, resolution='c', level=1,
                       fill=None, outline='white', x_offset=0, y_offset=0):
        """Add coastlines to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3, 4}
                Detail level of dataset
            fill : str or (R, G, B), optional
                Land color
            outline : str or (R, G, B), optional
                Coastline color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """

        self._add_feature(image, area_def, 'polygon', 'GSHHS',
                          resolution=resolution, level=level,
                          fill=fill, outline=outline, x_offset=x_offset,
                          y_offset=y_offset)

    def add_coastlines_to_file(self, filename, area_def, resolution='c',
                               level=1, fill=None, outline='white',
                               x_offset=0, y_offset=0):
        """Add coastlines to an image file.

        :Parameters:
            filename : str
                Image file
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3, 4}
                Detail level of dataset
            fill : str or (R, G, B)
                Land color
            outline : str or (R, G, B), optional
                Coastline color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """

        image = Image.open(filename)
        self.add_coastlines(image, area_def,
                            resolution=resolution, level=level,
                            fill=fill, outline=outline, x_offset=x_offset,
                            y_offset=y_offset)
        image.save(filename)

    def add_borders(self, image, area_def, resolution='c', level=1,
                    outline='white', x_offset=0, y_offset=0):
        """Add borders to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3}
                Detail level of dataset
            outline : str or (R, G, B), optional
                Border color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """

        self._add_feature(image, area_def, 'line', 'WDBII',
                          tag='border', resolution=resolution, level=level,
                          outline=outline, x_offset=x_offset,
                          y_offset=y_offset)

    def add_borders_to_file(self, filename, area_def, resolution='c', level=1,
                            outline='white', x_offset=0, y_offset=0):
        """Add borders to an image file.

        :Parameters:
            image : object
                Image file
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3}
                Detail level of dataset
            outline : str or (R, G, B), optional
                Border color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """
        image = Image.open(filename)
        self.add_borders(image, area_def, resolution=resolution,
                         level=level, outline=outline, x_offset=x_offset,
                         y_offset=y_offset)
        image.save(filename)

    def add_rivers(self, image, area_def, resolution='c', level=1,
                   outline='white', x_offset=0, y_offset=0):
        """Add rivers to a PIL image object.

        :Parameters:
            image : object
                PIL image object
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
                Detail level of dataset
            outline : str or (R, G, B), optional
                River color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """

        self._add_feature(image, area_def, 'line', 'WDBII',
                          tag='river', zero_pad=True, resolution=resolution,
                          level=level, outline=outline, x_offset=x_offset,
                          y_offset=y_offset)

    def add_rivers_to_file(self, filename, area_def, resolution='c', level=1,
                           outline='white', x_offset=0, y_offset=0):
        """Add rivers to an image file.

        :Parameters:
            image : object
                Image file
            proj4_string : str
                Projection of area as Proj.4 string
            area_extent : list
                Area extent as a list (LL_x, LL_y, UR_x, UR_y)
            resolution : str, optional {'c', 'l', 'i', 'h', 'f'}
                Dataset resolution to use
            level : int, optional {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}
                Detail level of dataset
            outline : str or (R, G, B), optional
                River color
            x_offset : float, optional
                Pixel offset in x direction
            y_offset : float, optional
                Pixel offset in y direction

        """

        image = Image.open(filename)
        self.add_rivers(image, area_def, resolution=resolution, level=level,
                        outline=outline, x_offset=x_offset, y_offset=y_offset)
        image.save(filename)

    def _get_font(self, outline, font_file, font_size):
        """Return a font."""
        if font_file is None:
            return ImageFont.load_default()
        return ImageFont.truetype(font_file, font_size)
