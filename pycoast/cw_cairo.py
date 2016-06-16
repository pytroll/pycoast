#!/usr/bin/env python
# -*- coding: utf-8 -*-
# pycoast, Writing of coastlines, borders and rivers to images in Python
#
# Copyright (C) 2011-2016
#    Esben S. Nielsen
#    Hróbjartur Þorsteinsson
#    Stefano Cerino
#    Katja Hungershofer
#    Panu Lahtinen
#    Sauli Joro
#    Adam Dybbroe
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

import pdb

import numpy as np
from PIL import Image
import logging

logger = logging.getLogger(__name__)

# Prefer cairocffi, but use API compatible pycairo if not available.
# Why do we prefer cairocffi? FIXME
# try:
#     import cairocffi as cairo
# except ImportError:
#     logger.error('cairocffi is needed')
#     raise


from cairo import FORMAT_ARGB32, Context, ImageSurface

from .cw_base import ContourWriterBase


class CairoDrawObject(object):

    def __init__(self, context, pilimage):
        self.context = context
        self.pilimage = pilimage


# Doesn't work with cairo, but works fine with cairocffi
# class CustomImageSurface(ImageSurface):

#     def __init__(self):
#         ImageSurface.__init__(self)
#         self.size = (self.get_width(), self.get_height())


class ContourWriterCairo(ContourWriterBase):

    """Adds countours from GSHHS and WDBII to images using the Cairo
       engine for high quality images.

    :Parameters:
    db_root_path : str
        Path to root dir of GSHHS and WDBII shapefiles
    """
    _draw_module = "Cairo"
    # This is a flag to make _add_grid aware of which text draw
    # routine from PIL, aggdraw or Cairo should be used (unfortunately
    # they are not fully compatible)

    # Map the color from string to R, G, B
    colormap = {'white': (255, 255, 255),
                'blue': (0, 0, 255),
                'black': (0, 0, 0),
                'grey': (179, 179, 179),
                'red': (255, 0, 0),
                'green': (0, 255, 0)}

    def _get_canvas(self, image):
        """Returns PIL draw canvass, convert from Cairo image context if
        necessary.
        """

        if isinstance(image, Context):
            # Create PIL image from cairo context:
            raise NotImplementedError(
                'Cairo context as input not supported!')
        elif isinstance(image, ImageSurface):
            # Create PIL image from cairo image surface:
            self._surface = image
            pil_img = self._cairo_to_pil()
            canvas = CairoDrawObject(Context(image), pil_img)
            return canvas

        elif isinstance(image, Image.Image):
            self._surface = _pil_to_cairo(image)
            # fmt = FORMAT_ARGB32
            # width, height = image.size
            # data = image.convert('RGBA')
            # data = np.asarray(data)
            # tmp = data.copy()
            # self._surface = ImageSurface.create_for_data(tmp, FORMAT_ARGB32,
            #                                              width, height)
            canvas = CairoDrawObject(Context(self._surface), image)
            return canvas
        else:
            raise ValueError("Unsupported image format.")

    def _engine_text_draw(self, img_ctx, (x_pos, y_pos), txt, font, **kwargs):

        #outline = kwargs.get('outline')
        # if type(outline) is str:
        #    r, g, b = self.colormap.get(outline)
        # else:
        #    r, g, b = outline

        # img_ctx.set_source_rgb(r, g, b)  # Solid color
        # img_ctx.set_line_width(kwargs['width'])
        # img_ctx.stroke()

        img_ctx.move_to(x_pos, y_pos)
        img_ctx.show_text(txt)

    def _draw_polygon(self, draw_obj, coordinates, **kwargs):
        """Draw polygon
        """

        alpha = kwargs.get('outline_opacity', 255)
        alpha /= 255.
        outline = kwargs.get('outline')
        r, g, b = get_rgb_values(outline, self.colormap)

        img_ctx = draw_obj.context

        x, y = coordinates[0]
        img_ctx.move_to(x, y)
        for x, y in coordinates[1:]:
            img_ctx.line_to(x, y)

        img_ctx.set_source_rgba(r, g, b, alpha)  # Solid color
        img_ctx.set_line_width(kwargs['width'])
        img_ctx.stroke()

    def _draw_rectangle(self, draw_obj, coordinates, **kwargs):
        """Draw rectangle
        """
        self._draw_polygon(draw_obj, coordinates, **kwargs)

    def _draw_ellipse(self, draw_obj, coordinates, **kwargs):
        """Draw ellipse
        """
        raise NotImplementedError('Ellipse drawing not implemented for Cairo '
                                  'backend')

    def _draw_text(self, draw_obj, position, txt, font, align='cc', **kwargs):
        """Draw text with Cairo module
        """
        img_ctx = draw_obj.context

        self._set_font(img_ctx, font)
        (xbearing, ybearing, txt_width, txt_height,
         x_advance, y_advance) = img_ctx.text_extents(txt)

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

        self._engine_text_draw(img_ctx, (x_pos, y_pos), txt, font, **kwargs)

    def _draw_text_box(self, draw_obj, text_position, text, font, outline,
                       box_outline, box_opacity):
        """Add text box in xy
        """
        img_ctx = draw_obj.context

        self._set_font(img_ctx, font)
        if box_outline is not None:
            (xbearing, ybearing, txt_width, txt_height,
             x_advance, y_advance) = img_ctx.text_extents(text)
            margin = 2
            xUL = text_position[0] - margin
            yUL = text_position[1]
            xLR = xUL + txt_width + (2 * margin)
            yLR = yUL + txt_height
            box = ((xUL, yUL), (xLR, yUL), (xLR, yLR), (xUL, yLR), (xUL, yUL))

            self._draw_rectangle(img_ctx, box, fill=box_outline,
                                 fill_opacity=box_opacity, outline=box_outline)

        self._draw_text(img_ctx, text_position, text, font, align="no")

    def _draw_line(self, draw_obj, coordinates, **kwargs):
        """Draw line
        """
        self._draw_polygon(draw_obj, coordinates, **kwargs)

    def _set_font(self, img_ctx, font):
        """Set font."""
        pass
        # raise NotImplementedError('Font things are not implemented.')

    def _finalize(self, draw_obj):
        """Finish the draw object (contains cairo image context object)
        """
        img_ctx = draw_obj.context
        img = self._cairo_to_pil()
        draw_obj.pilimage.paste(img)

        # draw_obj.pilimage.save('/tmp/kurt.png')
        # self._surface.write_to_png('/tmp/cairo_surface.png')

    def add_shapefile_shapes(self, img_ctx, area_def, filename,
                             feature_type=None, fill=None, fill_opacity=255,
                             outline='white', width=1, outline_opacity=255,
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
        width : float, optional
            line width
        outline_opacity : int, optional {0; 255}
            Opacity of lines
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """
        self._add_shapefile_shapes(img_ctx, area_def,
                                   filename, feature_type=feature_type,
                                   x_offset=x_offset, y_offset=y_offset,
                                   fill=fill, fill_opacity=fill_opacity,
                                   outline=outline, width=width,
                                   outline_opacity=outline_opacity)

    def add_shapefile_shape(self, img_ctx, area_def, filename, shape_id,
                            feature_type=None,
                            fill=None, fill_opacity=255, outline='white',
                            width=1, outline_opacity=255, x_offset=0,
                            y_offset=0):
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
        fill_opacity : int, optional {0; 255}
            Opacity of polygon fill
        outline : str or (R, G, B), optional
            line color
        width : float, optional
            line width
        outline_opacity : int, optional {0; 255}
            Opacity of lines
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """
        self._add_shapefile_shape(img_ctx, area_def, filename, shape_id,
                                  feature_type=feature_type,
                                  x_offset=x_offset, y_offset=y_offset,
                                  fill=fill, fill_opacity=fill_opacity,
                                  outline=outline, width=width,
                                  outline_opacity=outline_opacity)

    def add_line(self, img_ctx, area_def, lonlats,
                 fill=None, fill_opacity=255, outline='white', width=1,
                 outline_opacity=255, x_offset=0, y_offset=0):
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
        outline_opacity : int, optional {0; 255}
            Opacity of lines
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """
        self._add_line(img_ctx, area_def, lonlats, x_offset=x_offset,
                       y_offset=y_offset, fill=fill, fill_opacity=fill_opacity,
                       outline=outline, width=width,
                       outline_opacity=outline_opacity)

    def add_polygon(self, img_ctx, area_def, lonlats,
                    fill=None, fill_opacity=255, outline='white', width=1,
                    outline_opacity=255, x_offset=0, y_offset=0):
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
        fill_opacity : int, optional {0; 255}
            Opacity of polygon fill
        outline : str or (R, G, B), optional
            line color
        width : float, optional
            line width
        outline_opacity : int, optional {0; 255}
            Opacity of lines
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """
        self._add_polygon(img_ctx, area_def, lonlats, x_offset=x_offset,
                          y_offset=y_offset, fill=fill,
                          fill_opacity=fill_opacity, outline=outline,
                          width=width, outline_opacity=outline_opacity)

    def add_grid(self, img_surf, area_def, (Dlon, Dlat), (dlon, dlat),
                 font=None, write_text=True, fill=None, fill_opacity=255,
                 outline='white', width=1, outline_opacity=255,
                 minor_outline='white', minor_width=0.5,
                 minor_outline_opacity=255, minor_is_tick=True,
                 lon_placement='tb', lat_placement='lr'):
        """Add a lon-lat grid to a PIL image object

        :Parameters:
        image : object
            PIL image object
        proj4_string : str
            Projection of area as Proj.4 string
        (Dlon,Dlat): (float,float)
            Major grid line separation
        (dlon,dlat): (float,float)
            Minor grid line separation
        font: Aggdraw Font object, optional
            Font for major line markings
        write_text : boolean, optional
            Deterine if line markings are enabled
        fill_opacity : int, optional {0; 255}
            Opacity of text
        outline : str or (R, G, B), optional
            Major line color
        width : float, optional
            Major line width
        outline_opacity : int, optional {0; 255}
            Opacity of major lines
        minor_outline : str or (R, G, B), optional
            Minor line/tick color
        minor_width : float, optional
            Minor line width
        minor_outline_opacity : int, optional {0; 255}
            Opacity of minor lines/ticks
        minor_is_tick : boolean, optional
            Use tick minor line style (True) or full minor line style (False)
        """
        self._add_grid(img_surf, area_def, Dlon, Dlat, dlon, dlat,
                       font=font, write_text=write_text,
                       fill=fill, fill_opacity=fill_opacity, outline=outline,
                       width=width, outline_opacity=outline_opacity,
                       minor_outline=minor_outline, minor_width=minor_width,
                       minor_outline_opacity=minor_outline_opacity,
                       minor_is_tick=minor_is_tick,
                       lon_placement=lon_placement, lat_placement=lat_placement)

    def add_grid_to_file(self, filename, area_def, (Dlon, Dlat), (dlon, dlat),
                         font=None, write_text=True,
                         fill=None, fill_opacity=255,
                         outline='white', width=1, outline_opacity=255,
                         minor_outline='white', minor_width=0.5,
                         minor_outline_opacity=255, minor_is_tick=True,
                         lon_placement='tb', lat_placement='lr'):
        """Add a lon-lat grid to an image

        :Parameters:
        image : object
            PIL image object
        proj4_string : str
            Projection of area as Proj.4 string
        (Dlon,Dlat): (float,float)
            Major grid line separation
        (dlon,dlat): (float,float)
            Minor grid line separation
        font: Aggdraw Font object, optional
            Font for major line markings
        write_text : boolean, optional
            Deterine if line markings are enabled
        fill_opacity : int, optional {0; 255}
            Opacity of text
        outline : str or (R, G, B), optional
            Major line color
        width : float, optional
            Major line width
        outline_opacity : int, optional {0; 255}
            Opacity of major lines
        minor_outline : str or (R, G, B), optional
            Minor line/tick color
        minor_width : float, optional
            Minor line width
        minor_outline_opacity : int, optional {0; 255}
            Opacity of minor lines/ticks
        minor_is_tick : boolean, optional
            Use tick minor line style (True) or full minor line style (False)
        """

        img_surface = _read_image(filename)

        self.add_grid(img_surface, area_def, (Dlon, Dlat), (dlon, dlat),
                      font=font, write_text=write_text,
                      fill=fill, fill_opacity=fill_opacity,
                      outline=outline, width=width,
                      outline_opacity=outline_opacity,
                      minor_outline=minor_outline,
                      minor_width=minor_width,
                      minor_outline_opacity=minor_outline_opacity,
                      minor_is_tick=minor_is_tick,
                      lon_placement=lon_placement, lat_placement=lat_placement)
        _save_image(filename, img_surface)

    def add_coastlines(self, img_surface, area_def, resolution='c', level=1,
                       fill=None, fill_opacity=255, outline='white', width=1,
                       outline_opacity=255, x_offset=0, y_offset=0, **kwargs):
        """Add coastlines to a PIL image object

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
        fill_opacity : int, optional {0; 255}
            Opacity of land color
        outline : str or (R, G, B), optional
            Coastline color
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        self._add_feature(img_surface, area_def, 'polygon', 'GSHHS',
                          resolution=resolution, level=level,
                          fill=fill, fill_opacity=fill_opacity,
                          outline=outline, width=width,
                          outline_opacity=outline_opacity, x_offset=x_offset,
                          y_offset=y_offset, **kwargs)

    def add_coastlines_to_file(self, filename, area_def, resolution='c',
                               level=1, fill=None, fill_opacity=255,
                               outline='white', width=1, outline_opacity=255,
                               x_offset=0, y_offset=0, **kwargs):
        """Add coastlines to an image file

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
        fill : str or (R, G, B), optional
            Land color
        fill_opacity : int, optional {0; 255}
            Opacity of land color
        outline : str or (R, G, B), optional
            Coastline color
        width : float, o ptional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        img_surf = _read_image(filename)
        self.add_coastlines(img_surf, area_def, resolution=resolution,
                            level=level, fill=fill,
                            fill_opacity=fill_opacity, outline=outline,
                            width=width, outline_opacity=outline_opacity,
                            x_offset=x_offset, y_offset=y_offset, **kwargs)
        _save_image(filename, img_surf)

    def add_borders(self, img_surf, area_def, resolution='c', level=1,
                    outline='white', width=1, outline_opacity=255,
                    x_offset=0, y_offset=0, **kwargs):
        """Add borders to a PIL image object

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
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        self._add_feature(img_surf, area_def, 'line', 'WDBII', tag='border',
                          resolution=resolution, level=level, outline=outline,
                          width=width, outline_opacity=outline_opacity,
                          x_offset=x_offset, y_offset=y_offset, **kwargs)

    def add_borders_to_file(self, filename, area_def, resolution='c',
                            level=1, outline='white', width=1,
                            outline_opacity=255, x_offset=0, y_offset=0,
                            **kwargs):
        """Add borders to an image file

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
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        img_surf = _read_image(filename)
        self.add_borders(img_surf, area_def, resolution=resolution, level=level,
                         outline=outline, width=width,
                         outline_opacity=outline_opacity, x_offset=x_offset,
                         y_offset=y_offset, **kwargs)
        _save_image(filename, img_surf)

    def add_rivers(self, img_surf, area_def, resolution='c', level=1,
                   outline='white', width=1, outline_opacity=255,
                   x_offset=0, y_offset=0, **kwargs):
        """Add rivers to a PIL image object

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
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        self._add_feature(img_surf, area_def, 'line', 'WDBII', tag='river',
                          zero_pad=True, resolution=resolution, level=level,
                          outline=outline, width=width,
                          outline_opacity=outline_opacity, x_offset=x_offset,
                          y_offset=y_offset, **kwargs)

    def add_rivers_to_file(self, filename, area_def, resolution='c', level=1,
                           outline='white', width=1, outline_opacity=255,
                           x_offset=0, y_offset=0, **kwargs):
        """Add rivers to an image file

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
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """

        img_surf = _read_image(filename)
        self.add_rivers(img_surf, area_def, resolution=resolution, level=level,
                        outline=outline, width=width,
                        outline_opacity=outline_opacity, x_offset=x_offset,
                        y_offset=y_offset)
        _save_image(filename, img_surf)

    def _get_font(self, outline, font_file, font_size):
        """Return a font."""
        raise NotImplementedError("Get font not implemented")
        # return aggdraw.Font(outline, font_file, size=font_size)

    def _cairo_to_pil(self, out_mode='RGB'):
        '''Convert Cairo image context to PIL image.'''
        # Finalize the image to fix changes

        cairo_format = self._surface.get_format()
        if cairo_format == FORMAT_ARGB32:
            pil_mode = 'RGB'
            # Cairo buffer is supposed to be ARGB, but seems to be RGBA.
            # Convert this to RGB for PIL which supports
            # only RGB or RGBA. Thus, get rid of the last column (the alpha
            # layer)
            argb_array = np.fromstring(
                self._surface.get_data(), 'c').reshape(-1, 4)
            rgb_array = argb_array[:, 2::-1]
            #rgb_array = argb_array[:, :3]
            pil_data = rgb_array.reshape(-1).tostring()
        else:
            raise ValueError('Unsupported cairo format: %d' % cairo_format)
        pil_image = Image.frombuffer(pil_mode,
                                     (self._surface.get_width(),
                                      self._surface.get_height()),
                                     pil_data, "raw", pil_mode, 0, 1)
        return pil_image.convert(out_mode)


def _pil_to_cairo(pil_img):
    '''Convert PIL image to cairo surface'''

    img_rgba = pil_img.convert('RGBA')
    data = np.array(img_rgba.tobytes('raw', 'BGRA'))
    width, height = img_rgba.size
    surface = ImageSurface.create_for_data(data, FORMAT_ARGB32,
                                           width, height)
    return surface


def _cairo_to_pil(surface, out_mode='RGB'):
    '''Convert Cairo image context to PIL image.'''
    # Finalize the image to fix changes

    cairo_format = surface.get_format()
    if cairo_format == FORMAT_ARGB32:
        pil_mode = 'RGB'
        # Cairo has ARGB. Convert this to RGB for PIL which supports
        # only RGB or RGBA.
        argb_array = np.fromstring(surface.get_data(), 'c').reshape(-1, 4)
        rgb_array = argb_array[:, 2::-1]
        pil_data = rgb_array.reshape(-1).tostring()
    else:
        raise ValueError('Unsupported cairo format: %d' % cairo_format)
    pil_image = Image.frombuffer(pil_mode,
                                 (surface.get_width(), surface.get_height()),
                                 pil_data, "raw", pil_mode, 0, 1)
    return pil_image.convert(out_mode)


def get_rgb_values(param, colormap):
    # Convert color RGB values to be between 0 and 1:
    if type(param) is str:
        if param not in colormap:
            logger.warning("Color %s not supported", param)
        r, g, b = colormap.get(param, (255, 255, 255))
    else:
        r, g, b = param

    return r / 255., g / 255., b / 255.


def _read_image(fname):
    """Read the given file and create a cairo context."""
    img = Image.open(fname)
    return _pil_to_cairo(img)


def _save_image(fname, img_surf, out_mode='RGB'):
    """Save the given cairo context to a file."""
    img = _cairo_to_pil(img_surf, out_mode)
    img.save(fname)
