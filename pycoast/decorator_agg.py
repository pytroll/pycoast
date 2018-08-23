#
# Copyright (C) 2011, 2016  Hrobjartur Thorsteinsson
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
"""Tools for labeling and adding colour scales to images.
"""
from pycoast.decorator_base import DecoratorBase

try:
    from PIL import ImageDraw
except ImportError:
    print("ImportError: Missing module: ImageDraw")


class DecoratorAGG(DecoratorBase):

    def add_scale(self, colormap, **kwargs):
        self._add_scale(colormap, **kwargs)

    def _load_default_font(self):
        import aggdraw
        from pkg_resources import resource_filename
        font_path = resource_filename('pycoast.fonts', 'DejaVuSerif.ttf')
        return aggdraw.Font('black', font_path, size=16)

    def add_text(self, txt, **kwargs):
        self._add_text(txt, **kwargs)

    def add_logo(self, logo_path, **kwargs):
        self._add_logo(logo_path, **kwargs)

    def _get_canvas(self, image):
        """Returns AGG image object
        """
        import aggdraw
        return aggdraw.Draw(image)

    def _finalize(self, draw):
        """Flush the AGG image object
        """
        draw.flush()

    def _draw_text_line(self, draw, xy, text, font, fill='black'):
        draw.text(xy, text, font)

    def _draw_rectangle(self, draw, xys, outline='white', bg='white', bg_opacity=255, outline_width=1, outline_opacity=255, **kwargs):
        import aggdraw
        pen = aggdraw.Pen(
            outline, width=outline_width, opacity=outline_opacity)
        brush = aggdraw.Brush(bg, opacity=bg_opacity)
        # draw bg and outline
        # bg unaliased (otherwise gaps between successive bgs)

        if bg is not None:
            draw.setantialias(False)
            draw.rectangle(xys, None, brush)
            draw.setantialias(True)
        # adjust to correct for outline exceeding requested area,
        # due to outline width expanding outwards.
        xys[0] += outline_width / 2.0
        xys[1] += outline_width / 2.0
        xys[2] -= outline_width / 2.0
        xys[3] -= outline_width / 2.0
        if outline is not None:
            draw.rectangle(xys, pen, None)

    def _draw_line(self, draw, xys, **kwargs):
        import aggdraw
        if kwargs['line'] is None:
            pen = None
        else:
            pen = aggdraw.Pen(
                kwargs['line'], width=kwargs['line_width'], opacity=kwargs['line_opacity'])
        xys_straight = [item for t in xys for item in t]
        draw.line(xys_straight, pen)

    def _draw_polygon(self, draw, xys, outline=None, fill='white', fill_opacity=255, outline_width=1, outline_opacity=255):
        import aggdraw
        if outline is None:
            pen = None
        else:
            pen = aggdraw.Pen(
                outline, width=outline_width, opacity=outline_opacity)
        if fill is None:
            brush = None
        else:
            brush = aggdraw.Brush(fill, opacity=fill_opacity)
        xys_straight = [item for t in xys for item in t]
        draw.polygon(xys_straight, pen, brush)


#########################################
# float list generator
def _frange(x, y, jump):
    while x < y:
        yield x
        x += jump


def frange(x, y, jump):
    return [p for p in _frange(x, y, jump)]


class ShapeFileError(Exception):
    pass
