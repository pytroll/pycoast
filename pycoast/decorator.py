#
# Copyright (C) 2011  Hrobjartur Thorsteinsson
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

try:
    from PIL import ImageDraw
except ImportError:
    print("ImportError: Missing module: ImageDraw")


class Decorator(DecoratorBase):

    def add_scale(self, color_def, font=None, size=None, fill='black',
                  outline=None, outline_width=1, bg='white', extend=False, unit='', margins=None, minortick=0.0,
                  nan_color=(0, 0, 0), nan_check_color=(1, 1, 1), nan_check_size=0):
        """ Todo
        """
        self._add_scale(color_def, font=font, size=size, fill=fill,
                        outline=outline, outline_widht=outline_width, bg=bg, extend=extend, unit=unit, margins=margins, minortick=minortick,
                        nan_color=nan_color, nan_check_color=nan_check_color, nan_check_size=nan_check_size)

    def _load_default_font(self):
        return ImageFont.load_default()

    def add_text(self, txt, **kwargs):
        self._add_text(txt, **kwargs)

    def add_logo(self, logo_path, **kwargs):
        self._add_logo(logo_path, **kwargs)

    def _get_canvas(self, image):
        """Returns PIL image object
        """
        return ImageDraw.Draw(image)
