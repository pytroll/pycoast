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

import numpy as np
from PIL import Image
import ImageDraw
import shapefile
import pyproj


class ShapeFileError(Exception):
    pass


class ContourWriterBase(object):
    """Base class for contourwriters. Do not instantiate.
    
    :Parameters:
    db_root_path : str
        Path to root dir of GSHHS and WDBII shapefiles
    """

    def __init__(self, db_root_path):
        self.db_root_path = db_root_path

    def _add_feature(self, image, area_def, feature_type, 
                     db_name, tag=None, zero_pad=False, resolution='c', 
                     level=1, x_offset=0, y_offset=0, **kwargs):
        """Add a contour feature to image
        """

        try:
            proj4_string = area_def.proj4_string
            area_extent = area_def.area_extent
        except AttributeError:
            proj4_string = area_def[0]
            area_extent = area_def[1]

        draw = self._get_canvas(image)
        
        # Area and projection info
        x_size, y_size = image.size        
        prj = pyproj.Proj(proj4_string)
        
        
        # Calculate min and max lons and lats of interest        
        lon_min, lon_max, lat_min, lat_max = \
                _get_lon_lat_bounding_box(area_extent, x_size, y_size, prj) 
        
        # Iterate through detail levels        
        for shapes in self._iterate_db(db_name, tag, resolution, 
                                       level, zero_pad):

            # Iterate through shapes
            for i, shape in enumerate(shapes):
                # Check if polygon is possibly relevant
                s_lon_ll, s_lat_ll, s_lon_ur, s_lat_ur = shape.bbox
                if (lon_max < s_lon_ll or lon_min > s_lon_ur or 
                    lat_max < s_lat_ll or lat_min > s_lat_ur):
                    # Polygon is irrelevant
                    continue          
                
                # Get pixel index coordinates of shape
                
                index_arrays, is_reduced = _get_pixel_index(shape, area_extent, 
                                                            x_size, y_size, 
                                                            prj, 
                                                            x_offset=x_offset,
                                                            y_offset=y_offset)       
                
                # Skip empty datasets               
                if len(index_arrays) == 0:
                    continue

                # Make PIL draw the polygon or line
                for index_array in index_arrays:
                    if feature_type.lower() == 'polygon' and not is_reduced:
                        # Draw polygon if dataset has not been reduced
                        #draw.polygon(index_array.flatten().tolist(), fill=fill, 
                        #             outline=outline)
                        self._draw_polygon(draw, index_array.flatten().tolist(), **kwargs)
                    elif feature_type.lower() == 'line' or is_reduced:
                        # Draw line
                        self._draw_line(draw, index_array.flatten().tolist(), **kwargs)
                        #draw.line(index_array.flatten().tolist(), fill=outline)
                    else:
                        raise ValueError('Unknown contour type: %s' % feature_type)
                        
        self._finalize(draw)

    def _iterate_db(self, db_name, tag, resolution, level, zero_pad):
        """Iterate trough datasets
        """
        
        format_string = '%s_%s_'
        if tag is not None:
            format_string += '%s_'
            
        if zero_pad:
            format_string += 'L%02i.shp' 
        else:
            format_string += 'L%s.shp'
            
        for i in range(level):
            
            # One shapefile per level
            if tag is None:
                shapefilename = \
                        os.path.join(self.db_root_path, '%s_shp' % db_name, 
                                     resolution, format_string % 
                                     (db_name, resolution, (i + 1)))
            else:
                shapefilename = \
                        os.path.join(self.db_root_path, '%s_shp' % db_name, 
                                     resolution, format_string % 
                                     (db_name, tag, resolution, (i + 1)))
            try:
                s = shapefile.Reader(shapefilename)
                shapes = s.shapes()
            except AttributeError:
                raise ShapeFileError('Could not find shapefile %s' % shapefilename)
            yield shapes
            
    def _finalize(self, draw):
        """Do any need finalization of the drawing
        """
        
        pass    


class ContourWriter(ContourWriterBase):
    """Adds countours from GSHHS and WDBII to images
    
    :Parameters:
    db_root_path : str
        Path to root dir of GSHHS and WDBII shapefiles
    """
    
    def _get_canvas(self, image):
        """Returns PIL image object
        """
        
        return ImageDraw.Draw(image)
        
    def _draw_polygon(self, draw, coordinates, **kwargs):
        """Draw polygon
        """
        
        draw.polygon(coordinates, fill=kwargs['fill'], outline=kwargs['outline'])
        
    def _draw_line(self, draw, coordinates, **kwargs):
        """Draw line
        """
        
        draw.line(coordinates, fill=kwargs['outline'])
           
    def add_coastlines(self, image, area_def, resolution='c', level=1, 
                       fill=None, outline='white', x_offset=0, y_offset=0):
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
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction
        """
        
        image = Image.open(filename)
        self.add_rivers(image, area_def, resolution=resolution, level=level, 
                        outline=outline, x_offset=x_offset, y_offset=y_offset)
        image.save(filename)



class ContourWriterAGG(ContourWriterBase):
    """Adds countours from GSHHS and WDBII to images 
       using the AGG engine for high quality images.
    
    :Parameters:
    db_root_path : str
        Path to root dir of GSHHS and WDBII shapefiles
    """
    
    def _get_canvas(self, image):
        """Returns AGG image object
        """
        
        import aggdraw
        return aggdraw.Draw(image)
        
    def _draw_polygon(self, draw, coordinates, **kwargs):
        """Draw polygon
        """
        
        import aggdraw
        pen = aggdraw.Pen(kwargs['outline'], kwargs['width'], kwargs['outline_opacity'])
        if kwargs['fill'] is None:
            fill_opacity = 0
        else:
            fill_opacity = kwargs['fill_opacity']
        brush = aggdraw.Brush(kwargs['fill'], fill_opacity)
        draw.polygon(coordinates, pen, brush)
        
    def _draw_line(self, draw, coordinates, **kwargs):
        """Draw line
        """
        
        import aggdraw
        pen = aggdraw.Pen(kwargs['outline'], kwargs['width'], kwargs['outline_opacity'])
        draw.line(coordinates, pen)
        
    def _finalize(self, draw):
        """Flush the AGG image object
        """
        
        draw.flush()
           
    def add_coastlines(self, image, area_def, resolution='c', level=1, 
                       fill=None, fill_opacity=255, outline='white', width=1, 
                       outline_opacity=255, x_offset=0, y_offset=0):
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
        
        self._add_feature(image, area_def, 'polygon', 'GSHHS', 
                          resolution=resolution, level=level, 
                          fill=fill, fill_opacity=fill_opacity, 
                          outline=outline, width=width,
                          outline_opacity=outline_opacity, x_offset=x_offset,
                          y_offset=y_offset)
                              
    def add_coastlines_to_file(self, filename, area_def, resolution='c', 
                               level=1, fill=None, fill_opacity=255, 
                               outline='white', width=1, outline_opacity=255, 
                               x_offset=0, y_offset=0):
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
        width : float, optional
            Width of coastline
        outline_opacity : int, optional {0; 255}
            Opacity of coastline color
        x_offset : float, optional
            Pixel offset in x direction
        y_offset : float, optional
            Pixel offset in y direction      
        """
        
        image = Image.open(filename)
        self.add_coastlines(image, area_def, resolution=resolution, 
                            level=level, fill=fill, 
                            fill_opacity=fill_opacity, outline=outline, 
                            width=width, outline_opacity=outline_opacity, 
                            x_offset=x_offset, y_offset=y_offset)
        image.save(filename)

    def add_borders(self, image, area_def, resolution='c', level=1, 
                    outline='white', width=1, outline_opacity=255, 
                    x_offset=0, y_offset=0):
                            
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
        
        self._add_feature(image, area_def, 'line', 'WDBII', tag='border', 
                          resolution=resolution, level=level, outline=outline, 
                          width=width, outline_opacity=outline_opacity, 
                          x_offset=x_offset, y_offset=y_offset)

    def add_borders_to_file(self, filename, area_def, resolution='c', 
                            level=1, outline='white', width=1, 
                            outline_opacity=255, x_offset=0, y_offset=0):
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
        image = Image.open(filename)
        self.add_borders(image, area_def, resolution=resolution, level=level, 
                         outline=outline, width=width, 
                         outline_opacity=outline_opacity, x_offset=x_offset,
                         y_offset=y_offset)
        image.save(filename)
        
    def add_rivers(self, image, area_def, resolution='c', level=1, 
                   outline='white', width=1, outline_opacity=255, 
                   x_offset=0, y_offset=0):
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
        
        self._add_feature(image, area_def, 'line', 'WDBII', tag='river', 
                          zero_pad=True, resolution=resolution, level=level, 
                          outline=outline, width=width, 
                          outline_opacity=outline_opacity, x_offset=x_offset,
                          y_offset=y_offset)
                          
    def add_rivers_to_file(self, filename, area_def, resolution='c', level=1, 
                           outline='white', width=1, outline_opacity=255, 
                           x_offset=0, y_offset=0):
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
        
        image = Image.open(filename)
        self.add_rivers(image, area_def, resolution=resolution, level=level, 
                        outline=outline, width=width, 
                        outline_opacity=outline_opacity, x_offset=x_offset,
                        y_offset=y_offset)
        image.save(filename)


def _get_lon_lat_bounding_box(area_extent, x_size, y_size, prj):
    """Get extreme lon and lat values
    """
            
    x_ll, y_ll, x_ur, y_ur = area_extent
    x_range = np.linspace(x_ll, x_ur, num=x_size)
    y_range = np.linspace(y_ll, y_ur, num=y_size)
    
    lons_s1, lats_s1 = prj(np.ones(y_range.size) * x_ll, y_range, inverse=True)
    lons_s2, lats_s2 = prj(x_range, np.ones(x_range.size) * y_ur, inverse=True)
    lons_s3, lats_s3 = prj(np.ones(y_range.size) * x_ur, y_range, inverse=True)
    lons_s4, lats_s4 = prj(x_range, np.ones(x_range.size) * y_ll, inverse=True)
    
    angle_sum = 0
    prev = None
    for lon in np.concatenate((lons_s1, lons_s2, lons_s3[::-1], lons_s4[::-1])):
        if prev is not None:
            delta = lon - prev
            if abs(delta) > 180:
                delta = (abs(delta) - 360) * np.sign(delta)
            angle_sum += delta
        prev = lon
 
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
        lon_min = lons_s1.min()
        lon_max = lons_s3.max()
        lat_min = lats_s4.min()
        lat_max = lats_s2.max()
    else:
        # Pretty strange
        lat_min = -90
        lat_max = 90
        lon_min = -180
        lon_max = 180

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
    """Map coordinates of shape to image coordinates
    """
    
    # Get shape data as array and reproject    
    shape_data = np.array(shape.points)
    lons = shape_data[:, 0]
    lats = shape_data[:, 1]

    x_ll, y_ll, x_ur, y_ur = area_extent

    x, y = prj(lons, lats)
 
    #Handle out of bounds
    i = 0
    segments = []
    if 1e30 in x or 1e30 in y:
        # Split polygon in line segments within projection
        is_reduced = True
        if x[0] == 1e30 or y[0] == 1e30:
            in_segment = False
        else:
            in_segment = True
            
        for j in range(x.size):
            if (x[j] == 1e30 or y[j] == 1e30):
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
                        
