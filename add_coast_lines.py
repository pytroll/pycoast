import os

import numpy as np
from PIL import Image
import ImageDraw
import shapefile
import pyproj


class ContourWriter(object):

    def __init__(self, db_root_path):
        self.db_root_path = db_root_path

    def add_coastlines(self, image, proj4_string, area_extent, resolution='c', 
                       level=4, fill=None, outline=None):

        x_size, y_size = image.size
        x_ll, y_ll, x_ur, y_ur = area_extent
        prj = pyproj.Proj(proj4_string)
        draw = ImageDraw.Draw(image)
        lon_ll, lat_ll = prj(x_ll, y_ll, inverse=True)
        lon_ur, lat_ur = prj(x_ur, y_ur, inverse=True)
        
        for i in range(level):
            #s = shapefile.Reader('/home/esn/data/gshhs/GSHHS_shp/%s/GSHHS_%s_L%s.shp' % 
            #                      (resolution, resolution, (i + 1)))
            shapefilename = os.path.join(self.db_root_path, 'GSHHS_shp', 
                                     resolution, 'GSHHS_%s_L%s.shp' % 
                                     (resolution, (i + 1)))
            s = shapefile.Reader(shapefilename)
            shapes = s.shapes()

            for j, shape in enumerate(shapes):
                s_lon_ll, s_lat_ll, s_lon_ur, s_lat_ur = shape.bbox
                if (lon_ur < s_lon_ll - 10 or lon_ll - 10 > s_lon_ur or 
                    lat_ur < s_lat_ll - 10 or lat_ll - 10 > s_lat_ur):
                    # Polygon is irrelevant
                    continue           
                shape_data = np.array(shape.points)
                lons = shape_data[:, 0]
                lats = shape_data[:, 1]

                x, y = prj(lons, lats)

                
                l_x = (x_ur - x_ll) / x_size
                l_y = (y_ur - y_ll) / y_size

                n_x = ((-x_ll + x) / l_x).astype(np.int)
                n_y = ((y_ur - y) / l_y).astype(np.int)

                index_array = np.vstack((n_x, n_y)).T

                
                draw.polygon(index_array.flatten().tolist(), fill=fill, outline=outline)
        

img = Image.open('BMNG_clouds_201109181715_areaT2.png')
proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
cw = ContourWriter('/home/esn/data/gshhs')
cw.add_coastlines(img, proj4_string, area_extent, resolution='l', level=4)
img.show()

