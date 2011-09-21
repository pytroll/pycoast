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
                       level=1, fill=None, outline=None):
        """
        resolution: {'c', 'l', 'i', 'h', 'f'}
        level: {1, 2, 3, 4}
        """

        draw = ImageDraw.Draw(image)
        
        # Area and projection info
        x_size, y_size = image.size
        x_ll, y_ll, x_ur, y_ur = area_extent
        prj = pyproj.Proj(proj4_string)
        
        # Calculate min and max lons and lats of interest
        # TODO: handle dateline and pole areas
         
        x_range = np.linspace(x_ll, x_ur, num=x_size)
        y_range = np.linspace(y_ll, y_ur, num=y_size)
        
        lons_s1, lats_s1 = prj(np.ones(y_range.size) * x_ll, y_range, inverse=True)
        lons_s2, lats_s2 = prj(x_range, np.ones(x_range.size) * y_ur, inverse=True)
        lons_s3, lats_s3 = prj(np.ones(y_range.size) * x_ur, y_range, inverse=True)
        lons_s4, lats_s4 = prj(x_range, np.ones(x_range.size) * y_ll, inverse=True)
        
        lon_min = lons_s1.min()
        lon_max = lons_s3.max()
        lat_min = lats_s4.min()
        lat_max = lats_s2.max()
        
        # Iterate through detail levels        
        for i in range(level):
            
            # One shapefile per level
            shapefilename = os.path.join(self.db_root_path, 'GSHHS_shp', 
                                     resolution, 'GSHHS_%s_L%s.shp' % 
                                     (resolution, (i + 1)))
            s = shapefile.Reader(shapefilename)
            shapes = s.shapes()

            # Iterate through shapes
            for j, shape in enumerate(shapes):
            
                # Check if polygon is possibly relevant
                s_lon_ll, s_lat_ll, s_lon_ur, s_lat_ur = shape.bbox
                if (lon_max < s_lon_ll or lon_min > s_lon_ur or 
                    lat_max < s_lat_ll or lat_min > s_lat_ur):
                    # Polygon is irrelevant
                    continue           
                
                # Get shape data as array and reproject    
                shape_data = np.array(shape.points)
                lons = shape_data[:, 0]
                lats = shape_data[:, 1]

                x, y = prj(lons, lats)

                # Convert to pixel index coordinates                
                l_x = (x_ur - x_ll) / x_size
                l_y = (y_ur - y_ll) / y_size

                n_x = ((-x_ll + x) / l_x).astype(np.int)
                n_y = ((y_ur - y) / l_y).astype(np.int)

                index_array = np.vstack((n_x, n_y)).T

                # Make PIL draw the polygon
                draw.polygon(index_array.flatten().tolist(), fill=fill, outline=outline)
      
    def add_coastlines_to_file(self, filename, proj4_string, area_extent, resolution='c', 
                               level=1, fill=None, outline=None):
        """
        resolution: {'c', 'l', 'i', 'h', 'f'}
        level: {1, 2, 3, 4}
        """
        
        image = Image.open(filename)
        self.add_coastlines(image, proj4_string, area_extent, resolution=resolution, 
                            level=level, fill=fill, outline=outline)
        image.save(filename)
        
if __name__ == '__main__':
    img = Image.open('BMNG_clouds_201109181715_areaT2.png')
    proj4_string = '+proj=stere +lon_0=8.00 +lat_0=50.00 +lat_ts=50.00 +ellps=WGS84'
    area_extent = (-3363403.31,-2291879.85,2630596.69,2203620.1)
    cw = ContourWriter('/home/esn/data/gshhs')
    cw.add_coastlines(img, proj4_string, area_extent, resolution='l', level=4)
    img.show()

