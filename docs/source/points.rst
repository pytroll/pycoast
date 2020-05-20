Custom points with or without textbox
--------------------------------

Pycoast can add a symbol to points of interest on an image. The following examples show how
we might use the :meth:`~pycoast.cw_agg.ContourWriterAGG.add_points`  method to annotate the
points on an image.

First of all, we setup a PIL image with an area definition, then we add coastlines and
borders for reference.

    >>> from PIL import Image
    >>> from pycoast import ContourWriterAGG
    >>> img = Image.new('RGB', (1024, 1024), (255, 255, 255))
    >>> proj4_string = '+proj=laea +lat_0=90 +lon_0=0 +a=6371228.0 +units=m'
    >>> area_extent = (-5326849.0625, -5326849.0625, 5326849.0625, 5326849.0625)
    >>> area_def = AreaDefinition('nh', 'nh', 'nh', proj4_string, 1024, 1024, area_extent)
    >>> cw = ContourWriterAGG('/home/esn/data/gshhs')
    >>> cw.add_coastlines(img, area_def, outline='black', resolution='l', level=4)
    >>> cw.add_borders(img, area_def, outline='black', width=3, level=1, resolution='c')

Now we can add a circle, which is the default symbol, with default point size 6 at the
location of Berlin, the name of the lacation will marked in a text box with black borders
and the default text size is 12.

    >>> points_list = [((13.4050, 52.5200), 'Berlin')]
    >>> cw.add_points(pil_img, area, points_list=points_list, font_file=font_file)

We can also annotate the image with text only by setting the ptsize to 0.
The example below will add 'Rome' at the given location without a symbol.

    >>> points_list = [((12.4964, 41.9028), 'Rome')]
    >>> cw.add_points(pil_img, area, points_list=points_list,
    ...               font_file=font_file, font_size=16,
    ...               symbol='circle', ptsize=0,
    ...               box_outline='black', text_linewidth=1,
    ...               box_fill='yellow', box_opacity=200)

Similarly, assign the description as an empty string will only draw the symbol on the image.
The example below will draw a square symbol at the location of Paris.

    >>> points_list = [((2.3522, 48.8566), '')]
    >>> cw.add_points(pil_img, area, points_list=points_list,i
    ...               font_file=font_file,
    ...               symbol='square', ptsize=10,
    ...               outline='red', width=1,
    ...               fill='blue', fill_opacity=128)

Finally, we can fully costomize the annotation as the example below, which will add
a circle in black with linewide set to 2 and filled in red color with opacity equals 255;
the description will be 'London' in a textbox with blue borders and filled with green color
with opacity set to 128.

    >>> points_list = [((0.1278, 51.5074), 'London')]
    >>> cw.add_points(img, area_def, points_list=points_list,
    ...               font_file=font_file, font_size=14,
    ...               symbol='circle', ptsize=14,
    ...               outline='black', width=2,
    ...               fill='red', fill_opacity=255,
    ...               box_outline='blue', box_linewidth=1.5,
    ...               box_fill='green', box_opacity=128)
    >>> img.show()

.. image:: images/nh_points_agg.png

The :meth:`~pycoast.cw_agg.ContourWriterAGG.add_points` method accepts a list of longitude, latitude pairs, with an optional
Description string.

.. _PIL: http://www.pythonware.com/products/pil/
