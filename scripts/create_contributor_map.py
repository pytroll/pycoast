#!/usr/bin/env python3
# encoding: utf-8
"""Create a PNG world map with PyTroll contributor countries highlighted.

Some files may be downloaded or created to create the image. It is suggested
that this script be run from the directory it is stored in. Example:

    python create_contributor_map.py -u davidh-ssec

"""
import getpass
import logging
import os
import sys

import json
import geocoder
import numpy as np
from github import Github
from PIL import Image
from pyproj import Proj
from time import sleep
import pycountry
from collections import Counter

from pycoast import ContourWriter, ContourWriterAGG

LOG = logging.getLogger(__name__)
geocoders = [geocoder.google, geocoder.arcgis_reverse.ArcgisReverse]
MAX_CONTRIBUTORS = 10.


def get_all_pytroll_contributors(user):
    g = Github(user, getpass.getpass('GitHub Password: '))
    pytroll_org = g.get_organization('pytroll')
    LOG.info("Getting all PyTroll repositories...")
    pytroll_repos = [x for x in pytroll_org.get_repos() if x.name != u'aggdraw']
    LOG.info("Getting all PyTroll contributors...")
    all_pytroll_contributors = [
        u for r in pytroll_repos for u in r.get_contributors()]
    set_pytroll_contributors = {u.login: u for u in all_pytroll_contributors}
    return set_pytroll_contributors


def get_country(locations, cache_file='countries.json'):
    """Return 3-digit country code for each location provided.
    """
    # cache country results because most APIs have query limits
    if os.path.isfile(cache_file):
        json_cache = json.load(open(cache_file, 'r'))
    else:
        json_cache = {}

    need_sleep = False
    for loc in locations:
        for gcoder in geocoders:
            if loc in json_cache:
                yield json_cache[loc]
                need_sleep = False
                break
            else:
                need_sleep = True
                LOG.debug("Finding: %s", loc)
                result = gcoder(loc)
                if not result.ok or result.country is None:
                    LOG.debug("Bad result using %r, country %r", gcoder, result.country)
                    continue
                country = result.country
                if len(country) == 2:
                    country = pycountry.countries.get(alpha_2=country).alpha_3
                else:
                    try:
                        country = pycountry.countries.get(alpha_3=country).alpha_3
                    except KeyError:
                        LOG.error("Invalid country code: %s", country)
                        continue
                if country is not None:
                    json_cache[loc] = country
                yield country
                break
        else:
            try:
                # last resort
                country = pycountry.countries.get(name=loc).alpha_3
                json_cache[loc] = country
                yield country
                continue
            except KeyError:
                pass
            LOG.warning("Could not find country for {}".format(loc))
            yield None
        if need_sleep:
            sleep(1.0)

    json.dump(json_cache, open(cache_file, 'w'), indent=4, sort_keys=True)


def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Create a world map image of PyTroll contributor locations")
    parser.add_argument('--map-proj', default="+proj=moll",
                        help='PROJ.4 map projection string to use for the output image')
    parser.add_argument('--output-size', default=(1200, 600), nargs=2, type=int,
                        help='\'width height\' of output image')
    parser.add_argument('--contributor-color', default=(255, 127, 127),
                        help='Color of dots for each contributor')
    parser.add_argument('--bg-color', default=None, nargs=3, type=int,
                        help='Background color of image \'R G B\' 0-255 (default: transparent)')
    parser.add_argument('-o', '--output', default='contributors_map.png',
                        help='Output filename to save the image to')
    parser.add_argument('-u', '--github-user', required=True,
                        help='github username')
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO)

    # Check if we have the shapefile we want
    if not os.path.isfile('ne_110m_admin_0_countries.shp'):
        # Download the NaturalEarth shapefile we are going to use
        import urllib.request
        import shutil
        import zipfile
        url = "http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/110m/cultural/ne_110m_admin_0_countries.zip"
        zip_fn = 'ne_110m_admin_0_countries.zip'
        LOG.info("Downloading NaturalEarth Shapefile")
        with urllib.request.urlopen(url) as response, open(zip_fn, 'wb') as out_file:
            shutil.copyfileobj(response, out_file)
        zip_ref = zipfile.ZipFile(zip_fn, 'r')
        zip_ref.extractall('.')
        zip_ref.close()

    im = Image.new('RGBA', args.output_size, color=tuple(args.bg_color) if args.bg_color is not None else None)
    cw = ContourWriterAGG()
    p = Proj(args.map_proj)
    a, _ = p(-180, 0)
    _, b = p(0, -90)
    c, _ = p(180, 0)
    _, d = p(0, 90)

    extents = (a, b, c, d)
    area_def = (args.map_proj, extents)

    LOG.info("Getting PyTroll contributors...")
    contributors = get_all_pytroll_contributors(args.github_user)
    LOG.info("Getting all PyTroll contributor locations...")
    all_user_locations = []
    for user in contributors.values():
        if user.location is None or user.location == '':
            LOG.info("User location not specified:  ID: '{}';\tName: '{}';\tURL: '{}'".format(user.id, user.name, user.url))
            continue
        all_user_locations.append(user.location)
    all_countries = list(get_country(all_user_locations))
    number_contributors = Counter(all_countries)
    del number_contributors[None]
    all_countries = list(number_contributors.keys())

    def shape_generator():
        import shapefile
        reader = shapefile.Reader("ne_110m_admin_0_countries.shp")
        name_idx = [idx for idx, f in enumerate(reader.fields) if f[0] == 'name'][0]
        for sr in reader.shapeRecords():
            name = sr.record[name_idx]
            if not isinstance(name, str):
                try:
                    name = name.decode()
                except UnicodeDecodeError:
                    name = name.decode('latin-1')
            shape_country = list(get_country([name]))[0]

            if shape_country is not None and shape_country in all_countries:
                LOG.info("Contributor country: %s, Code: %s, Number of Contributors: %d", name, shape_country, number_contributors[shape_country])
                num_cont = min(number_contributors[shape_country], MAX_CONTRIBUTORS)
                kwargs = {
                    'fill': args.contributor_color,
                    'fill_opacity': int((num_cont / 10.) * 55 + 200)
                }
                yield sr.shape, kwargs
            else:
                LOG.debug("No contributor's from Country: %s, Code: %s", name, shape_country)
                yield sr.shape, None

    LOG.info("Applying contributor locations to map image...")
    cw.add_shapes(im, area_def, 'polygon', shape_generator(),
                  outline=(0, 0, 0), outline_opacity=255, width=1,
                  fill=(127, 127, 127), fill_opacity=127)
    im.save(args.output, dpi=(300, 300))


if __name__ == "__main__":
    sys.exit(main())
