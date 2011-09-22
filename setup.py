from setuptools import setup

import imp

version = imp.load_source('pycoast.version', 'pycoast/version.py')

setup(name='pycoast',
      version=version.__version__,
      description='Writing of coastlines, borders and rivers to images in Python',
      author='Esben S. Nielsen',
      author_email='esn@dmi.dk',
      packages = ['pycoast'],      
      install_requires=['PIL', 'pyshp'], 
      zip_safe = False
      )

