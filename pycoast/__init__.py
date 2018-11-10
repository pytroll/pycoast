#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .version import get_versions
__version__ = get_versions()['version']
del get_versions

from .cw_pil import ContourWriterPIL
from .cw_agg import ContourWriterAGG
ContourWriter = ContourWriterAGG

