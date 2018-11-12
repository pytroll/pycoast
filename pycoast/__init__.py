#!/usr/bin/env python
# -*- coding: utf-8 -*-
from .version import get_versions
__version__ = get_versions()['version']
del get_versions

from .cw_pil import ContourWriterPIL
from .cw_agg import ContourWriterAGG


class ContourWriter(ContourWriterPIL):
    def __init__(self, *args, **kwargs):
        import warnings
        warnings.warn("'ContourWriter' has been deprecated please use "
                      "'ContourWriterPIL' or 'ContourWriterAGG' instead", DeprecationWarning)
        super(ContourWriter, self).__init__(*args, **kwargs)