#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Pycoast package for adding geographic-based decorations to images."""
from pycoast.cw_base import get_resolution_from_area  # noqa

from .cw_agg import ContourWriterAGG  # noqa
from .cw_pil import ContourWriterPIL
from .version import get_versions

__version__ = get_versions()["version"]
del get_versions


class ContourWriter(ContourWriterPIL):
    """Writer wrapper for deprecation warning.

    .. deprecated:: 1.2.0

        Use :class:`~pycoast.cw_pil.ContourWriterPIL` or :class:`~pycoast.cw_agg.ContourWriterAGG` instead.

    """

    def __init__(self, *args, **kwargs):
        import warnings

        warnings.warn(
            "'ContourWriter' has been deprecated please use "
            "'ContourWriterPIL' or 'ContourWriterAGG' instead",
            DeprecationWarning,
        )
        super(ContourWriter, self).__init__(*args, **kwargs)
