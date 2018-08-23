
# Import all the classes so that usage stays the same as previously.
# Default to PIL based ContourWriter in case the import of aggdraw or
# cairo fail.

from pycoast.version import __version__
from pycoast.cw_pil import ContourWriter
try:
    from pycoast.cw_agg import ContourWriterAGG
    from pycoast.decorator_agg import DecoratorAGG
except ImportError:
    # FIXME: This seems really wrong
    ContourWriterAGG = ContourWriter

