
# Import all the classes so that usage stays the same as previously.
# Default to PIL based ContourWriter in case the import of aggdraw or
# cairo fail.

from .cw_pil import ContourWriter
try:
    from .cw_agg import ContourWriterAGG
except ImportError:
    ContourWriterAGG = ContourWriter
try:
    from .cw_cairo import ContourWriterCairo
except ImportError:
    ContourWriterCairo = ContourWriterAGG
