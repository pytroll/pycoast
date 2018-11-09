
# Import all the classes so that usage stays the same as previously.
# Default to PIL based ContourWriter in case the import of aggdraw or
# cairo fail.

from .version import __version__
from .cw_pil import ContourWriter as ContourWriterPIL
from .cw_agg import ContourWriterAGG
ContourWriter = ContourWriterAGG
