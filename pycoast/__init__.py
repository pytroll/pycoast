
# Import all the classes so that usage stays the same as previously.

from .cw_pil import ContourWriter as ContourWriterPIL
from .cw_agg import ContourWriterAGG
ContourWriter = ContourWriterAGG

from .version import get_versions
__version__ = get_versions()['version']
del get_versions
