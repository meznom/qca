__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

from .models import QcaBond, QcaFixedCharge, QcaGrandCanonical
from .layout import Layout, Wire, NonuniformWire
from ._qca import ElectronsPerCell
