__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

from ._qca import QcaFixedCharge,QcaGrandCanonical,ElectronsPerCell
from .models import QcaBond
from .layout import Layout, Wire, NonuniformWire
