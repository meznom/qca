__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

from .models import QcaBond, QcaFixedCharge, QcaGrandCanonical, QcaError
from .layout import Layout, Wire, NonuniformWire, WireWithTwoDriverCells, \
                    NonuniformWireWithTwoDriverCells, InfiniteWire
from .selfconsistency import SelfConsistency
from ._qca import ElectronsPerCell