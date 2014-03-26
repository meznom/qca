__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

from .models import QcaBond, QcaFixedCharge, QcaGrandCanonical, QcaIsing, QcaError
from .layout import Layout, Wire, NonuniformWire, WireWithTwoDriverCells, \
                    NonuniformWireWithTwoDriverCells, InfiniteWire, AngleWire, \
                    KinkyWire
from .selfconsistency import SelfConsistency
from ._qca import ElectronsPerCell
from . import test
