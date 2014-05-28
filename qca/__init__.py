# Copyright (c) 2013,2014 Burkhard Ritter
# This code is distributed under the two-clause BSD License.

__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

from .models import QcaBond, QcaFixedCharge, QcaGrandCanonical, QcaIsing, QcaError
from .layout import Layout, Wire, NonuniformWire, WireWithTwoDriverCells, \
                    NonuniformWireWithTwoDriverCells, InfiniteWire, AngleWire, \
                    KinkyWire, MajorityGate
from .selfconsistency import SelfConsistency
from ._qca import ElectronsPerCell
from . import test
