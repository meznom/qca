__version__ = 'unknown'

try:
    from ._version import __version__
except ImportError:
    pass

#from ._qca import Layout,QcaBond,QcaFixedCharge,QcaGrandCanonical,ElectronsPerCell
from ._qca import Layout,QcaFixedCharge,QcaGrandCanonical,ElectronsPerCell
from .models import QcaBond
