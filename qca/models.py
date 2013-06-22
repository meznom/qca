import _qca
from . import __version__
from collections import OrderedDict

class QcaBond(_qca.QcaBond):
    def __init__(self):
        _qca.QcaBond.__init__(self)
        self.program = 'QcaEd'
        self.version = __version__

    def __getstate__(self):
        i = OrderedDict()
        i['parameters'] = OrderedDict()
        i['parameters']['model'] = 'QcaBond'
        i['parameters']['t'] = self.t
        i['parameters']['V0'] = self.V0
        i['parameters']['beta'] = self.beta
        return i
