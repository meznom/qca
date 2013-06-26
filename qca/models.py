import _qca
from . import __version__
from .layout import Layout
from collections import OrderedDict

class QcaBond(_qca.QcaBond):
    def __init__(self):
        _qca.QcaBond.__init__(self)
        self.l = Layout()
        self.program = 'QcaEd'
        self.version = __version__
        self.results = OrderedDict()

    def measureAll(self):
        self.results['P'] = [self.measurePolarization(i) for i in range(self.N_p())]
        self.results['N'] = [self.measureParticleNumber(i) for i in range(self.N_p())]

    def __getstate__(self):
        i = OrderedDict()
        i['parameters'] = OrderedDict()
        i['parameters']['model'] = 'QcaBond'
        i['parameters']['t'] = self.t
        i['parameters']['V0'] = self.V0
        i['parameters']['mu'] = self.mu
        i['parameters']['beta'] = self.beta
        i['parameters']['layout'] = self.l
        # ignore td, Vext, epsilonr, lambdaD, epsilon0, q
        i['results'] = self.results
        return i
