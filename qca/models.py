import _qca
from . import __version__
from .layout import Layout
from collections import OrderedDict

class QcaCommon(object):
    def __init__(self):
        self.program = 'QcaEd'
        self.version = __version__
        self.results = OrderedDict()
        self.l = Layout()
        self.model = ''

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l):
        self._l = l
        self._primitive_layout = l.primitive_layout

    @property
    def T(self):
        return 1.0/self.beta

    @T.setter
    def T(self, T_):
        self.beta = 1.0 / T_

    def init(self):
        pass

    def run(self, changedTonly=False):
        if not changedTonly:
            self.update()
        self.results['P'] = [self.measurePolarization(i) for i in range(self.N_p)]
        self.results['N'] = [self.measureParticleNumber(i) for i in range(self.N_p)]

    def __getstate__(self):
        i = OrderedDict()
        i['parameters'] = OrderedDict()
        i['parameters']['model'] = self.model
        i['parameters']['t'] = self.t
        i['parameters']['V0'] = self.V0
        i['parameters']['mu'] = self.mu
        i['parameters']['T'] = self.T
        i['parameters']['layout'] = self.l
        # ignore td, Vext, epsilonr, lambdaD, epsilon0, q
        i['results'] = self.results
        return i

class QcaBond(_qca.QcaBond, QcaCommon):
    def __init__(self):
        _qca.QcaBond.__init__(self)
        QcaCommon.__init__(self)
        self.model = 'QcaBond'

class QcaFixedCharge(_qca.QcaFixedCharge, QcaCommon):
    def __init__(self):
        _qca.QcaFixedCharge.__init__(self)
        QcaCommon.__init__(self)
        self.model = 'QcaFixedCharge'

class QcaGrandCanonical(_qca.QcaGrandCanonical, QcaCommon):
    def __init__(self):
        _qca.QcaGrandCanonical.__init__(self)
        QcaCommon.__init__(self)
        self.model = 'QcaGrandCanonical'
