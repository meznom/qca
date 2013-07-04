import _qca
from . import __version__
from .layout import Layout
from collections import OrderedDict

class QcaError(Exception):
    pass

class QcaCommon(object):
    def __init__(self, qcaSystem):
        self.program = 'QcaEd'
        self.version = __version__
        self.results = OrderedDict()
        self.s = qcaSystem
        self.l = Layout()
        self.model = ''

    @property
    def t(self):
        return self.s.t

    @t.setter
    def t(self, t_):
        self.s.t = t_

    @property
    def V0(self):
        return self.s.V0

    @V0.setter
    def V0(self, V0_):
        self.s.V0 = V0_

    @property
    def mu(self):
        return self.s.mu

    @mu.setter
    def mu(self, mu_):
        self.s.mu = mu_

    @property
    def l(self):
        return self._l

    @l.setter
    def l(self, l_):
        self._l = l_
        self.s.l = l_.primitive_layout

    @property
    def beta(self):
        return self.s.beta

    @beta.setter
    def beta(self, beta_):
        self.s.beta = beta_

    @property
    def T(self):
        return 1.0/self.beta

    @T.setter
    def T(self, T_):
        self.beta = 1.0 / T_

    @property
    def N_p(self):
        return self.s.N_p

    @property
    def N_sites(self):
        return self.s.N_sites

    def init(self):
        pass

    def energies(self):
        return self.s.energies()

    def run(self, changedTonly=False):
        if not changedTonly:
            self.s.update()
        self.results['P'] = [self.s.measurePolarization(i) for i in range(self.N_p)]
        self.results['N'] = [self.s.measureParticleNumber(i) for i in range(self.N_p)]

    def __getstate__(self):
        i = OrderedDict()
        i['info'] = OrderedDict()
        i['info']['program'] = self.program
        i['info']['version'] = self.version
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

    def __setstate__(self, i):
        # relies on derived classes calling __init__ first
        if (self.program != i['info']['program'] or
            self.version != i['info']['version'] or
            self.model != i['parameters']['model']):
            raise QcaError('QCA program versions or QCA models do not agree')
        self.t = i['parameters']['t']
        self.V0 = i['parameters']['V0']
        self.mu = i['parameters']['mu']
        self.T = i['parameters']['T']
        self.l = i['parameters']['layout']
        self.results = i['results']

    def __eq__(self, s):
        d1 = self.__dict__.copy()
        d2 = s.__dict__.copy()
        del d1['s']
        del d2['s']
        return d1 == d2

class QcaBond(QcaCommon):
    def __init__(self):
        QcaCommon.__init__(self, _qca.QcaBond())
        self.model = 'QcaBond'

    def __setstate__(self, i):
        self.__init__()
        QcaCommon.__setstate__(self, i)

class QcaFixedCharge(QcaCommon):
    def __init__(self):
        QcaCommon.__init__(self, _qca.QcaFixedCharge())
        self.model = 'QcaFixedCharge'

    def __setstate__(self, i):
        self.__init__()
        QcaCommon.__setstate__(self, i)

class QcaGrandCanonical(QcaCommon):
    def __init__(self):
        QcaCommon.__init__(self, _qca.QcaGrandCanonical())
        self.model = 'QcaGrandCanonical'
    
    def __setstate__(self, i):
        self.__init__()
        QcaCommon.__setstate__(self, i)
