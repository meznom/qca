import _qca
from collections import OrderedDict

class Layout(object):
    def __init__(self):
        self._pl = _qca.PrimitiveLayout()

    @property
    def primitive_layout(self):
        return self._pl

    def __getstate__(self):
        i = OrderedDict()
        i['r_sites'] = self._pl.r_sites
        i['r_charges'] = self._pl.r_charges
        i['charges'] = self._pl.charges
        i['epc'] = self._pl.epc
        return i

class Wire(Layout):
    def __init__(self, N_, V1_, boa_, P_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.boa = boa_
        self.P = P_
        self._pl.wire(N_, 1.0/V1_, boa_ * 1.0 / V1_, P_)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'wire'
        i['N'] = self.N
        i['V1'] = self.V1
        i['boa'] = self.boa
        i['P'] = self.P
        i.update(Layout.__getstate__(self))
        return i

class NonuniformWire(Layout):
    def __init__(self, N_, V1_, boas_, P_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.boas = boas_
        self.P = P_
        a = 1.0 / V1_
        bs = [boa * a for boa in boas_]
        self._pl.nonuniformWire(N_, a, bs, P_)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'nonuniformwire'
        i['N'] = self.N
        i['V1'] = self.V1
        i['boas'] = self.boas
        i['P'] = self.P
        i.update(Layout.__getstate__(self))
        return i
