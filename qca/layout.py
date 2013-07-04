import _qca
from collections import OrderedDict

class Layout(object):
    def __init__(self):
        self._pl = _qca.Layout()

    @property
    def primitive_layout(self):
        return self._pl

    def __getstate__(self):
        i = OrderedDict()
        i['r_sites'] = self._pl.r_sites
        i['r_charges'] = self._pl.r_charges
        i['charges'] = self._pl.charges
        i['epc'] = str(self._pl.epc)
        return i

    def __eq__(self, l):
        d1 = self.__dict__.copy()
        d2 = l.__dict__.copy()
        pl1 = d1['_pl']
        pl2 = d2['_pl']
        del d1['_pl']
        del d2['_pl']
        return (d1 == d2 and
                pl1.r_sites == pl2.r_sites and
                pl1.r_charges == pl2.r_charges and
                pl1.charges == pl2.charges and
                pl1.epc == pl2.epc)

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

    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['boa'], i['P'])

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

    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['boas'], i['P'])

class WireWithTwoDriverCells(Layout):
    def __init__(self, N_, V1_, boa_, P1_, P2_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.boa = boa_
        self.P1 = P1_
        self.P2 = P2_

        a = 1.0 / self.V1
        b = self.boa * a
        for i in range(self.N):
            self._pl.addCell((a+b)*i, 0, a)
        self._pl.addDriverCell(-b-a, 0, a, self.P1)
        self._pl.addDriverCell(self.N*(b+a), 0, a, self.P2)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'wire_with_two_driver_cells'
        i['N'] = self.N
        i['V1'] = self.V1
        i['boa'] = self.boa
        i['P1'] = self.P1
        i['P2'] = self.P2
        i.update(Layout.__getstate__(self))
        return i
    
    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['boa'], i['P1'], i['P2'])

class NonuniformWireWithTwoDriverCells(Layout):
    def __init__(self, N_, V1_, boas_, P1_, P2_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.boas = boas_
        self.P1 = P1_
        self.P2 = P2_

        assert len(self.boas) == self.N
        
        a = 1.0 / self.V1
        bs = [boa * a for boa in self.boas]
        x_off = 0
        self._pl.addDriverCell(x_off, 0, a, self.P1)
        for b in bs:
            x_off += b+a
            self._pl.addCell(x_off, 0, a)
        self._pl.addDriverCell(x_off, 0, a, self.P2)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'nonuniformwire_with_two_driver_cells'
        i['N'] = self.N
        i['V1'] = self.V1
        i['boas'] = self.boas
        i['P1'] = self.P1
        i['P2'] = self.P2
        i.update(Layout.__getstate__(self))
        return i
    
    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['boas'], i['P1'], i['P2'])
