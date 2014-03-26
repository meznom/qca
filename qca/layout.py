import _qca
from collections import OrderedDict
import math

class Layout(object):
    def __init__(self):
        self._pl = _qca.Layout()

    @property
    def primitive_layout(self):
        return self._pl

    def __getstate__(self):
        i = OrderedDict()
        # This can be very verbose, disabled for now.
        # i['r_sites'] = self._pl.r_sites
        # i['r_charges'] = self._pl.r_charges
        # i['charges'] = self._pl.charges
        # i['epc'] = str(self._pl.epc)
        return i

    def __setstate__(self, i):
        # we do not deserialize properly
        # TODO: properly reconstruct state
        self.__init__()

    def coma_getstate(self):
        return self.__getstate__()

    def coma_setstate(self, i):
        self.__setstate__(i)

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

        assert len(self.boas) == self.N+1
        
        a = 1.0 / self.V1
        bs = [boa * a for boa in self.boas]
        x_off = 0
        self._pl.addDriverCell(x_off, 0, a, self.P1)
        x_off += bs[0]+a
        for b in bs[1:]:
            self._pl.addCell(x_off, 0, a)
            x_off += b+a
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

class InfiniteWire(Layout):
    def __init__(self, N_, N_dead_, V1_, boa_, P_):
        Layout.__init__(self)
        self.N = N_
        self.N_dead = N_dead_
        self.V1 = V1_
        self.boa = boa_
        self.P = P_
        self.construct_wire()
    
    def construct_wire(self):
        a = 1.0 / self.V1
        b = self.boa * a
        for i in range(0, self.N_dead):
            self._pl.addDriverCell(i*(b+a), 0, a, self.P)
        for i in range(self.N_dead, self.N_dead + self.N):
            self._pl.addCell(i*(b+a), 0, a)
        for i in range(self.N_dead + self.N, 2*self.N_dead + self.N):
            self._pl.addDriverCell(i*(b+a), 0, a, self.P)
    
    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'infinite_wire'
        i['N'] = self.N
        i['N_dead'] = self.N_dead
        i['V1'] = self.V1
        i['boa'] = self.boa
        i['P'] = self.P
        i.update(Layout.__getstate__(self))
        return i
    
    def __setstate__(self, i):
        self.__init__(i['N'], i['N_dead'], i['V1'], i['boa'], i['P'])

class AngleWire(Layout):
    def __init__(self, N_, V1_, doa_, theta_, P_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.doa = doa_
        self.theta = theta_
        self.P = P_

        a = 1.0 / self.V1
        d = self.doa * a
        t = self.theta*math.pi/180.0
        self._pl.addDriverCell(0, 0, a, self.P)
        for i in range(1,self.N+1):
            self._pl.addCell(d*math.cos(t)*i, d*math.sin(t)*i, a)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'angle_wire'
        i['N'] = self.N
        i['V1'] = self.V1
        i['doa'] = self.doa
        i['theta'] = self.theta
        i['P'] = self.P
        i.update(Layout.__getstate__(self))
        return i

    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['doa'], i['theta'], i['P'])

class KinkyWire(Layout):
    def __init__(self, N_, V1_, doa_, P_, kinks_):
        Layout.__init__(self)
        self.N = N_
        self.V1 = V1_
        self.doa = doa_
        self.P = P_
        self.kinks = kinks_
        
        a = 1.0 / self.V1
        d = self.doa * a
        theta,x,y = 0,0,0
        thetas,xs,ys = [],[],[]
        for i in range(self.N):
            if self.kinks.has_key(i):
                theta += self.kinks[i] * math.pi / 180.0
            x += d * math.cos(theta)
            y += d * math.sin(theta)
            thetas.append(theta)
            xs.append(x)
            ys.append(y)

        self._pl.addDriverCell(0, 0, a, self.P)
        for x,y in zip(xs,ys):
            self._pl.addCell(x, y, a)

    def __getstate__(self):
        i = OrderedDict()
        i['type'] = 'kinky_wire'
        i['N'] = self.N
        i['V1'] = self.V1
        i['doa'] = self.doa
        i['P'] = self.P
        i['kinks'] = self.kinks
        i.update(Layout.__getstate__(self))
        return i

    def __setstate__(self, i):
        self.__init__(i['N'], i['V1'], i['doa'], i['P'], i['kinks'])
