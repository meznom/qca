import unittest
import qca
from coma import Measurement

class TestQCA(unittest.TestCase):
    def test_run(self):
        s = qca.QcaBond()

        s.l.wire(2,1,3,1)

        s.V0 = 1E6
        s.beta = 1

        s.update()

        P_1 = s.measurePolarization(0)
        P_2 = s.measurePolarization(1)

        print('Polarizations: {}, {}'.format(P_1, P_2))

    def test_serialization(self):
        s = qca.QcaBond()
        m = Measurement(1, '.')
        m.start()
        m.end()
        m.save(s)
        print(m)
