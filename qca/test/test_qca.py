# Copyright (c) 2013,2014 Burkhard Ritter
# This code is distributed under the two-clause BSD License.

import unittest
import os
import shutil
import pickle
import qca
from coma import FileMeasurement, Experiment
import numpy as np
import matplotlib.pyplot as plt

class TestQCA(unittest.TestCase):
    def setUp(self):
        self.f = os.path.join(os.path.dirname(__file__), 'testmeasurement.xml')
        self.d = os.path.join(os.path.dirname(__file__), 'test_create_an_experiment')
        if True:
            if os.path.exists(self.f):
                os.remove(self.f)
            if os.path.exists(self.d):
                shutil.rmtree(self.d)

    def tearDown(self):
        if True:
            if os.path.exists(self.f):
                os.remove(self.f)
            if os.path.exists(self.d):
                shutil.rmtree(self.d)

    def test_polarization_for_different_layouts(self):
        # these tests come from qcaTest.cpp, 
        # "reuse_same_system_multiple_times_with_different_layouts"
        s = qca.QcaBond()
        s.l = qca.Wire(2,100,2,1)
        s.beta = 1
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.863919, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.603095, 5)

        # for changed beta, update() is not necessary
        s.beta = 2
        s.init()
        s.run(changedTonly=True)
        self.assertAlmostEqual(s.results['P'][0], 0.986214, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.813721, 5)

        # change layout, for example, inter-cell spacing
        s.l = qca.Wire(2, 50, 1, 1)
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.978439, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.00399499, 5)

        # change layout completely: different number of cells
        s.l = qca.Wire(3, 1.0/0.03, 2, 1);
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.668657, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.455502, 5)
        self.assertAlmostEqual(s.results['P'][2], 0.115793, 5)
        
        # again, change layout completely
        s.l = qca.Wire(1, 100, 3, 0.4);
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.234819, 5)

        # test fixed charge system as well
        s = qca.QcaFixedCharge()
        s.l = qca.Wire(2, 100, 2, 0.1)
        s.V0 = 1000
        s.beta = 100
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.940904, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.765352, 5)

        # now from "limits_of_nonuniform_wire"
        # A non-uniform with uniform inter-cell spacing is just a uniform wire
        s1 = qca.QcaBond()
        s1.l = qca.NonuniformWire(3, 100, [2,2,2], 1)
        s1.beta = 1
        s1.init()
        s1.run()
        
        s2 = qca.QcaBond()
        s2.l = qca.Wire(3, 100, 2, 1);
        s2.beta = 1
        s2.init()
        s2.run()

        self.assertAlmostEqual(s1.results['P'][0], s2.results['P'][0], 5)
        self.assertAlmostEqual(s1.results['P'][1], s2.results['P'][1], 5)
        self.assertAlmostEqual(s1.results['P'][2], s2.results['P'][2], 5)

        # In the limit of a very large inter-cell spacing, we should see no
        # polarization response. First we put the driver cell far away.
        s = qca.QcaBond()
        s.l = qca.NonuniformWire(3, 100, [1000, 2, 2], 1)
        s.beta = 1
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0, 5)
        self.assertAlmostEqual(s.results['P'][1], 0, 5)
        self.assertAlmostEqual(s.results['P'][2], 0, 5)

        # Now we move the right-most cell (the output cell) far away.
        s.l = qca.NonuniformWire(3, 100, [2,2,1000], 1)
        s.beta = 1
        s.init()
        s.run()
        self.assertGreater(s.results['P'][0], 0.3)
        self.assertGreater(s.results['P'][1], 0.3)
        self.assertAlmostEqual(s.results['P'][2], 0, 5)

        # from "simple_sanity_checks_for_qca_grand_canonical_system"
        s = qca.QcaGrandCanonical()
        s.l = qca.Wire(1, 160, 3, 0)
        s.V0 = 1000
        s.mu = 300
        s.beta = 1000000
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0, 10)

        s.l = qca.Wire(1, 160, 3, 0.1)
        s.beta = 100000
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['N'][0][4], 2, 10)
        self.assertGreater(s.results['P'][0], 0.1)
        self.assertLess(s.results['P'][0], 1)

    def test_retrieve_energies(self):
        s = qca.QcaBond()
        s.l = qca.Wire(1,1,2,0)
        s.init()
        s.run()
        es = s.energies()
        self.assertEqual(len(es), 6)
        
        s.l = qca.Wire(2,1,2,0)
        s.init()
        s.run()
        es = s.energies()
        self.assertEqual(len(es), 36)

        s = qca.QcaFixedCharge()
        s.l = qca.Wire(1,1,2,0)
        s.init()
        s.run()
        es = s.energies()
        self.assertEqual(len(es), 28)
        
        s.l = qca.Wire(2,1,2,0)
        s.init()
        s.run()
        es = s.energies()
        self.assertEqual(len(es), 28*28)

        s = qca.QcaGrandCanonical()
        s.l = qca.Wire(1,1,2,0)
        s.init()
        s.run()
        es = s.energies()
        self.assertEqual(len(es), 256)

    def test_compensation_charge(self):
        s = qca.QcaBond()
        s.l = qca.Wire(1,40,2,1)
        s.T = 1
        s.V0 = 1E6
        
        Ps = []
        s.q = 0
        s.init()
        s.run()
        Ps.append(s.results['P'][0])

        s.q = 0.25
        s.init()
        s.run()
        Ps.append(s.results['P'][0])

        s.q = 0.5
        s.init()
        s.run()
        Ps.append(s.results['P'][0])

        print(Ps)

        self.assertGreater(Ps[2], Ps[1])
        self.assertGreater(Ps[1], Ps[0])

    def test_pickling_and_unpickling(self):
        s = qca.QcaBond()
        s.t = 2
        s.V0 = 1E6
        s.mu = 3
        s.T = 0.1
        s.q = 0.1
        s.l = qca.NonuniformWire(3,100,[2.0,2.1,2.2],0.7)
        s.init()
        s.run()

        ss = pickle.dumps(s)
        s2 = pickle.loads(ss)
        self.assertEqual(s, s2)

        s2.init()
        s2.run()
        self.assertEqual(s,s2)

        with self.assertRaises(qca.QcaError):
            s.version = 'foo'
            ss = pickle.dumps(s)
            s2 = pickle.loads(ss)

    def test_demonstrate_recommended_usage(self):
        s = qca.QcaBond()
        s.l = qca.Wire(2,100,2.5,1)
        s.t = 1
        s.V0 = 1E6
        s.mu = 0 # does not do anything for bond model
        s.T = 1

        s.init()
        s.run()
        print('T = {}\nP_1 = {}\nP_2 = {}\n'.format(s.T, s.results['P'][0], s.results['P'][1]))

        s.T = 0.1
        s.init()
        s.run(changedTonly=True)
        print('T = {}\nP_1 = {}\nP_2 = {}\n'.format(s.T, s.results['P'][0], s.results['P'][1]))
    
    def test_create_one_measurement(self):
        s = qca.QcaBond()
        s.l = qca.NonuniformWire(2,100,[2.5,3],1)
        s.V0 = 1E6
        s.beta = 1
        
        m = FileMeasurement(self.f)
        m.start()
        
        s.init()
        s.run()
        
        m.end()
        m.save(s)

    def test_create_an_experiment_and_create_a_plot(self):
        e = Experiment(self.d,1)
        e.description = 'Test creating an experiment for the qca module.'
        e.reset()
        e.start()

        s1 = qca.QcaBond()
        s1.l = qca.Wire(2, 100, 2, 1)
        s1.V0 = 1E6
        s1.T = 1
        s1.init()
        s1.run()
        
        s2 = qca.QcaFixedCharge()
        s2.l = qca.Wire(2, 100, 2, 1)
        s2.V0 = 1E6
        s2.T = 1
        s2.init()
        s2.run()
        
        for T in np.linspace(0.0001, 4, 200):
            m = e.new_measurement()
            m.start()
            s1.T = T
            s1.init()
            s1.run(changedTonly=True)
            m.end()
            m.save(s1)

            m = e.new_measurement()
            m.start()
            s2.T = T
            s2.init()
            s2.run(changedTonly=True)
            m.end()
            m.save(s2)

        e.end()
        e.save()

        Ts_b,Ts_f,Ps_b,Ps_f = [],[],[],[]
        for m in e.measurements():
            if m['parameters/model'] == 'QcaBond':
                Ts_b.append(m['parameters/T'])
                Ps_b.append(m['results/P'][1])
            elif m['parameters/model'] == 'QcaFixedCharge':
                Ts_f.append(m['parameters/T'])
                Ps_f.append(m['results/P'][1])
        # print(Ts_b,Ps_b)
        # print(Ts_f,Ps_f)
        
        fig = plt.figure(figsize=(10,10))
        p = fig.add_subplot(1,1,1)
        p.plot(Ts_f,Ps_f,label='Fixed')
        p.plot(Ts_b,Ps_b,label='Bond')
        p.legend()
        p.set_xlabel('temperature')
        p.set_ylabel('polarization $P_2$')
        p.text(0.98,0.80, '$V_1 = {}$'.format(s1.l.V1),
                fontsize='12',
                horizontalalignment='right', 
                verticalalignment='top', 
                transform=p.transAxes)
        fig.savefig(os.path.join(self.d,'bond_and_fixedcharge.pdf'))

    def test_can_construct_qca_ising(self):
        s = qca.QcaIsing()
        s.t = 1
        s.l = qca.Wire(3,100,2,1)
        s.beta = 1
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][0], 0.846923, 5)
        self.assertAlmostEqual(s.results['P'][1], 0.747275, 5)
        self.assertAlmostEqual(s.results['P'][2], 0.653792, 5)
        self.assertNotEqual(s.tprime, 0)
        self.assertAlmostEqual(s.tprime, 0.13656854249492380195, 8)

    def test_angle_wire(self):
        s = qca.QcaBond()

        s.l = qca.AngleWire(3, 100, 4.0, 0, 1)
        s.V0 = 1E6
        s.T = 1
        s.init()
        s.run()
        r1 = s.results.copy()

        s.l = qca.AngleWire(3, 100, 4.0, 90, 1)
        s.V0 = 1E6
        s.T = 1
        s.init()
        s.run()
        r2 = s.results.copy()

        s.l = qca.Wire(3, 100, 3.0, 1)
        s.V0 = 1E6
        s.T = 1
        s.init()
        s.run()
        r3 = s.results.copy()

        print('\nP_1 = {}\nP_2 = {}\nP_3 = {}'.format(r1['P'][0], r1['P'][1], r1['P'][2]))

        self.assertAlmostEqual(r1['P'][0], r2['P'][0])
        self.assertAlmostEqual(r1['P'][1], r2['P'][1])
        self.assertAlmostEqual(r1['P'][2], r2['P'][2])

        self.assertAlmostEqual(r1['P'][0], r3['P'][0])
        self.assertAlmostEqual(r1['P'][1], r3['P'][1])
        self.assertAlmostEqual(r1['P'][2], r3['P'][2])
    
    def test_kinky_wire(self):
        s = qca.QcaBond()

        # Compare normal wire and kinky wire without kinks
        s.l = qca.Wire(3,100,3.0,1)
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        r_ = s.results.copy()
        print('\nKinky wire, equal to a normal wire')
        print('P_1 = {}\nP_2 = {}\nP_3 = {}'
              .format(r_['P'][0], r_['P'][1], r_['P'][2]))

        kinks = [{},{0:0},{1:0},{1:0,2:0},{0:90},{3:10},{5:90}]
        for k in kinks:
            s.l = qca.KinkyWire(3,100,4.0,1,k)
            s.V0 = 1E6
            s.T = 1
            s.q = 0.5
            s.init()
            s.run()
            r = s.results.copy()
            self.assertAlmostEqual(r_['P'][0], r['P'][0])
            self.assertAlmostEqual(r_['P'][1], r['P'][1])
            self.assertAlmostEqual(r_['P'][2], r['P'][2])


        # Compare angle wire and kinky wire without kinks
        s.l = qca.AngleWire(3,100,4.0,45.0,1)
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        r_ = s.results.copy()
        print('\nKinky wire, equal to a 45 degree wire')
        print('P_1 = {}\nP_2 = {}\nP_3 = {}'
              .format(r_['P'][0], r_['P'][1], r_['P'][2]))

        kinks = [((0,45),),((0,45),(2,0)),((0,-45),),((0,135),),((0,45),(3,10))]
        for k in kinks:
            s.l = qca.KinkyWire(3,100,4.0,1,k)
            s.V0 = 1E6
            s.T = 1
            s.q = 0.5
            s.init()
            s.run()
            r = s.results.copy()
            self.assertAlmostEqual(r_['P'][0], r['P'][0])
            self.assertAlmostEqual(r_['P'][1], r['P'][1])
            self.assertAlmostEqual(r_['P'][2], r['P'][2])

        # Compare normal wire and kinky wire with one kink
        s.l = qca.Wire(3,100,3.0,1)
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        r_ = s.results.copy()
        print('\nNormal wire')
        print('P_1 = {}\nP_2 = {}\nP_3 = {}'
              .format(r_['P'][0], r_['P'][1], r_['P'][2]))

        s.l = qca.KinkyWire(3,100,4.0,1,{1: 10})
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        r = s.results.copy()
        print('\nKinky wire with 10 degree kink after the first cell')
        print('P_1 = {}\nP_2 = {}\nP_3 = {}'
              .format(r['P'][0], r['P'][1], r['P'][2]))
        self.assertAlmostEqual(r_['P'][0], r['P'][0], 2)
        self.assertGreater(r_['P'][1], r['P'][1])
        self.assertGreater(r_['P'][2], r['P'][2])

        s.l = qca.KinkyWire(3,100,4.0,1,{2: -10})
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        r = s.results.copy()
        print('\nKinky wire with -10 degree kink after the second cell')
        print('P_1 = {}\nP_2 = {}\nP_3 = {}'
              .format(r['P'][0], r['P'][1], r['P'][2]))
        self.assertAlmostEqual(r_['P'][0], r['P'][0],2)
        self.assertAlmostEqual(r_['P'][1], r['P'][1],2)
        self.assertGreater(r_['P'][2], r['P'][2])

    def test_majority_gate(self):
        s = qca.QcaIsing()
        
        s.l = qca.MajorityGate(1,200,3.0,0,0,0)
        s.V0 = 1E6
        s.T = 1
        s.q = 0.5
        s.init()
        s.run()
        self.assertAlmostEqual(s.results['P'][-1], 0)

        Is = [(1,1,1),(-1,-1,-1),
              (1,-1,1),(-1,1,1),(-1,-1,1), # or gate with I3 fixed
              (-1,1,1),(1,1,-1),(-1,1,-1)] # or gate with I2 fixed
        Os = [1,-1,
              1,1,-1,
              1,1,-1]
        for I,O in zip(Is,Os):
            s.l = qca.MajorityGate(1,200,3.0,*I)
            s.V0 = 1E6
            s.T = 1
            s.q = 0.5
            s.init()
            s.run()
            # print(s.results['P'])
            # test the output polarization
            self.assertGreater(s.results['P'][-1] * O, 0.5)
