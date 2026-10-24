import unittest
import numpy as np
from clockfield import Model, evolve, energy, measure

class PhysicsTests(unittest.TestCase):
    def test_conservation_and_reversal(self):
        rng = np.random.default_rng(41)
        q = rng.normal(size=(2,64)) + 1j*rng.normal(size=(2,64))
        p = rng.normal(size=(2,64)) + 1j*rng.normal(size=(2,64))
        m = Model(n=64)
        qt,pt = evolve(m,q,p,np.array([0.,3.7]))
        self.assertAlmostEqual(energy(m,qt[0],pt[0]),energy(m,qt[1],pt[1]),places=9)
        qr,pr = evolve(m,qt[1],pt[1],np.array([-3.7]))
        np.testing.assert_allclose(qr[0],q,atol=1e-12)
        np.testing.assert_allclose(pr[0],p,atol=1e-12)
    def test_clocks_and_positive_control(self):
        a = measure(Model())
        b = measure(Model(h=.65))
        for key in ('oscillation','exchange','pulse'):
            self.assertAlmostEqual(b[key]/a[key],.65,delta=.002)
        self.assertAlmostEqual(a['oscillation'],1/(2*np.pi),delta=1e-5)
        self.assertAlmostEqual(a['exchange'],(np.sqrt(2.6)-1)/(2*np.pi),delta=1e-5)
    def test_selective_control_is_detected(self):
        a,b = measure(Model()),measure(Model(c=.65))
        self.assertAlmostEqual(b['oscillation']/a['oscillation'],1,delta=.001)
        self.assertLess(b['pulse']/a['pulse'],.8)
    def test_pulse_matches_group_velocity(self):
        m=Model(); k=np.pi/2
        w=np.sqrt(m.mass**2+(2*m.c*np.sin(k*m.dx/2)/m.dx)**2)
        expected=m.h*m.c**2*np.sin(k*m.dx)/(m.dx*w)
        self.assertAlmostEqual(measure(m)["pulse_speed"],expected,delta=.005)

    def test_invalid_model(self):
        with self.assertRaises(ValueError): Model(c=0)

if __name__ == '__main__': unittest.main()
