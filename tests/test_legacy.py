import unittest
import numpy as np
from legacy_clockfield import Parameters, acceleration, assay

class LegacyTests(unittest.TestCase):
    def test_equilibrium(self):
        p=Parameters(tau=.4)
        np.testing.assert_allclose(acceleration(np.full(64,p.v,complex),p),0,atol=1e-12)
    def test_two_dimensional_slice(self):
        p=Parameters(tau=.2)
        rng=np.random.default_rng(41)
        q=p.v+1e-3*(rng.normal(size=64)+1j*rng.normal(size=64))
        f=np.tile(q,(8,1)); beta=abs(f)**2
        lap=(np.roll(f,1,0)+np.roll(f,-1,0)+np.roll(f,1,1)+np.roll(f,-1,1)-4*f)/p.dx**2
        expected=((1+p.tau*beta)**-4)*(lap/(1+p.tau*beta)+p.mu2*f-p.lam*beta*f)/p.gamma0**2
        np.testing.assert_allclose(acceleration(q,p),expected[0],atol=1e-13)
    def test_linearized_radial_force(self):
        p=Parameters(tau=.4); e=1e-7
        derivative=(acceleration(np.full(64,p.v+e,complex),p)-acceleration(np.full(64,p.v-e,complex),p))/(2*e)
        np.testing.assert_allclose(derivative,-2*p.mu2,rtol=1e-7)
    def test_measured_clock_ratio(self):
        a=assay(Parameters(tau=0))
        b=assay(Parameters(tau=.4))
        ratio=(b['pulse_speed']/a['pulse_speed'])/(b['radial_frequency']/a['radial_frequency'])
        self.assertAlmostEqual(ratio,1/np.sqrt(1+.4*1.4/.55),delta=.003)

if __name__=='__main__': unittest.main()
