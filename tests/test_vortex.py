"""Checks that catch force, component-coupling and eigenvalue-sign mistakes."""
import unittest
import numpy as np
from scipy.sparse.linalg import eigsh
from vortex_clockfield import VortexModel, radial_vortex, radial_modes, measured_omega


class VortexTests(unittest.TestCase):
    def test_full_time_evolution_matches_vacuum_lattice_mode(self):
        m = VortexModel(tau=0, half_width=3, dx=.5)
        q = np.full(m.shape, m.v, complex)
        mode = np.cos(np.pi*m.X/6)*np.cos(np.pi*m.Y/6)+0j
        mode[[0,-1],:]=0; mode[:,[0,-1]]=0
        expected = np.sqrt(2*m.mu2+8*np.sin(np.pi/(2*(m.n+1)))**2/m.dx**2)
        trace = m.evolve(q, 1e-5*mode, 6*2*np.pi/expected, dt=.01, projection=mode)
        self.assertLess(abs(measured_omega(trace)['omega']/expected-1), 1e-4)

    def test_global_phase_rotation_preserves_force(self):
        m = VortexModel(tau=.4, half_width=3, dx=.5)
        q = m.initial(); phase=np.exp(.73j)
        np.testing.assert_allclose(m.acceleration(q*phase), m.acceleration(q)*phase, atol=1e-13)

    def test_force_matches_original_complex_equation(self):
        m = VortexModel(tau=.4, half_width=3, dx=.5)
        rng = np.random.default_rng(17)
        q = m.initial() + .03*(rng.normal(size=m.shape)+1j*rng.normal(size=m.shape))
        beta = abs(q[1:-1,1:-1])**2
        lap = (q[2:,1:-1]+q[:-2,1:-1]+q[1:-1,2:]+q[1:-1,:-2]-4*q[1:-1,1:-1])/m.dx**2
        expected = (1+m.tau*beta)**-4*(lap/(1+m.tau*beta)+m.mu2*q[1:-1,1:-1]-m.lam*beta*q[1:-1,1:-1])
        np.testing.assert_allclose(m.acceleration(q)[1:-1,1:-1], expected, atol=1e-13)

    def test_hessian_matches_both_component_derivatives(self):
        m = VortexModel(tau=.1, half_width=3, dx=.5)
        q = m.equilibrium()[0]
        rng = np.random.default_rng(9)
        z = rng.normal(size=2*m.interior_size)
        dq = m.unpack(z)
        eps = 1e-6
        numerical = m.pack((m.static_force(q+eps*dq)-m.static_force(q-eps*dq))/(2*eps))
        np.testing.assert_allclose(numerical, -m.hessian(q)@z, rtol=1e-7, atol=1e-7)

    def test_complex_kink_negative_mode_is_not_abs_hidden(self):
        # Exact tau=0 transverse Pöschl–Teller eigenvalue: -mu2/2.
        m = VortexModel(tau=0, half_width=8, dx=.25)
        q = m.v*np.tanh(np.sqrt(m.mu2/2)*m.X)+0j
        w = eigsh(m.hessian(q), k=1, which='SA', return_eigenvectors=False)[0]
        # Finite y box adds approximately pi^2/(16^2).
        self.assertLess(w, -.6)

    def test_equilibrium_and_radial_discretization_converge(self):
        m = VortexModel(tau=.1, half_width=4, dx=.25)
        q, residual = m.equilibrium()
        self.assertLess(residual, 1e-8)
        self.assertAlmostEqual(m.winding(q), 1., places=10)
        a = radial_modes(radial_vortex(.1, radius=12, dr=.05))
        b = radial_modes(radial_vortex(.1, radius=12, dr=.025))
        self.assertGreater(a['omega2'][0], 0)
        self.assertLess(abs(a['omega2'][0]/b['omega2'][0]-1), .01)


if __name__ == '__main__':
    unittest.main()
