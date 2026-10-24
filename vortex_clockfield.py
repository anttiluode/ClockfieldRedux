"""Undamped legacy complex field: static vortices, full Hessian and dynamics.

Static residual R = lap(phi) + F(|phi|²) phi, F(b)=(1+tau*b)(mu2-lam*b).
The original acceleration is a(phi)*R, a=(1+tau*|phi|²)^-5 > 0.
At equilibrium, -dR is symmetric in real Cartesian components. Therefore the
linear dynamical operator a*K is similar to sqrt(a)*K*sqrt(a). The signs of
K's eigenvalues decide linear stability. This does NOT assert a conserved
energy for the original nonlinear second-order dynamics.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import diags, eye, kron, bmat
from scipy.sparse.linalg import spsolve, eigsh
from scipy.linalg import eigh_tridiagonal


MU2, LAM = 1.4, .55


def radial_vortex(tau, radius=16., dr=.025):
    """Cell-centred radial finite volumes, winding one, f(radius)=v."""
    if tau < 0 or radius <= 0 or dr <= 0:
        raise ValueError('Positive radius/dr and nonnegative tau required')
    n = round(radius/dr)
    if n < 8:
        raise ValueError('Radial grid needs at least 8 cells')
    dr = radius/n
    r = (np.arange(n)+.5)*dr
    v = np.sqrt(MU2/LAM)
    lower = (r-dr/2)/(r*dr**2)
    upper = (r+dr/2)/(r*dr**2)
    diagonal = -lower-upper-1/r**2
    diagonal[-1] -= upper[-1]  # half-cell Dirichlet outer boundary
    bc = np.zeros(n); bc[-1] = 2*upper[-1]*v
    D = diags([lower[1:], diagonal, upper[:-1]], [-1,0,1], format='csc')
    f = v*np.tanh(r*np.sqrt(MU2/2))
    for _ in range(40):
        beta = f*f
        F = (1+tau*beta)*(MU2-LAM*beta)
        Fp = tau*MU2-LAM-2*tau*LAM*beta
        residual = D@f+bc+F*f
        if np.max(abs(residual)) < 1e-10:
            break
        delta = spsolve(D+diags(F+2*beta*Fp), -residual)
        scale = 1.
        while scale > 1e-7:
            trial = f+scale*delta
            b = trial*trial
            rr = D@trial+bc+(1+tau*b)*(MU2-LAM*b)*trial
            if np.max(abs(rr)) < np.max(abs(residual)):
                f = trial
                break
            scale *= .5
        else:
            raise RuntimeError('Radial Newton line search failed')
    else:
        raise RuntimeError('Radial equilibrium did not converge')
    return dict(tau=tau, r=r, f=f, dr=dr, radius=radius, v=v,
                diagonal=diagonal, upper=upper, lower=lower,
                residual=float(np.max(abs(residual))))


def radial_modes(profile, count=4):
    """Amplitude modes only; full complex stability is tested separately in 2D."""
    p = profile; tau=p['tau']; beta=p['f']**2
    F=(1+tau*beta)*(MU2-LAM*beta)
    Fp=tau*MU2-LAM-2*tau*LAM*beta
    a=(1+tau*beta)**-5
    d=a*(-p['diagonal']-F-2*beta*Fp)
    off=-np.sqrt(a[:-1]*a[1:]*p['upper'][:-1]*p['lower'][1:])
    w2, z = eigh_tridiagonal(d, off, select='i', select_range=(0,count-1))
    # z = sqrt(r/a) u, so physical field perturbation is sqrt(a/r) z.
    modes = np.sqrt(a/p['r'])[:,None]*z
    modes /= np.max(abs(modes), axis=0)
    threshold=(1+tau*p['v']**2)**-4*2*MU2
    r50=float(np.interp(.5*p['v'], p['f'], p['r']))
    weight=p['r']*modes[:,0]**2
    localization=float(np.sum(weight[p['r']<4*r50])/np.sum(weight))
    return dict(omega2=w2, modes=modes, threshold2=threshold,
                core_fraction=localization, r50=r50)


@dataclass
class VortexModel:
    tau: float = 0.
    half_width: float = 8.
    dx: float = .25
    mu2: float = MU2
    lam: float = LAM

    def __post_init__(self):
        if self.tau < 0 or min(self.half_width,self.dx,self.mu2,self.lam) <= 0:
            raise ValueError('Invalid model scales')
        # Radial helper uses the declared fixed potential.
        if self.mu2 != MU2 or self.lam != LAM:
            raise ValueError('This assay fixes mu2=1.4 and lambda=0.55')
        intervals=round(2*self.half_width/self.dx)
        if intervals < 8 or intervals % 2:
            raise ValueError('Require an even number of intervals >= 8')
        self.dx=2*self.half_width/intervals
        self.axis=np.linspace(-self.half_width,self.half_width,intervals+1)
        self.X,self.Y=np.meshgrid(self.axis,self.axis)
        self.r=np.hypot(self.X,self.Y)
        self.phase=np.exp(1j*np.arctan2(self.Y,self.X))
        self.shape=self.r.shape; self.n=intervals-1; self.interior_size=self.n**2
        self.v=np.sqrt(self.mu2/self.lam)
        self.gamma0=(1+self.tau*self.v**2)**-2
        self.xi=(1+self.tau*self.v**2)**-.5/np.sqrt(2*self.mu2)
        d=diags([-np.ones(self.n-1),2*np.ones(self.n),-np.ones(self.n-1)],[-1,0,1])
        self.minus_lap=(kron(eye(self.n),d)+kron(d,eye(self.n))).tocsc()/self.dx**2

    def initial(self):
        p=radial_vortex(self.tau, radius=2*self.half_width, dr=.025)
        f=np.interp(self.r, np.r_[0,p['r'],p['radius']], np.r_[0,p['f'],self.v])
        q=f*self.phase
        # Explicit fixed winding-one vacuum boundary, not an absorber.
        q[[0,-1],:]=self.v*self.phase[[0,-1],:]
        q[:,[0,-1]]=self.v*self.phase[:,[0,-1]]
        return q

    def pack(self, q):
        z=q[1:-1,1:-1].ravel()
        return np.r_[z.real,z.imag]

    def unpack(self, z):
        q=np.zeros(self.shape,complex); k=self.interior_size
        q[1:-1,1:-1]=(z[:k]+1j*z[k:]).reshape(self.n,self.n)
        return q

    def static_force(self, q):
        u=q[1:-1,1:-1]; beta=abs(u)**2
        lap=(q[2:,1:-1]+q[:-2,1:-1]+q[1:-1,2:]+q[1:-1,:-2]-4*u)/self.dx**2
        out=np.zeros_like(q)
        out[1:-1,1:-1]=lap+(1+self.tau*beta)*(self.mu2-self.lam*beta)*u
        return out

    def acceleration(self, q):
        return (1+self.tau*abs(q)**2)**-5*self.static_force(q)

    def hessian(self, q):
        u=q[1:-1,1:-1].ravel(); beta=abs(u)**2
        F=(1+self.tau*beta)*(self.mu2-self.lam*beta)
        Fp=self.tau*self.mu2-self.lam-2*self.tau*self.lam*beta
        xx=self.minus_lap-diags(F+2*Fp*u.real**2)
        yy=self.minus_lap-diags(F+2*Fp*u.imag**2)
        xy=diags(-2*Fp*u.real*u.imag)
        return bmat([[xx,xy],[xy,yy]],format='csc')

    def equilibrium(self):
        q=self.initial()
        for _ in range(20):
            rr=self.static_force(q); error=float(np.max(abs(rr))/self.v)
            if error < 1e-10:
                return q,error
            delta=self.unpack(spsolve(self.hessian(q), self.pack(rr)))
            scale=1.
            while scale > 1e-7:
                trial=q+scale*delta
                if np.max(abs(self.static_force(trial))) < np.max(abs(rr)):
                    q=trial; break
                scale*=.5
            else:
                raise RuntimeError('Square Newton line search failed')
        raise RuntimeError('Square equilibrium did not converge')

    def stability(self, q):
        K=self.hessian(q)
        w,z=eigsh(K,k=3,which='SA',tol=1e-9,v0=np.random.default_rng(0).normal(size=K.shape[0]))
        residuals=[np.linalg.norm(K@z[:,i]-w[i]*z[:,i]) for i in range(len(w))]
        return dict(hessian_eigenvalues=w.tolist(),eigen_residuals=residuals,
                    stable=bool(w[0]>=-1e-7 and max(residuals)<1e-7))

    def winding(self, q):
        # Counterclockwise contour inside the fixed boundary.
        i=max(2,self.shape[0]//4); j=self.shape[0]-1-i
        contour=np.r_[q[i,i:j],q[i:j,j],q[j,j:i:-1],q[j:i:-1,i]]
        return float(np.sum(np.angle(np.roll(contour,-1)*np.conj(contour)))/(2*np.pi))

    def radii(self, q):
        mid=self.shape[0]//2
        amplitude=abs(q[mid,mid:])/self.v
        positive=self.axis[mid:]
        return {f'r{int(100*f)}':float(np.interp(f,amplitude,positive)) for f in [.5,.9]}

    def evolve(self, equilibrium, disturbance, duration, dt=.025, projection=None):
        """Velocity Verlet, zero initial velocity; boundaries remain fixed."""
        q=equilibrium+disturbance
        velocity=np.zeros_like(q); force=self.acceleration(q)
        steps=int(np.ceil(duration/dt)); dt=duration/steps
        stride=max(1,round(.1/dt))
        t=[]; signal=[]; rms=[]; wind=[]
        weight=(1+self.tau*abs(equilibrium)**2)**5
        norm=np.sum(weight*abs(projection)**2) if projection is not None else 1.
        initial_rms=np.sqrt(np.mean(abs(disturbance[1:-1,1:-1])**2))
        for step in range(steps+1):
            if step%stride==0 or step==steps:
                delta=q-equilibrium
                t.append(step*dt)
                rms.append(float(np.sqrt(np.mean(abs(delta[1:-1,1:-1])**2))))
                wind.append(self.winding(q))
                signal.append(float(np.sum(weight*np.real(np.conj(projection)*delta))/norm) if projection is not None else 0.)
            if step==steps: break
            q+=dt*velocity+.5*dt*dt*force
            new=self.acceleration(q); velocity+=.5*dt*(force+new); force=new
            if not np.all(np.isfinite(q)):
                raise RuntimeError('Nonfinite field')
        return dict(time=np.asarray(t),signal=np.asarray(signal),rms=np.asarray(rms),
                    winding=np.asarray(wind),initial_rms=float(initial_rms),dt=dt)


def measured_omega(trace):
    t=trace['time']; y=trace['signal']
    indices=np.flatnonzero((y[:-1]<0)&(y[1:]>=0))
    crossings=t[indices]-y[indices]*(t[indices+1]-t[indices])/(y[indices+1]-y[indices])
    return dict(omega=float(2*np.pi/np.mean(np.diff(crossings))) if len(crossings)>=4 else None,
                crossings=len(crossings),period_scatter=float(np.std(np.diff(crossings))/np.mean(np.diff(crossings))) if len(crossings)>=4 else None)
