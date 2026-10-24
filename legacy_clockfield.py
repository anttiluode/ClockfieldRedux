"""Clock assay on the original nonlinear Clockfield force, without damping.

The main experiment isolates the conservative acceleration terms, but does not
assert a conserved Hamiltonian for this field-dependent non-divergence-form PDE.
An independent function reproduces the original damped scalar recurrence.
"""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Parameters:
    tau: float = 0.
    mu2: float = 1.4
    lam: float = .55
    length: float = 128.
    n: int = 256
    ds: float = .025
    epsilon: float = 1e-4

    def __post_init__(self):
        if self.tau<0 or self.n<64 or self.n%2 or min(self.mu2,self.lam,self.length,self.ds,self.epsilon)<=0:
            raise ValueError('Require tau >=0, even n>=64, and positive scales')
    @property
    def v(self): return np.sqrt(self.mu2/self.lam)
    @property
    def gamma0(self): return (1+self.tau*self.v**2)**-2
    @property
    def ceff0(self): return (1+self.tau*self.v**2)**-.5
    @property
    def dx(self): return self.length/self.n


def acceleration(phi,p):
    """Full nonlinear source force, expressed in s=Gamma0*t units."""
    beta=abs(phi)**2
    lap=(np.roll(phi,1,axis=-1)+np.roll(phi,-1,axis=-1)-2*phi)/p.dx**2
    return (1+p.tau*beta)**-4/p.gamma0**2*(lap/(1+p.tau*beta)+p.mu2*phi-p.lam*beta*phi)


def analytic_signal(y):
    n=len(y); h=np.zeros(n); h[0]=h[n//2]=1; h[1:n//2]=2
    return np.fft.ifft(np.fft.fft(y)*h)


def upward_crossings(t,y):
    i=np.flatnonzero((y[:-1]<0)&(y[1:]>=0))
    return t[i]-y[i]*(t[i+1]-t[i])/(y[i+1]-y[i])


def assay(p, traces=False):
    x=np.arange(p.n)*p.dx
    distance=(x-p.length/4+p.length/2)%p.length-p.length/2
    phase=p.epsilon*p.v*np.exp(-distance**2/(2*8**2))*np.cos(2*np.pi*12*x/p.length)
    phase-=phase.mean()
    q=np.full((2,p.n),p.v,complex)
    q[0]+=p.epsilon*p.v; q[1]+=1j*phase
    velocity=np.zeros_like(q)
    k=2*np.pi*np.fft.fftfreq(p.n,d=p.dx)
    # Right-moving lattice phase branch. The actual nonlinear dynamics thereafter
    # freely modifies both real and imaginary field components.
    omega=2*p.ceff0*np.sin(k*p.dx/2)/p.dx
    velocity[1]=1j*np.fft.ifft(-1j*omega*np.fft.fft(phase)).real
    force=acceleration(q,p)
    duration=max(45.,60/p.ceff0)
    steps=int(np.ceil(duration/p.ds)); stride=max(1,round(.1/p.ds))
    time=[]; radial=[]; center=[]; coherence=[]
    angle=np.exp(2j*np.pi*x/p.length)
    max_excursion=0.
    for step in range(steps+1):
        if step%stride==0:
            intensity=abs(analytic_signal(q[1].imag))**2
            z=np.sum(intensity*angle)/np.sum(intensity)
            time.append(step*p.ds/p.gamma0)
            radial.append(q[0].real.mean()-p.v)
            center.append(np.angle(z)); coherence.append(abs(z))
            max_excursion=max(max_excursion,float(np.max(abs(q-p.v))))
        if step==steps: break
        q+=p.ds*velocity+.5*p.ds**2*force
        next_force=acceleration(q,p)
        velocity+=.5*p.ds*(force+next_force)
        force=next_force
        if not np.all(np.isfinite(q)): raise ValueError('Nonfinite field; reduce ds')
    t=np.asarray(time); radial=np.asarray(radial)
    crossing=upward_crossings(t,radial)
    if len(crossing)<4: raise ValueError('Insufficient radial cycles')
    if min(coherence)<.5: raise ValueError('Packet centroid lost localization')
    pos=np.unwrap(center)*p.length/(2*np.pi)
    fit=np.polyfit(t,pos,1); residual=pos-np.polyval(fit,t)
    result=dict(radial_frequency=float(1/np.mean(np.diff(crossing))),pulse_speed=float(fit[0]),
                gamma0=p.gamma0,ceff0=p.ceff0,radial_crossings=len(crossing),
                centroid_coherence=float(min(coherence)),max_field_excursion=max_excursion,
                centroid_fit_rms=float(np.sqrt(np.mean(residual**2))),duration=float(t[-1]))
    if traces: return result,dict(time=t,radial=radial,position=pos-pos[0])
    return result


def original_damping_audit(p):
    """Exact source-style uniform radial recurrence, in original time units."""
    dt=.02; damping=.005
    q=p.v*(1+p.epsilon); previous=q
    values=[q-p.v]
    for _ in range(3000):
        beta=q*q
        force=(1+p.tau*beta)**-4*(p.mu2*q-p.lam*beta*q)
        nxt=2*q-previous+force*dt*dt-damping*(q-previous)
        previous,q=q,nxt; values.append(q-p.v)
    time=np.arange(len(values))*dt
    crossings=upward_crossings(time,np.array(values))
    omega=p.gamma0*np.sqrt(2*p.mu2); friction=damping/dt
    return dict(tau=p.tau,dt=dt,steps=3000,source_damping=damping,
                approximate_continuous_friction=friction,
                linearized_regime='underdamped' if omega>friction/2 else 'overdamped',
                upward_crossings=len(crossings),
                measured_frequency=float(1/np.mean(np.diff(crossings))) if len(crossings)>=4 else None,
                final_relative_displacement=float((q-p.v)/(p.v*p.epsilon)))
