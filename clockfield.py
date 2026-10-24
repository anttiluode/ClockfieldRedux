"""Conservative periodic two-component Klein–Gordon lattice, dimensionless units."""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Model:
    n: int = 256
    length: float = 128.
    mass: float = 1.
    c: float = 1.
    coupling: float = .8
    h: float = 1.

    def __post_init__(self):
        if self.n < 64 or self.n % 2 or min(self.length,self.mass,self.c,self.coupling,self.h) <= 0:
            raise ValueError('Require even n >= 64 and positive physical parameters')

    @property
    def dx(self): return self.length/self.n

    def frequencies(self):
        k = 2*np.pi*np.fft.fftfreq(self.n,d=self.dx)
        base = self.mass**2 + (2*self.c*np.sin(k*self.dx/2)/self.dx)**2
        return np.sqrt(np.array([base,base+2*self.coupling]))


def modes(a):
    return np.stack((a[0]+a[1],a[0]-a[1]))/np.sqrt(2)


def evolve(m,q,p,times):
    """Exact lattice evolution sampled at times (including negative times).

    Complex coordinates represent two independent real copies of the same model.
    h multiplies the ENTIRE Hamiltonian: qdot=h*p, pdot=-h*K*q.
    """
    w = m.frequencies()
    qk,pk = np.fft.fft(modes(q)),np.fft.fft(modes(p))
    phase = np.asarray(times)[:,None,None]*m.h*w
    co,si = np.cos(phase),np.sin(phase)
    qm = np.fft.ifft(qk*co+pk*si/w,axis=-1)
    pm = np.fft.ifft(pk*co-qk*w*si,axis=-1)
    return np.stack((qm[:,0]+qm[:,1],qm[:,0]-qm[:,1]),axis=1)/np.sqrt(2), np.stack((pm[:,0]+pm[:,1],pm[:,0]-pm[:,1]),axis=1)/np.sqrt(2)


def energy(m,q,p):
    grad = (np.roll(q,-1,axis=-1)-q)/m.dx
    return float(m.h*m.dx/2*(np.sum(abs(p)**2+m.mass**2*abs(q)**2+m.c**2*abs(grad)**2)+m.coupling*np.sum(abs(q[0]-q[1])**2)))


def frequency(t,y):
    """Frequency from interpolated upward zero crossings, not a spectral formula."""
    i = np.flatnonzero((y[:-1]<0)&(y[1:]>=0))
    if len(i)<4: raise ValueError('Too few cycles for frequency measurement')
    crossing = t[i]-y[i]*(t[i+1]-t[i])/(y[i+1]-y[i])
    return float(1/np.mean(np.diff(crossing)))


def measure(m,samples=2048,return_traces=False):
    # Fixed observation windows in reference units, adjusted only for the known
    # imposed h control. No measured frequency is supplied to the estimator.
    t = np.linspace(0,160/m.h,samples)
    q = np.ones((2,m.n),complex)
    p = np.zeros_like(q)
    qo,po = evolve(m,q,p,t)
    osc = qo[:,0,0].real
    q = np.zeros_like(q); q[0]=1
    # Analytic-signal initialization excites both normal branches positively.
    p = modes(np.fft.ifft(-1j*m.frequencies()*np.fft.fft(modes(q))))
    qe,pe = evolve(m,q,p,t)
    exchange = abs(qe[:,0,0])**2-abs(qe[:,1,0])**2
    x = np.arange(m.n)*m.dx
    distance = (x-m.length/4+m.length/2)%m.length-m.length/2
    packet = np.exp(-distance**2/(2*8.**2))*np.exp(1j*(np.pi/2)*x)
    q = np.stack((packet,packet))/np.sqrt(2)
    p = modes(np.fft.ifft(-1j*m.frequencies()*np.fft.fft(modes(q))))
    tp = np.linspace(0,70/m.h,samples//2)
    qp,pp = evolve(m,q,p,tp)
    intensity = np.sum(abs(qp)**2,axis=1)
    moment = (intensity @ np.exp(2j*np.pi*x/m.length))/intensity.sum(axis=1)
    if np.min(abs(moment))<.5: raise ValueError('Packet too delocalized for circular centroid')
    position = np.unwrap(np.angle(moment))*m.length/(2*np.pi)
    speed = float(np.polyfit(tp,position,1)[0])
    result = dict(oscillation=frequency(t,osc),exchange=frequency(t,exchange),pulse=speed/m.length,
                  pulse_speed=speed,centroid_coherence=float(min(abs(moment))),
                  relative_energy_drift=abs(energy(m,qp[-1],pp[-1])/energy(m,q,p)-1))
    if return_traces: return result,dict(t=t,osc=osc,exchange=exchange,tp=tp,position=position,x=x,intensity=intensity)
    return result
