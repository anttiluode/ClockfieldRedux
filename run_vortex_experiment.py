"""Run the declared vortex assay; CPU only. Outputs retain failed gates."""
import argparse
import json
from pathlib import Path
import platform
import time
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from vortex_clockfield import VortexModel, radial_vortex, radial_modes, measured_omega


def relative(a,b):
    return float(abs(a/b-1)) if a is not None and b is not None else None


def percentage(value, digits=3):
    return 'unresolved' if value is None else f'{100*value:.{digits}f}%'


def summary_html(lines):
    """Render this runner's small, trusted Markdown subset without dependencies."""
    import html
    import re
    def inline(text):
        text=html.escape(text)
        text=re.sub(r'\*\*(.*?)\*\*',r'<strong>\1</strong>',text)
        return re.sub(r'\[([^]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',text)
    parts=[]; table=False
    for line in lines:
        if line.startswith('|'):
            cells=[c.strip() for c in line.strip('|').split('|')]
            if all(set(c)<=set('-: ') for c in cells): continue
            tag='td' if table else 'th'
            if not table: parts.append('<div class="scroll"><table>'); table=True
            parts.append('<tr>'+''.join(f'<{tag}>{inline(c)}</{tag}>' for c in cells)+'</tr>')
        else:
            if table: parts.append('</table></div>'); table=False
            if line.startswith('# '): continue
            if line.startswith('## '): parts.append('<h2>'+inline(line[3:])+'</h2>')
            elif line: parts.append('<p>'+inline(line)+'</p>')
    if table: parts.append('</table></div>')
    return ''.join(parts)


def mode_on_square(m,p,s):
    u=np.interp(m.r, np.r_[0,p['r'],p['radius']], np.r_[0,s['modes'][:,0],0])
    mode=u*m.phase
    mode[[0,-1],:]=0; mode[:,[0,-1]]=0
    return mode


def run_case(tau,dx,half_width,radial,extra_checks=False):
    m=VortexModel(tau=tau,dx=dx,half_width=half_width)
    q,error=m.equilibrium()
    stability=m.stability(q)
    p,s=radial; omega=float(np.sqrt(s['omega2'][0])); vacuum=np.sqrt(s['threshold2'])
    mode=mode_on_square(m,p,s)
    duration=8*2*np.pi/vacuum
    trace=m.evolve(q,1e-4*m.v*mode,duration,projection=mode)
    clock=measured_omega(trace)
    rng=np.random.default_rng(41)
    noise=(rng.normal(size=m.shape)+1j*rng.normal(size=m.shape))*np.exp(-m.r**2/8)
    noise[[0,-1],:]=0; noise[:,[0,-1]]=0
    noise*=1e-3*m.v/np.max(abs(noise))
    survival=m.evolve(q,noise,duration)
    survival_result=dict(seed=41,max_rms_growth=float(max(survival['rms'])/survival['initial_rms']),
                         max_winding_error=float(max(abs(survival['winding']-1))),duration=duration)
    radii=m.radii(q)
    result=dict(tau=tau,dx=m.dx,half_width=half_width,grid=m.shape[0],static_residual=error,
                **radii,xi=m.xi,r50_over_xi=radii['r50']/m.xi,r90_over_xi=radii['r90']/m.xi,
                stability=stability,clock=clock,vacuum_omega=float(vacuum),
                measured_clock_ratio=clock['omega']/vacuum if clock['omega'] is not None else None,
                radial_prediction_omega=omega,prediction_error=relative(clock['omega'],omega),
                survival=survival_result,dt=trace['dt'],amplitude=1e-4,
                trace={k:trace[k].tolist() for k in ['time','signal','rms','winding']},
                noise_trace={k:survival[k].tolist() for k in ['time','rms','winding']})
    if extra_checks:
        fine=m.evolve(q,1e-4*m.v*mode,duration,dt=.0125,projection=mode)
        small=m.evolve(q,5e-5*m.v*mode,duration,projection=mode)
        result['refinements']=dict(half_dt=measured_omega(fine),half_amplitude=measured_omega(small))
        result['refinements']['half_dt_error']=relative(result['refinements']['half_dt']['omega'],clock['omega'])
        result['refinements']['half_amplitude_error']=relative(result['refinements']['half_amplitude']['omega'],clock['omega'])
    result['checks']=dict(equilibrium=error<1e-8,linear_stability=stability['stable'],
                          survival=survival_result['max_rms_growth']<5 and survival_result['max_winding_error']<1e-8,
                          resolved_oscillation=clock['crossings']>=4,
                          prediction=clock['omega'] is not None and result['prediction_error']<.02)
    if extra_checks:
        for key in ['half_dt','half_amplitude']:
            err=result['refinements'][key+'_error']
            result['checks'][key]=err is not None and err<.01
    return result,(m,q)


def radial_summary(tau):
    profiles=[radial_vortex(tau,radius=R,dr=h) for R,h in [(16,.025),(32,.025),(16,.0125)]]
    spectra=[radial_modes(p) for p in profiles]
    s=spectra[0]
    box=relative(np.sqrt(spectra[1]['omega2'][0]),np.sqrt(s['omega2'][0]))
    grid=relative(np.sqrt(spectra[2]['omega2'][0]),np.sqrt(s['omega2'][0]))
    candidate=bool(0<s['omega2'][0]<.99*s['threshold2'] and box<.02 and grid<.02)
    rows=[]
    for p,sp in zip(profiles,spectra):
        rows.append(dict(radius=p['radius'],dr=p['dr'],omega2=sp['omega2'].tolist(),
                         omega2_over_threshold2=(sp['omega2']/sp['threshold2']).tolist(),
                         core_fraction=sp['core_fraction'],static_residual=p['residual']))
    return dict(tau=tau,candidate_bound_clock=candidate,box_error=box,grid_error=grid,
                radial_clock_ratio=float(np.sqrt(s['omega2'][0]/s['threshold2'])),runs=rows,
                profile=dict(r=profiles[0]['r'].tolist(),amplitude=(profiles[0]['f']/profiles[0]['v']).tolist(),
                             mode=s['modes'][:,0].tolist())),(profiles[0],s)


def reports(out,data,visual):
    base=[r for r in data['cases'] if r['dx']==.25 and r['half_width']==8]
    lines=['# Vortex: object, ruler, clock','',
           'Full complex legacy field, undamped. Fixed winding-one outer boundary. '
           'Finite-box stability and finite-duration survival are not infinite-domain stability proofs.','',
           '| tau | r50 / xi | r90 / xi | Radial mode / vacuum | Bound clock candidate? | Full complex linear stability |',
           '|---|---:|---:|---:|---|---|']
    for r,s in zip(base,data['radial']):
        lines.append(f"| {r['tau']:g} | {r['r50_over_xi']:.4f} | {r['r90_over_xi']:.4f} | {s['radial_clock_ratio']:.4f} | {'yes' if s['candidate_bound_clock'] else 'not resolved'} | {'pass' if r['stability']['stable'] else 'FAIL'} |")
    lines+=['','The spectral ratio in the table is a radial-sector calculation. At tau=0.4 the lowest '
            'reported radial mode is extended and above the vacuum threshold; it is not an internal clock. '
            'Full 2D time measurements and all convergence checks follow.','',
            '| tau | dx | half-width | Measured omega / vacuum | Prediction error | Noise RMS growth |',
            '|---|---:|---:|---:|---:|---:|']
    for r in data['cases']:
        ratio='unresolved' if r['measured_clock_ratio'] is None else f"{r['measured_clock_ratio']:.5f}"
        err='unresolved' if r['prediction_error'] is None else f"{100*r['prediction_error']:.3f}%"
        lines.append(f"| {r['tau']:g} | {r['dx']:g} | {r['half_width']:g} | {ratio} | {err} | {r['survival']['max_rms_growth']:.3f} |")
    lines+=['','## Checks','']
    for c in data['comparisons']:
        lines.append(f"- tau={c['tau']:g}, {c['kind']}: r50 {percentage(c['r50'])}, r90 {percentage(c['r90'])}, "
                     f"oscillation {percentage(c['omega'])} — {'pass' if c['pass'] else 'FAIL'} (2% limit).")
    for r in base:
        lines.append(f"- tau={r['tau']:g}: half-step {percentage(r['refinements']['half_dt_error'],4)}, "
                     f"half-amplitude {percentage(r['refinements']['half_amplitude_error'],4)} (1% limit).")
    lines+=['',f"All declared numerical checks: **{'PASS' if data['all_numerical_checks'] else 'NOT ALL PASS'}**.",
            'A successful stability check is not a successful universal-clock test. Changing tau changes '
            'dimensionless internal structure. A characteristic wave length is not automatically a material ruler.','',
            'An above-threshold radial mode is a box oscillation, not a bound-clock candidate; its '
            'frequency agreement gate is still recorded, never silently waived.','',
            'See [VORTEX_MODEL.md](../VORTEX_MODEL.md), [protocol](../VORTEX_PROTOCOL.md), '
            'and [raw results](results.json).']
    summary='\n'.join(lines)+'\n'; (out/'SUMMARY.md').write_text(summary)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(15,8.5),layout='constrained')
    m,q=visual
    im=axes[0,0].imshow(abs(q)/m.v,origin='lower',extent=[-8,8,-8,8],cmap='magma',vmin=0,vmax=1)
    axes[0,0].set(title='A complex vortex (tau = 0.1)',xlabel='x / lattice unit',ylabel='y / lattice unit')
    fig.colorbar(im,ax=axes[0,0],label='Field amplitude / vacuum')
    for s in data['radial']:
        p=s['profile']; color={0:'#147d92',.1:'#d47b15',.4:'#8c58a6'}[s['tau']]
        axes[0,1].plot(p['r'],p['amplitude'],label=f"tau={s['tau']:g}",color=color)
        axes[0,2].plot(p['r'],np.abs(p['mode']),label=f"tau={s['tau']:g}",color=color)
    axes[0,1].set(xlim=(0,5),ylim=(0,1.05),title='An object supplies more than one length',xlabel='Radius / lattice unit',ylabel='Field amplitude / vacuum')
    axes[0,1].axhline(.5,color='gray',lw=.7,ls=':'); axes[0,1].axhline(.9,color='gray',lw=.7,ls=':'); axes[0,1].legend()
    axes[0,2].set(title='Breathing mode loses localization',xlabel='Radius / lattice unit',ylabel='Mode amplitude (peak = 1)'); axes[0,2].legend()
    for r in base:
        color={0:'#147d92',.1:'#d47b15',.4:'#8c58a6'}[r['tau']]
        t=np.asarray(r['trace']['time'])*r['vacuum_omega']/(2*np.pi)
        axes[1,0].plot(t,np.asarray(r['trace']['signal'])/(1e-4*np.sqrt(1.4/.55)),color=color,label=f"tau={r['tau']:g}")
        nt=np.asarray(r['noise_trace']['time'])*r['vacuum_omega']/(2*np.pi)
        nr=np.asarray(r['noise_trace']['rms'])
        axes[1,2].plot(nt,nr/nr[0],color=color,label=f"tau={r['tau']:g}")
    axes[1,0].set(title='Full 2D field: measured oscillations',xlabel='Time / vacuum radial period',ylabel='Mode projection / initial amplitude'); axes[1,0].legend()
    taus=[r['tau'] for r in base]
    axes[1,1].plot(taus,[r['r50_over_xi'] for r in base],'o-',label='Core radius / xi')
    axes[1,1].plot(taus,[r['r90_over_xi'] for r in base],'s-',label='90% radius / xi')
    axes[1,1].set(title='Internal ratios change',xlabel='Global coupling tau',ylabel='Dimensionless length'); axes[1,1].legend()
    axes[1,2].set(title='Independent complex-noise survival test',xlabel='Time / vacuum radial period',ylabel='Perturbation RMS / initial RMS'); axes[1,2].legend()
    fig.suptitle('ClockfieldRedux — a vortex can survive while its clock changes',fontsize=17)
    fig.savefig(out/'vortex.svg'); fig.savefig(out/'preview.png',dpi=110); plt.close(fig)
    svg=(out/'vortex.svg').read_text(); svg=svg[svg.index('<svg'):]
    body=summary_html(lines)
    (out/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>ClockfieldRedux vortex assay</title><style>body{font:17px/1.55 system-ui;max-width:1250px;margin:32px auto;padding:0 20px;background:#f6f5f1;color:#192b33}svg{width:100%;height:auto}p{overflow-wrap:anywhere}h1{line-height:1.15}.scroll{overflow-x:auto}table{border-collapse:collapse;width:100%;font-size:15px;background:white}th,td{padding:10px 12px;border-bottom:1px solid #d6dfdf;text-align:left}th{background:#e0ebed}a{color:#08637e}</style><h1>A vortex: object, ruler, clock</h1><p>Measured full-complex-field experiment. Finite boxes, fixed outer winding, no gravity claim.</p>'+svg+body+'</html>')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',default='vortex_results')
    args=parser.parse_args(); out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    start=time.perf_counter()
    data=dict(model='undamped legacy complex field',mu2=1.4,lam=.55,
              environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
              radial=[],cases=[],comparisons=[])
    visual=None
    for tau in [0.,.1,.4]:
        rs,radial=radial_summary(tau); data['radial'].append(rs)
        runs=[]
        for dx,R in [(.25,8.),(.125,8.),(.25,16.)]:
            print(f'tau={tau:g}, dx={dx:g}, half-width={R:g}',flush=True)
            result,v=run_case(tau,dx,R,radial,extra_checks=(dx==.25 and R==8))
            data['cases'].append(result); runs.append(result)
            if tau==.1 and dx==.25 and R==8: visual=v
            print('  checks:',result['checks'],flush=True)
        for other,kind in zip(runs[1:],['half grid spacing','double box']):
            c=dict(tau=tau,kind=kind,**{k:relative(other[k],runs[0][k]) for k in ['r50','r90']},
                   omega=relative(other['clock']['omega'],runs[0]['clock']['omega']))
            c['pass']=all(c[k] is not None and c[k]<.02 for k in ['r50','r90','omega'])
            data['comparisons'].append(c)
    data['all_numerical_checks']=bool(all(all(r['checks'].values()) for r in data['cases']) and all(c['pass'] for c in data['comparisons']) and all(r['box_error']<.02 and r['grid_error']<.02 for r in data['radial']))
    data['elapsed_seconds']=time.perf_counter()-start
    # Full-resolution traces, compact encoding to stay below GitHub's 1 MB file API limit.
    (out/'results.json').write_text(json.dumps(data,separators=(',',':'),allow_nan=False)+'\n')
    print(reports(out,data,visual))
    # Nonzero exit records a gate failure; files remain available for inspection.
    return 0 if data['all_numerical_checks'] else 1


if __name__=='__main__':
    raise SystemExit(main())
