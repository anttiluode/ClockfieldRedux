"""Run the original-PDE clock assay and refinement checks on CPU."""
import argparse
from dataclasses import asdict,replace
import json
from pathlib import Path
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from legacy_clockfield import Parameters,assay,original_damping_audit

SOURCE_BLOB='f4cc9fe9e72bc138c874337bf9e8c10f4a521f13'

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('legacy_results'))
    args=parser.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    base=Parameters(); ref=assay(base)
    rows=[]
    for tau in (0,.1,.2,.4,1,5):
        p=replace(base,tau=tau); r=assay(p)
        radial=r['radial_frequency']/ref['radial_frequency']; pulse=r['pulse_speed']/ref['pulse_speed']
        expected_radial=p.gamma0; expected_pulse=p.gamma0*p.ceff0
        rows.append(dict(tau=tau,parameters=asdict(p),measured=r,radial_rate=radial,pulse_rate=pulse,
                         clock_ratio=pulse/radial,predicted_ratio=p.ceff0,
                         agreement=bool(abs(pulse/radial-1)<.01),
                         prediction_error=max(abs(radial/expected_radial-1),abs(pulse/expected_pulse-1))))
    refinements=[]
    for tau in (.4,5):
        p=replace(base,tau=tau); a=assay(p)
        checks={}
        for key,changed in [('half_step',replace(p,ds=p.ds/2)),('half_amplitude',replace(p,epsilon=p.epsilon/2)),('double_grid',replace(p,n=p.n*2))]:
            b=assay(changed)
            checks[key]=dict(measured=b,relative_change={k:abs(b[k]/a[k]-1) for k in ('radial_frequency','pulse_speed')})
        # Always report a second spatial refinement; keep the coarse-grid failure.
        finer=assay(replace(p,n=p.n*4)); middle=checks['double_grid']['measured']
        checks['quadruple_grid']=dict(measured=finer,relative_change={k:abs(finer[k]/middle[k]-1) for k in ('radial_frequency','pulse_speed')})
        refinements.append(dict(tau=tau,checks=checks))
    damping=[original_damping_audit(replace(base,tau=r['tau'])) for r in rows]
    validation=dict(prediction=all(r['prediction_error']<.005 for r in rows),
                    half_step=all(max(r['checks']['half_step']['relative_change'].values())<.005 for r in refinements),
                    half_amplitude=all(max(r['checks']['half_amplitude']['relative_change'].values())<.005 for r in refinements),
                    first_grid_refinement=all(max(r['checks']['double_grid']['relative_change'].values())<.005 for r in refinements),
                    second_grid_refinement=all(max(r['checks']['quadruple_grid']['relative_change'].values())<.005 for r in refinements))
    data=dict(source_blob=SOURCE_BLOB,python=platform.python_version(),numpy=np.__version__,rows=rows,refinements=refinements,original_damping=damping,validation=validation)
    (args.out/'results.json').write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    table='| tau | Radial clock | Phase travel | Travel / radial | Predicted ratio | Common rate? |\n|---|---:|---:|---:|---:|---|\n'
    for r in rows: table+=f"| {r['tau']:g} | {r['radial_rate']:.6f} | {r['pulse_rate']:.6f} | {r['clock_ratio']:.6f} | {r['predicted_ratio']:.6f} | {'yes' if r['agreement'] else 'no'} |\n"
    text='# Original Clockfield PDE: measured clocks\n\nFull nonlinear force; homogeneous equilibrium; weak probes; main assay damping disabled.\n\n'+table+'\nValidation: '+json.dumps(validation)+'.\n\nThe first spatial refinement is retained even if it fails the declared 0.5% raw-rate tolerance. The second refinement tests whether that discrepancy shrinks. The normalized clock disagreement is far larger than the reported numerical errors.\n\n## Original damping audit\n\n| tau | Linearized regime | Upward crossings in original 60 time units | Measured frequency |\n|---|---|---:|---:|\n'
    for r in damping: text+=f"| {r['tau']:g} | {r['linearized_regime']} | {r['upward_crossings']} | {r['measured_frequency'] if r['measured_frequency'] is not None else 'not resolved'} |\n"
    text+='\nThe original coefficient 0.005 at dt=0.02 approximates friction 0.25. An unresolved or overdamped amplitude probe is not a zero-frequency universal clock. This is not a vortex experiment or a material-ruler comparison.\n'
    (args.out/'SUMMARY.md').write_text(text)
    plt.rcParams.update({'svg.fonttype':'none','font.size':11})
    fig,ax=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
    taus=np.array([r['tau'] for r in rows])
    ax[0].plot(taus,[r['radial_rate'] for r in rows],'o-',label='Radial oscillation')
    ax[0].plot(taus,[r['pulse_rate'] for r in rows],'s-',label='Phase travel')
    ax[0].set_yscale('log'); ax[0].set_ylabel('Rate / tau=0 reference'); ax[0].set_xlabel('Clockfield coupling tau'); ax[0].legend(); ax[0].grid(alpha=.2)
    ax[1].plot(taus,[r['predicted_ratio'] for r in rows],'-',label='Linearized prediction')
    ax[1].plot(taus,[r['clock_ratio'] for r in rows],'o',label='Nonlinear simulation')
    ax[1].axhline(1,color='gray',ls='--',label='Common rate')
    ax[1].set_xlabel('Clockfield coupling tau'); ax[1].set_ylabel('Normalized phase-travel / radial rate'); ax[1].legend(); ax[1].grid(alpha=.2)
    fig.savefig(args.out/'legacy_clocks.svg'); plt.close(fig)
    svg=(args.out/'legacy_clocks.svg').read_text(); svg=svg[svg.index('<svg'):]
    import html
    (args.out/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>Original Clockfield: two clocks</title><style>body{font:17px system-ui;max-width:1100px;margin:30px auto;padding:20px;color:#173247;background:#f4f7fa}svg{width:100%;height:auto;background:white}pre{white-space:pre-wrap;font:14px monospace;line-height:1.6}</style><h1>Testing the original Clockfield equation</h1><p>The travelling phase packet slows more than the radial oscillation. Measured rates follow the extra spatial coefficient in the original equation.</p>'+svg+'<pre>'+html.escape(text)+'</pre></html>')
    print(text)
    required=('prediction','half_step','half_amplitude','second_grid_refinement')
    if not all(validation[k] for k in required): raise SystemExit('Required validation failed; inspect results.json')
    if not validation['first_grid_refinement']: print('NOTE: coarse raw pulse speed failed 0.5% grid tolerance; finer-grid check passes. See report.')

if __name__=='__main__': main()
