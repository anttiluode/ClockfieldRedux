"""Run deterministic controls and write a portable report. No network or GPU needed."""
import argparse
from dataclasses import asdict,replace
import json
from pathlib import Path
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from clockfield import Model,measure

CLOCKS = ('oscillation','pulse','exchange')

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('results'))
    args = parser.parse_args(); args.out.mkdir(parents=True,exist_ok=True)
    base=Model(); reference,trace=measure(base,return_traces=True)
    cases=[('reference',base)]
    for a in (.85,.65,.45): cases.append((f'imposed h={a}',replace(base,h=a)))
    for a in (.85,.65,.45): cases.append((f'stiffness c={a}',replace(base,c=a)))
    for s in (.15,.30,.45): cases.append((f'candidate s={s}',replace(base,mass=1+s,c=1-s)))
    rows=[]
    for name,m in cases:
        measured=measure(m); refined=measure(m,samples=4096)
        rates={k:measured[k]/reference[k] for k in CLOCKS}
        spread=max(rates.values())/min(rates.values())-1
        rows.append(dict(name=name,parameters=asdict(m),measured=measured,rates=rates,
                         spread=spread,agreement=bool(spread<.01),
                         sampling_error=max(abs(refined[k]/measured[k]-1) for k in CLOCKS)))
    grids={}
    for name,m in (cases[0],cases[-1]):
        values=[dict(n=n,pulse=measure(replace(m,n=n))['pulse']) for n in (256,512,1024,2048)]
        grids[name]=dict(values=values,last_relative_change=abs(values[-1]['pulse']/values[-2]['pulse']-1))
    checks=dict(energy=all(r['measured']['relative_energy_drift']<1e-10 for r in rows),
                sampling=all(r['sampling_error']<.002 for r in rows),
                spatial=all(g['last_relative_change']<.01 for g in grids.values()),
                imposed_control=all(r['agreement'] and max(abs(v-r['parameters']['h']) for v in r['rates'].values())<.002 for r in rows[1:4]),
                selective_control=all(not r['agreement'] for r in rows[4:7]))
    payload=dict(model='two-component periodic Hamiltonian lattice v1',python=platform.python_version(),numpy=np.__version__,rows=rows,spatial=grids,checks=checks)
    (args.out/'results.json').write_text(json.dumps(payload,indent=2,allow_nan=False)+'\n')
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,1,figsize=(11,9),layout='constrained')
    xx=np.arange(len(rows)); colors=['#147d92','#d87520','#823cad']
    for i,k in enumerate(CLOCKS): axes[0].bar(xx+(i-1)*.24,[r['rates'][k] for r in rows],.23,label=k,color=colors[i])
    axes[0].set_xticks(xx,[r['name'] for r in rows],rotation=28,ha='right')
    axes[0].set_ylabel('Rate / own reference'); axes[0].legend()
    axes[0].set_title('A common slowdown is reproduced only by the imposed control')
    for name,m in [cases[0],cases[2],cases[-1]]:
        _,tr=measure(m,return_traces=True)
        axes[1].plot(tr['tp'],tr['position']-tr['position'][0],label=name)
    axes[1].set_xlabel('External evolution parameter t'); axes[1].set_ylabel('Packet displacement (fixed lattice ruler)'); axes[1].legend(); axes[1].grid(alpha=.2)
    fig.savefig(args.out/'clocks.svg'); plt.close(fig)
    table='\n'.join('| '+r['name']+' | '+' | '.join(f"{r['rates'][k]:.4f}" for k in CLOCKS)+f" | {'yes' if r['agreement'] else 'no'} |" for r in rows)
    summary='# Measured clock comparison\n\nDeterministic linear toy model; no emergent-gravity claim.\n\n| Case | Oscillation | Pulse | Exchange | Agreement within 1% |\n|---|---:|---:|---:|---|\n'+table+'\n\nValidation: '+json.dumps(checks)+'.\n\nThe imposed Hamiltonian multiplier makes clocks agree by construction. All tested candidate backgrounds fail the common-rate criterion. Pulse rates use a fixed lattice ruler. See results.json for raw measurements and refinement checks.\n'
    (args.out/'SUMMARY.md').write_text(summary)
    htmlrows=''.join('<tr><td>'+r['name']+'</td>'+''.join(f'<td>{r["rates"][k]:.4f}</td>' for k in CLOCKS)+f'<td>{"yes" if r["agreement"] else "no"}</td></tr>' for r in rows)
    svg=(args.out/'clocks.svg').read_text(); svg=svg[svg.index('<svg'):]
    (args.out/'report.html').write_text('<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width"><title>ClockfieldRedux — clock comparison</title><style>body{font:17px system-ui;max-width:1100px;margin:40px auto;padding:20px;color:#193047;background:#f4f7fa}h1{font-size:40px}svg{width:100%;height:auto;background:white;border-radius:12px}table{border-collapse:collapse;width:100%;background:white}td,th{padding:12px;text-align:right;border-bottom:1px solid #ddd}td:first-child,th:first-child{text-align:left}p{line-height:1.6}</style><h1>Do the clocks agree?</h1><p>Three clocks, one conservative wave lattice. A universal slowdown is an experimental question. Here it occurs in the explicitly imposed control; the tested candidate backgrounds change the mechanisms differently.</p>'+svg+'<table><tr><th>Background</th><th>Oscillation</th><th>Pulse</th><th>Exchange</th><th>Within 1%?</th></tr>'+htmlrows+'</table><p>The pulse uses a fixed lattice ruler. This first experiment uses homogeneous parameters and linear waves; it contains no self-supported mass, dynamical ruler, entropy mechanism, or gravitational field.</p><p>Numerical validation: '+json.dumps(checks)+'</p></html>')
    print(summary)
    if not all(checks.values()): raise SystemExit('Numerical validation failed; inspect results.json')

if __name__=='__main__': main()
