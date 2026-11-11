#!/usr/bin/env python3
"""Regenerate the logistic-map OPS/OFS experiment in the Chaos manuscript.

IMPORTANT: the seed is a NEW, documented seed. The original set of initial
conditions was not recovered, so do not expect bitwise reproduction of the
previously submitted figure; the experiment's stated protocol is reproduced.
"""
from __future__ import annotations
import argparse, csv, json, platform
from pathlib import Path
import numpy as np
import scipy
from ordfisher.core import ordinal_probabilities, spectra_from_probabilities
from ordfisher.models import logistic_series
from ordfisher import __version__


def write_csv(path,rows,columns):
    with Path(path).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=columns)
        w.writeheader();w.writerows(rows)


def run(*,out:Path,seed:int=20261105,repetitions:int=10,
        retained:int=120000,burn:int=5000,
        r_min:float=3.5,r_max:float=4.0,r_step:float=.01,
        m:int=4,tau:int=1,use_numba:bool=True):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    n_r=round((r_max-r_min)/r_step)+1
    if n_r<1 or abs((r_min+(n_r-1)*r_step)-r_max)>1e-8:
        raise ValueError('r_step must exactly span [r_min,r_max]')
    parameters=np.linspace(r_min,r_max,n_r)
    children=np.random.SeedSequence(seed).spawn(len(parameters))
    runs=[];summary=[]
    for kk,r in enumerate(parameters):
        generator=np.random.default_rng(children[kk])
        x0s=generator.uniform(0,1,size=repetitions)
        sub=[]
        for j,x0 in enumerate(x0s):
            series=logistic_series(float(r),float(x0),burn=burn,retained=retained,use_numba=use_numba)
            p,diag=ordinal_probabilities(series,m=m,tau=tau,engine='numpy',ties='stable',return_diagnostics=True)
            result=spectra_from_probabilities(p,m=m)
            f_cent,p_cent=result.centroids()[1],result.centroids()[0]
            row={"r":float(r),"realization":j,"x0":float(x0),
                 "d_spec":result.ops_ofs_tv(),
                 "centroid_difference":abs(f_cent-p_cent),
                 "fisher_energy":result.fisher_energy,
                 "shannon_entropy":result.shannon_entropy,
                 "missing_patterns":int(np.count_nonzero(p==0)),
                 "n_tied_windows":diag['n_tied_windows']}
            runs.append(row);sub.append(row)
        def summarise(k):
            data=np.array([z[k] for z in sub],dtype=float)
            return float(data.mean()),float(data.std(ddof=1)) if len(data)>1 else 0.0
        d_m,d_sd=summarise('d_spec')
        c_m,c_sd=summarise('centroid_difference')
        i_m,i_sd=summarise('fisher_energy')
        summary.append({"r":float(r),"n_realizations":len(sub),
                        "d_spec_mean":d_m,"d_spec_sd":d_sd,
                        "centroid_difference_mean":c_m,"centroid_difference_sd":c_sd,
                        "fisher_energy_mean":i_m,"fisher_energy_sd":i_sd,
                        "total_tied_windows":sum(z['n_tied_windows'] for z in sub)})
        print(f'r={r:.2f}: OPS/OFS TV = {d_m:.6f} +/- {d_sd:.6f}',flush=True)
    write_csv(out/'logistic_runs.csv',runs,list(runs[0]))
    write_csv(out/'logistic_summary.csv',summary,list(summary[0]))
    metadata={"experiment":"Logistic-map OPS/OFS comparison", "seed":seed,
              "ordfisher_version":__version__,
              "new_seed_not_original_figures":True,
              "n_realizations":repetitions,"burn":burn,"retained":retained,
              "m":m,"tau":tau,"r_min":r_min,"r_max":r_max,"r_step":r_step,
              "n_r":n_r,"tie_policy":"stable (earlier time index wins exact ties)",
              "sd":"sample standard deviation (ddof=1)",
              "numpy":np.__version__,"scipy":scipy.__version__,"python":platform.python_version(),
              "note":"Analytic non-overlapping confidence intervals are not claimed: ordinal windows overlap."}
    (out/'logistic_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        x=np.array([e['r'] for e in summary])
        y=np.array([e['d_spec_mean'] for e in summary]);sy=np.array([e['d_spec_sd'] for e in summary])
        z=np.array([e['centroid_difference_mean'] for e in summary]);sz=np.array([e['centroid_difference_sd'] for e in summary])
        fig,ax=plt.subplots(figsize=(8.8,4.7))
        ax.plot(x,y,label=r'OPS/OFS $d_{\mathrm{spec}}$',linewidth=1.45)
        ax.fill_between(x,y-sy,y+sy,alpha=.18)
        ax.plot(x,z,linestyle='--',label=r'$|\Delta\nu_c|$',linewidth=1.3)
        ax.fill_between(x,z-sz,z+sz,alpha=.13)
        ax.set(xlabel='Logistic parameter r',ylabel='Profile comparison',title='Logistic map: m=4, tau=1')
        ax.legend();fig.tight_layout()
        fig.savefig(out/'logistic_figure.png',dpi=240)
        fig.savefig(out/'logistic_figure.pdf')
        plt.close(fig)
    except ImportError:
        print('matplotlib not installed; CSV and metadata have been saved')
    chosen={f'{r:.2f}':s for r in (3.58,3.70,3.93) for s in summary if np.isclose(s['r'],r)}
    print('Selected values:',json.dumps(chosen,indent=2))
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,default=Path('results/logistic'))
    ap.add_argument('--seed',type=int,default=20261105)
    ap.add_argument('--repetitions',type=int,default=10)
    ap.add_argument('--retained',type=int,default=120000)
    ap.add_argument('--burn',type=int,default=5000)
    ap.add_argument('--r-min',type=float,default=3.5)
    ap.add_argument('--r-max',type=float,default=4.0)
    ap.add_argument('--r-step',type=float,default=.01)
    ap.add_argument('--no-numba',action='store_true')
    args=ap.parse_args()
    run(out=args.out,seed=args.seed,repetitions=args.repetitions,
        retained=args.retained,burn=args.burn,r_min=args.r_min,r_max=args.r_max,
        r_step=args.r_step,use_numba=not args.no_numba)
