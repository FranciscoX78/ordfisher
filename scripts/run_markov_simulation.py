#!/usr/bin/env python3
"""Independently simulate the two stationary Markov chains from rational Q.

Uses new fixed seeds, so Monte Carlo means are not required to exactly match
historical figures. Their exact 4-pattern laws and spectra DO match.
"""
import argparse,csv,json,platform
from pathlib import Path
from itertools import combinations
import numpy as np
from ordfisher import __version__
from ordfisher.markov import verify_exact_certificate,simulate_markov,theoretical_probabilities
from ordfisher.core import (ordinal_probabilities,spectra_from_probabilities,
                             total_variation,cayley_geometry)

def run(*,out:Path,seed:int=20261105,realizations:int=30,windows:int=32768):
    if realizations<2 or windows<10:raise ValueError('Need at least two repetitions and ten windows')
    verify_exact_certificate()
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    rng_seeds=np.random.SeedSequence(seed).spawn(realizations*2)
    rows=[];est={"plus":[],"minus":[]}
    for which,name in enumerate(('plus','minus')):
        for k in range(realizations):
            rng=np.random.default_rng(rng_seeds[which*realizations+k])
            this_seed=int(rng.integers(0,2**62,dtype=np.int64))
            series=simulate_markov(name,windows+3,seed=this_seed)
            p,diag=ordinal_probabilities(series,m=4,tau=1,return_diagnostics=True)
            r=spectra_from_probabilities(p,m=4)
            est[name].append(r)
            rows.append({'process':name,'realization':k,'seed':this_seed,
                         'n_windows':diag['n_windows'],
                         'normalized_entropy':r.shannon_entropy,
                         'fisher_energy':r.fisher_energy,
                         'n_missing_patterns':int(np.count_nonzero(p==0))})
        print(f'{name}: simulated {realizations} independent trajectories',flush=True)
    def tvpairs(a,b,same=False,attr='ofs_profile'):
        if same:
            return np.array([total_variation(getattr(a[i],attr),getattr(a[j],attr))
                        for i,j in combinations(range(len(a)),2)])
        return np.array([total_variation(getattr(x,attr),getattr(y,attr)) for x in a for y in b])
    summary={"seed":seed,"ordfisher_version":__version__,"realizations_per_process":realizations,
             "overlapping_windows_per_trajectory":windows,
             "continuous_observations_per_trajectory":windows+3,
             "new_exact_certificate_not_original":True,
             "new_simulation_seed_not_original":True}
    for prop in ('ofs_profile','ops_profile'):
        prefix='OFS' if prop=='ofs_profile' else 'OPS'
        for label,x in [('between',tvpairs(est['plus'],est['minus'],attr=prop)),
                        ('same_plus',tvpairs(est['plus'],est['plus'],same=True,attr=prop)),
                        ('same_minus',tvpairs(est['minus'],est['minus'],same=True,attr=prop))]:
            summary[prefix+'_'+label+'_mean_tv']=float(x.mean())
            summary[prefix+'_'+label+'_n_pairs']=len(x)
    a=spectra_from_probabilities(theoretical_probabilities(1),m=4)
    b=spectra_from_probabilities(theoretical_probabilities(-1),m=4)
    summary['exact_between_ofs_tv']=total_variation(a.ofs_profile,b.ofs_profile)
    summary['exact_between_ops_tv']=total_variation(a.ops_profile,b.ops_profile)
    with (out/'markov_realizations.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]))
        w.writeheader();w.writerows(rows)
    (out/'markov_simulation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2),flush=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f=a.frequencies[1:]
    fig,ax=plt.subplots(1,2,figsize=(10.6,4.0))
    ax[0].plot(f,a.ops_profile,marker='o',label='Process A')
    ax[0].plot(f,b.ops_profile,marker='s',linestyle='--',label='Process B')
    ax[0].set_title('OPS: exact matching spectra')
    ax[1].plot(f,a.ofs_profile,marker='o',label='Process A')
    ax[1].plot(f,b.ofs_profile,marker='s',linestyle='--',label='Process B')
    ax[1].set_title(f'OFS: TV = {summary["exact_between_ofs_tv"]:.4f}')
    for xx in ax:
        xx.set(xlabel='Cayley graph frequency',ylabel='Normalized spectral power')
        xx.legend()
    fig.tight_layout();fig.savefig(out/'markov_spectra.png',dpi=240);fig.savefig(out/'markov_spectra.pdf');plt.close(fig)
    return summary

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=Path('results/markov'))
    ap.add_argument('--seed',type=int,default=20261105)
    ap.add_argument('--realizations',type=int,default=30)
    ap.add_argument('--windows',type=int,default=32768)
    run(**vars(ap.parse_args()))
