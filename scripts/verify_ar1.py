#!/usr/bin/env python3
"""Exact Gaussian-orthant m=4 AR(1) computations and plot for appendix."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from ordfisher.models import ar1_gaussian_ordinal_probabilities,equal_fisher_ar1_pair
from ordfisher.core import spectra_from_probabilities,check_spectral_identities

def run(out:Path):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    pair=equal_fisher_ar1_pair()
    assert abs(pair['phi_positive']-0.6536965295755)<1e-9
    assert abs(pair['tv_fisher_energy_profiles']-0.41275842196)<1e-9
    assert abs(pair['tv_normalized_ofs_profiles']-0.32230089999)<1e-9
    (out/'ar1_equal_energy_pair.json').write_text(json.dumps(pair,indent=2)+'\n')
    rows=[]
    for j in range(-99,100):
        phi=j/100
        p=ar1_gaussian_ordinal_probabilities(phi)
        check_spectral_identities(p,m=4)
        z=spectra_from_probabilities(p,m=4)
        rows.append({'phi':phi,'fisher_energy':z.fisher_energy,
                     'd_spec':z.ops_ofs_tv() if z.ofs_profile is not None else '',
                     'shannon_entropy':z.shannon_entropy})
    with (out/'ar1_sweep.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    x=[r['phi'] for r in rows]
    y=[np.nan if r['d_spec']=='' else r['d_spec'] for r in rows]
    fig,ax=plt.subplots(figsize=(7.5,4.5))
    ax.plot(x,y,label='OFS/OPS TV distance')
    ax.set(xlabel='Gaussian AR(1) parameter phi',ylabel='Normalized-profile distance')
    fig.tight_layout();fig.savefig(out/'ar1_sweep.png',dpi=240);fig.savefig(out/'ar1_sweep.pdf')
    plt.close(fig)
    print(json.dumps(pair,indent=2))
    return pair

if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--out',type=Path,default=Path('results/ar1'))
    run(**vars(ap.parse_args()))
