#!/usr/bin/env python3
"""Verify the rational Markov certificate with integer arithmetic.

The certificate is a newly constructed exact realization of the manuscript's
constraints, not a recovered copy of the original unpublished arrays.
"""
import argparse, json, csv
from pathlib import Path
from ordfisher.markov import (verify_exact_certificate,ordinal_transition_probability,
                               exact_Q, theoretical_probabilities)
from ordfisher.core import spectra_from_probabilities,check_spectral_identities,cayley_geometry

def run(out=None):
    result=verify_exact_certificate(verbose=True)
    transitions={"plus":str(ordinal_transition_probability(1,'021','120')),
                 "minus":str(ordinal_transition_probability(-1,'021','120'))}
    assert transitions=={"plus":"56/125","minus":"13/250"}
    result['exact_ordinal_transition_021_to_120']=transitions
    for sign in (1,-1):
        check_spectral_identities(theoretical_probabilities(sign),m=4)
    if out is not None:
        dest=Path(out);dest.parent.mkdir(parents=True,exist_ok=True)
        dest.write_text(json.dumps(result,indent=2,default=float)+'\n',encoding='utf-8')
        graph=cayley_geometry(4)
        a=spectra_from_probabilities(theoretical_probabilities(+1),m=4)
        b=spectra_from_probabilities(theoretical_probabilities(-1),m=4)
        with (dest.parent/'markov_exact_spectra.csv').open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['frequency','OPS_plus','OPS_minus','OFS_plus','OFS_minus',
                             'normalized_OPS_plus','normalized_OPS_minus',
                             'normalized_OFS_plus','normalized_OFS_minus'])
            for i,v in enumerate(a.frequencies):
                pp=lambda r,attr: 0.0 if i==0 else float(getattr(r,attr)[i-1])
                writer.writerow([v,a.ops[i],b.ops[i],a.ofs[i],b.ofs[i],
                                 pp(a,'ops_profile'),pp(b,'ops_profile'),
                                 pp(a,'ofs_profile'),pp(b,'ofs_profile')])
        with (dest.parent/'markov_exact_ordinal_probabilities.csv').open('w',newline='') as f:
            writer=csv.writer(f)
            writer.writerow(['pattern_lexicographic','p_plus','p_minus'])
            for i,p in enumerate(graph.patterns):
                writer.writerow([''.join(map(str,p)),a.probabilities[i],b.probabilities[i]])
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=Path('results/markov_certificate_verification.json'))
    run(**vars(p.parse_args()))
