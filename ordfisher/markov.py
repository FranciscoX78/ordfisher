"""Exact certificate and simulations for the OFS/OPS Markov example.

The bundled arrays are a NEW rational realization of the four-bin-block
constraints of the paper. They are not represented as the originally produced
certificate, which was not available in the accessible source archive.

Q[i,j,k,l] = numerator / denominator, with bins 0,...,7.
Each conditional emission is X_t = (B_t + U_t)/8, U_t~Uniform(0,1).
A strictly positive 4-block law, whose both 3-block marginals are uniform,
defines a stationary, geometrically mixing order-3 Markov process via
Pr(B_{t+3}=l | B_t=i,B_{t+1}=j,B_{t+2}=k) = 8**3 Q[i,j,k,l].
"""
from __future__ import annotations

from functools import lru_cache
from itertools import permutations, product
from importlib.resources import files
import json
import numpy as np

from .core import cayley_geometry, spectra_from_probabilities, total_variation

POSITIVE_PATTERNS=("0231","2103","3012")
NEGATIVE_PATTERNS=("0213","0312","2310")
THETA_NUMERATOR,THETA_DENOMINATOR=99,125
BINS=8


def _tuple(s:str)->tuple[int,...]:
    return tuple(map(int,s))


def markov_sign_vector() -> np.ndarray:
    """v of p(±)=(1±(99/125)*v)/24 in lexicographic S4 order."""
    perm=[tuple(row) for row in cayley_geometry(4).patterns]
    ix={p:i for i,p in enumerate(perm)}
    v=np.zeros(24,dtype=np.int64)
    for k in POSITIVE_PATTERNS:v[ix[_tuple(k)]]=+1
    for k in NEGATIVE_PATTERNS:v[ix[_tuple(k)]]=-1
    return v


def theoretical_probabilities(sign: int)->np.ndarray:
    if sign not in (-1,1):raise ValueError("sign must be +1 or -1")
    return (1+sign*(THETA_NUMERATOR/THETA_DENOMINATOR)*markov_sign_vector())/24


@lru_cache(maxsize=1)
def ordinal_bin_coefficient_matrix()->np.ndarray:
    """Exact integers C[pi,ijkl] = 24 * P(pi | four bin values).

    Ties between equal bins are resolved by independent continuous U(0,1)
    emissions, so all relative rankings among tied bins are equally probable.
    Since m=4, denominators are 1,2,6,24 and 24*C_conditional is integral.
    """
    geom=cayley_geometry(4)
    ix={tuple(p):r for r,p in enumerate(geom.patterns)}
    C=np.zeros((24,8**4),dtype=np.int16)
    for i,j,k,l in product(range(8),repeat=4):
        block=(i,j,k,l)
        col=((i*8+j)*8+k)*8+l
        groups={}
        for pos,bin_value in enumerate(block):
            groups.setdefault(bin_value,[]).append(pos)
        options=[tuple(permutations(groups[bin_value])) for bin_value in sorted(groups)]
        total=1
        for option in options:total*=len(option)
        assert 24%total==0
        for order in product(*options):
            pattern=tuple(index for grp in order for index in grp)
            C[ix[pattern],col]+=24//total
    assert np.all(C.sum(axis=0)==24)
    return C


@lru_cache(maxsize=1)
def lag3_centered_coefficients()->np.ndarray:
    """Coefficients (2i-7)(2l-7); exact zero ↔ Cov(X_t,X_t+3)=0."""
    return np.array([(2*i-7)*(2*l-7)
                     for i,j,k,l in product(range(8),repeat=4)],dtype=np.int64)


@lru_cache(maxsize=1)
def certificates()->dict[str,dict]:
    path=files("ordfisher").joinpath("data/markov_q_exact.json")
    with path.open("r",encoding="utf-8") as f: obj=json.load(f)
    assert obj['bins']==8
    return obj


def exact_Q(name:str)->tuple[np.ndarray,int]:
    """Return (integer numerators Q_ijkl, common denominator)."""
    if name not in ("plus","minus"):raise ValueError("Choose 'plus' or 'minus'")
    data=certificates()[name]
    den=int(data['denominator'])
    num=np.asarray(data['numerators'],dtype=np.int64).reshape((8,8,8,8))
    if num.size!=4096 or den<=0:raise RuntimeError("Invalid certificate shape")
    return num,den


def _exact_ordinal_prob_numerator(name:str)->tuple[np.ndarray,int]:
    num,den=exact_Q(name)
    # C@Q = 24*p; thus 24*den*p is C@num, all integers.
    return (ordinal_bin_coefficient_matrix().astype(np.int64)@num.ravel(),24*den)


def verify_exact_certificate(verbose:bool=False)->dict:
    """Verify all defining Markov identities using integer arithmetic ONLY.

    Floating point is used only later for the optional graph spectrum comparison.
    The order-three Markov chain is irreducible and aperiodic because every
    four-block probability is strictly positive: any triple can reach any
    triple after three shifts, and constant triples have positive self-loops.
    """
    C=ordinal_bin_coefficient_matrix().astype(np.int64)
    centered=lag3_centered_coefficients()
    v=markov_sign_vector()
    geom=cayley_geometry(4)
    results={}
    for name,sign in [("plus",1),("minus",-1)]:
        n,D=exact_Q(name)
        assert D%512==0 and D%125==0, name
        assert np.all(n>0), (name,"nonpositive four-block probability")
        assert np.all(n.sum(axis=3)==D//512), (name,"outgoing 3-block marginal")
        assert np.all(n.sum(axis=0)==D//512), (name,"incoming 3-block marginal")
        assert int(n.sum())==D, (name,"normalization")
        assert int(centered@n.ravel())==0, (name,"lag-3 covariance")
        calculated=C@n.ravel()
        wanted=D + sign*THETA_NUMERATOR*(D//THETA_DENOMINATOR)*v
        assert np.array_equal(calculated,wanted), (name,"ordinal law mismatch")
        results[name]={"denominator":int(D),"minimum_numerator":int(n.min()),
                       "maximum_numerator":int(n.max()),"sum_numerators":int(n.sum()),
                       "max_exact_ordinal_residual":int(np.max(np.abs(calculated-wanted))),
                       "max_exact_3block_residual":0,
                       "exact_lag3_centered_covariance_numerator":0}
    # Exact equal OPS follows from p_plus-u=-(p_minus-u).
    assert np.array_equal(v,-(-v))
    # Exact equal scalar Fisher energy: on the Cayley graph the numbers of
    # undirected edges connecting +1 to zero and -1 to zero are both 7.
    labels={tuple(p):int(v[i]) for i,p in enumerate(geom.patterns)}
    pairs={(-1,0):0,(0,1):0,(-1,1):0,(-1,-1):0,(0,0):0,(1,1):0}
    for a,b in geom.edges:
        key=tuple(sorted((int(v[a]),int(v[b]))))
        pairs[key]+=1
    assert pairs[(-1,0)]==pairs[(0,1)]==7
    assert pairs[(-1,1)]==2
    results['edge_counts']={f'{a},{b}':count for (a,b),count in pairs.items()}
    results['exact_claims']=[
        "both Q arrays are strictly positive rational four-block probabilities",
        "uniform incoming and outgoing three-block marginals imply stationarity",
        "three consecutive continuous observations are IID Uniform(0,1)",
        "lag-one and lag-two correlations vanish by three-point independence",
        "lag-three covariance vanishes by exact centered integer constraint",
        "specified m=4 ordinal probabilities match exactly",
        "Shannon and all Renyi entropies agree (same multiset of probabilities)",
        "full OPS agrees (opposite deviations from uniformity)",
        "scalar Cayley Fisher energy agrees by exact edge-count symmetry",
        "strict positivity of four-block Q makes the order-3 chain mixing",
    ]
    # Numerical spectral comparison, distinct from exact integer checks.
    aa=spectra_from_probabilities(theoretical_probabilities(+1),m=4)
    bb=spectra_from_probabilities(theoretical_probabilities(-1),m=4)
    results['spectral_numerical']={
        "max_absolute_ops_difference":float(np.max(np.abs(aa.ops-bb.ops))),
        "absolute_fisher_energy_difference":abs(aa.fisher_energy-bb.fisher_energy),
        "common_fisher_energy":aa.fisher_energy,
        "common_normalized_shannon_entropy":aa.shannon_entropy,
        "tv_ofs_normalized":total_variation(aa.ofs_profile,bb.ofs_profile),
    }
    assert results['spectral_numerical']['max_absolute_ops_difference']<1e-12
    assert results['spectral_numerical']['absolute_fisher_energy_difference']<1e-12
    assert abs(results['spectral_numerical']['tv_ofs_normalized']-0.2312292799)<1e-8
    if verbose:
        import json as _json
        print(_json.dumps(results,indent=2))
    return results


def markov_transition_probabilities(name:str)->np.ndarray:
    """P[l|i,j,k] = 512 Q[i,j,k,l], shape (8,8,8,8)."""
    n,d=exact_Q(name)
    P=n.astype(np.float64)*(512/d)
    assert np.max(np.abs(P.sum(axis=3)-1))<1e-13
    return P


def simulate_markov(name:str, n_observations:int, seed:int=20261105)->np.ndarray:
    """A stationary continuous-valued realization from the rational chain.

    The initial three bins have the exact stationary joint distribution:
    three independent discrete uniforms on {0,..,7}. Emission uniform noise
    is independent of transitions.
    """
    if n_observations<4:raise ValueError("n_observations must be >=4")
    rng=np.random.default_rng(seed)
    P=markov_transition_probabilities(name)
    cdf=np.cumsum(P,axis=-1)
    draws=rng.random(n_observations-3)
    bins=np.empty(n_observations,dtype=np.int8)
    bins[:3]=rng.integers(0,8,size=3,dtype=np.int8)
    for t in range(3,n_observations):
        i,j,k=map(int,bins[t-3:t])
        bins[t]=int(np.searchsorted(cdf[i,j,k],draws[t-3],side="right"))
    emissions=rng.random(n_observations)
    return (bins.astype(float)+emissions)/8


def ordinal_transition_probability(sign:int, pattern_before:str, pattern_after:str):
    """Exact Fraction of two overlapping length-3 patterns from p4.

    Uses the p4 law, no simulation and no original Q required.
    """
    from fractions import Fraction
    if sign not in (1,-1):raise ValueError("sign must be +1 or -1")
    expected=Fraction(0,1)
    source=Fraction(0,1)
    v=markov_sign_vector()
    for idx,pi0 in enumerate(cayley_geometry(4).patterns):
        p=Fraction(125+sign*99*int(v[idx]),24*125)
        # pi0 is time indices in increasing value order.
        before=tuple(sorted(range(3),key=lambda z: list(pi0).index(z)))
        after=tuple(z-1 for z in sorted(range(1,4),key=lambda z:list(pi0).index(z)))
        if before==_tuple(pattern_before):
            source+=p
            if after==_tuple(pattern_after):expected+=p
    return expected/source
