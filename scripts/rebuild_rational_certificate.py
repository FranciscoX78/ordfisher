#!/usr/bin/env python3
"""Construct NEW strictly positive rational Q+ and Q- for Markov example.

Method: feasibility LP in floats -> rational rounding preserving exact three-block
marginals -> exact rational correction along four-cycle transportation moves.
The final saved certificate is verified with INTEGER ARITHMETIC, not by solver
floating-point residuals. scipy and sympy are needed for reconstruction.

The distributed fixed certificate is sufficient for the manuscript; running this
script is optional. Different SciPy/HiGHS versions may choose different rational
certificates that satisfy the same exact constraints.
"""
from __future__ import annotations
from itertools import product
import json
from math import lcm
from pathlib import Path
import numpy as np
from scipy.sparse import lil_matrix,csr_matrix
from scipy.linalg import qr
from scipy.optimize import linprog
from sympy import Matrix
from ordfisher.markov import (ordinal_bin_coefficient_matrix,
                              lag3_centered_coefficients,markov_sign_vector)


def rebuild(output:Path):
    C=ordinal_bin_coefficient_matrix().astype(np.int64)
    cov=lag3_centered_coefficients()
    v=markov_sign_vector()
    n=4096
    A=lil_matrix((1049,n),dtype=float)
    for i,j,k,l in product(range(8),repeat=4):
        col=((i*8+j)*8+k)*8+l
        A[(i*8+j)*8+k,col]=1
        A[512+(j*8+k)*8+l,col]=1
    A[1024:1048,:]=csr_matrix(C/24)
    A[1048,:]=csr_matrix(cov.reshape((1,-1)).astype(float))
    A=A.tocsr()
    base_target=np.zeros(1049)
    base_target[:1024]=1/512
    # 2x2 transportation moves at fixed intermediate bins (j,k): preserve
    # all incoming and outgoing 3-block marginals EXACTLY.
    cycles=[]
    directions=np.empty((25,64*49),dtype=np.int64)
    c=0
    for j,k in product(range(8),repeat=2):
        for i,l in product(range(7),repeat=2):
            a=(((i*8+j)*8+k)*8+l)
            b=(((i*8+j)*8+k)*8+7)
            d=(((7*8+j)*8+k)*8+l)
            e=(((7*8+j)*8+k)*8+7)
            directions[:24,c]=C[:,a]-C[:,b]-C[:,d]+C[:,e]
            directions[24,c]=int(cov[a]-cov[b]-cov[d]+cov[e])
            cycles.append((a,b,d,e));c+=1
    R,piv=qr(directions.astype(float),pivoting=True,mode='economic')[1:3]
    rank=int(np.sum(np.abs(np.diag(R))>1e-7))
    cols=piv[:rank].tolist()
    R2,pivrows=qr(directions[:,cols].T.astype(float),pivoting=True,mode='economic')[1:3]
    rows=pivrows[:rank].tolist()
    invertible=Matrix(directions[np.ix_(rows,cols)].tolist())
    assert invertible.rank()==rank
    D=10**10
    out={"description":"New exact rational four-block laws; not the original unpublished arrays.",
         "bins":8,"length":4,
         "axes":"Q[i,j,k,l] lexicographic flattening, index=512*i+64*j+8*k+l"}
    for sign,name in [(+1,"plus"),(-1,"minus")]:
        target=base_target.copy()
        target[1024:1048]=(1+sign*99/125*v)/24
        res=linprog(np.zeros(n),A_eq=A,b_eq=target,
                    bounds=(2e-6,None),method='highs')
        if not res.success:raise RuntimeError(f"No feasible LP for {name}: {res.message}")
        # Rounding the upper-left 7x7 cells of every 8x8 transportation
        # matrix leaves exactly solvable marginal constraints for last row/col.
        floating=res.x.reshape((8,8,8,8))
        numer=np.zeros((8,8,8,8),dtype=np.int64)
        for j,k in product(range(8),repeat=2):
            mat=numer[:,j,k,:]
            mat[:7,:7]=np.rint(floating[:7,j,k,:7]*D).astype(np.int64)
            mat[:7,7]=D//512-mat[:7,:7].sum(axis=1)
            mat[7,:7]=D//512-mat[:7,:7].sum(axis=0)
            mat[7,7]=D//512-mat[7,:7].sum()
            assert np.all(mat.sum(axis=1)==D//512)
            assert np.all(mat.sum(axis=0)==D//512)
        assert numer.min()>0
        wanted=D+sign*99*(D//125)*v
        rhs=np.append(wanted-C@numer.ravel(),-cov@numer.ravel())
        coeff=invertible.inv()*Matrix([int(rhs[k]) for k in rows])
        multiplier=1
        for z in coeff:multiplier=lcm(multiplier,int(z.q))
        adjustment=np.zeros(n,dtype=np.int64)
        for j,weight in zip(cols,coeff):
            x=int(weight*multiplier)
            a,b,d,e=cycles[j]
            adjustment[a]+=x; adjustment[b]-=x
            adjustment[d]-=x; adjustment[e]+=x
        num_final=numer.ravel()*multiplier+adjustment
        den_final=int(D*multiplier)
        if np.min(num_final)<=0:raise RuntimeError('Rational correction broke positivity')
        if not np.array_equal(C@num_final,den_final+sign*99*(den_final//125)*v):
            raise RuntimeError('Ordinal law not exactly satisfied')
        if int(cov@num_final)!=0:raise RuntimeError('lag3 covariance mismatch')
        mm=num_final.reshape((8,8,8,8))
        if not (np.all(mm.sum(axis=3)==den_final//512) and
                np.all(mm.sum(axis=0)==den_final//512)):
            raise RuntimeError('stationary three-block marginals failed')
        out[name]={"denominator":den_final,"numerators":list(map(int,num_final))}
        print(f'{name}: common denominator {den_final}; minimum numerator {int(num_final.min())}')
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(out,separators=(',',':'))+'\n',encoding='utf-8')
    print('Saved:',output)
    return output

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--out',type=Path,default=Path('results/rebuilt_markov_q_exact.json'))
    args=p.parse_args()
    rebuild(args.out)
