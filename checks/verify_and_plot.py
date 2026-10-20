"""Reproducible numerical checks and the manuscript's analytic scale plot.

Run from any directory: python checks/verify_and_plot.py
Dependencies: Python 3, numpy, scipy, matplotlib.
Numerical checks support, but do not replace, the proofs.
"""
from pathlib import Path
import itertools
import json
import math
from collections import Counter
import numpy as np
from scipy.special import gamma, zeta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]

def constants(s):
    c = 1.0 / zeta(s, 1)
    A = gamma(1 - 1/s) * c**(1/s)
    B = (2**(1/s)-1) * A
    D = c / ((s-1)*(2-s)) if 1 < s < 2 else None
    return c, A, B, D

def mean_direct(n, s):
    w = np.arange(1, n, dtype=float)**(-s)
    H = np.cumsum(w)
    return 1.0 + sum(-np.expm1(np.log1p(-w[i]/H[i:]).sum())
                     for i in range(1, n-1))

def mean_series(n, s, order=24):
    """Product formula evaluated by a log series, O(order*n).

    The error returned bounds log-series truncation in exact arithmetic;
    it is NOT a rigorous floating-point interval.
    """
    w = np.arange(1, n, dtype=float)**(-s)
    H = np.cumsum(w)
    logq = np.zeros(n-1)
    for k in range(1, order+1):
        suffix = np.cumsum((H**(-k))[::-1])[::-1]
        logq -= (w**k) * suffix / k
    q = np.exp(logq)
    q[0] = 0.0  # root certainly has child 2
    k = order+1
    suffix = np.cumsum((H**(-k))[::-1])[::-1]
    tail = w[1:]**k * suffix[1:] / (k*(1-w[1:]/H[1:]))
    truncation_bound = float(np.sum(q[1:]*np.minimum(tail, 1.0)))
    return float(-np.expm1(logq[1:]).sum()+1), truncation_bound

def enumerate_moments(n, s):
    w = np.arange(1, n, dtype=float)**(-s)
    H = np.cumsum(w)
    total = mean = second = 0.0
    for parents in itertools.product(*(range(t) for t in range(1,n))):
        pr = math.prod(w[i]/H[t] for t,i in enumerate(parents))
        count = len(set(parents))
        total += pr
        mean += pr*count
        second += pr*count*count
    return total, mean, second-mean*mean

def coupling_check(n, s, rng):
    t = np.arange(2,n+1)
    X = rng.zipf(s, n-1)
    w = np.arange(1,n,dtype=float)**(-s)
    H = np.cumsum(w)
    Y = np.searchsorted(H, rng.random(n-1)*H)+1
    bad = X >= t
    P = np.where(bad,Y,X)
    M = int(bad.sum())
    G = set(X[~bad]); B = set(X[bad]); C = set(Y[bad])
    K = len(set(X)); I = len(set(P))
    R = M-len(B-G)+len(C-G)
    assert K-I == M-R and R >= 0 and abs(K-I) <= M
    cx, cp = Counter(X), Counter(P)
    energy = sum((int(cx[x]==1)-int(cp[p]==1))**2 for x,p in zip(X,P))
    assert energy <= 3*M
    return energy, M

def plot_scales():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'savefig.bbox':'tight'})
    s = np.linspace(1.28,1.98,500)
    star = 1+1/np.sqrt(2)
    fig, ax = plt.subplots(figsize=(6.3,2.3))
    ax.axvspan(1.28,star,color='#e6a468',alpha=.15,zorder=0)
    ax.axvspan(star,1.98,color='#4b94b3',alpha=.10,zorder=0)
    ax.plot(s,2-s,color='#b75c23',lw=2,label=r'Birth restriction: $2-s$')
    ax.plot(s,1/(2*s),color='#176b89',lw=2,label=r'Gaussian scale: $1/(2s)$')
    ax.plot(s,2/s-1,color='#777777',lw=1.4,ls='--',label=r'Remainder bound: $2/s-1$')
    ax.axvline(star,color='#444444',ls=':',lw=1)
    ax.scatter([star],[1/(2*star)],color='#202020',s=20,zorder=4)
    ax.text(star+.012,.70,r'$s_*=1+1/\sqrt{2}$',fontsize=9)
    ax.set(xlim=(1.28,1.98),ylim=(.04,.77),xlabel=r'Power $s$',ylabel=r'Exponent of $n$')
    ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(.985,.77),fontsize=8)
    ax.grid(axis='y',alpha=.15)
    fig.tight_layout(pad=.5)
    (ROOT/'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures'/'centering_scales.pdf')
    fig.savefig(ROOT/'figures'/'centering_scales.png',dpi=200)
    plt.close(fig)

def main():
    rng = np.random.default_rng(20261020)
    report = {'seed':20261020, 'enumeration':[], 'coupling_trials':0,
              'mean_expansion':[], 'notes':[
                  'Series error estimates exclude floating-point roundoff.',
                  'Mean ratios may converge slowly; finite-n checks are not proofs.',
                  'The figure plots proved exponents, not simulation estimates.']}
    for s in (1.25,1.6,1+1/np.sqrt(2),2.5):
        for n in (4,7):
            mass,mu,var = enumerate_moments(n,s)
            direct = mean_direct(n,s)
            series,error = mean_series(n,s)
            assert abs(mass-1)<1e-12
            assert abs(mu-direct)<1e-11
            assert abs(series-direct)<1e-10
            report['enumeration'].append(dict(n=n,s=s,mean=mu,variance=var,
                                               direct_mean=direct,series_mean=series))
        for _ in range(250):
            coupling_check(300,s,rng)
            report['coupling_trials'] += 1
    for s in (1.5,1.6,1+1/np.sqrt(2),1.85):
        c,A,B,D = constants(s)
        for n in (10000,100000,1000000):
            mu,error = mean_series(n,s)
            ratio = (A*n**(1/s)-mu)/(D*n**(2-s))
            report['mean_expansion'].append(dict(n=n,s=s,mean=mu,
                bias_ratio=ratio,series_truncation_bound=error,
                gaussian_scale=math.sqrt(B)*n**(1/(2*s))))
    plot_scales()
    out=ROOT/'checks'/'verification_results.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'enumerated_cases':len(report['enumeration']),
                      'coupling_trials':report['coupling_trials'],
                      'large_n_cases':len(report['mean_expansion']),
                      'all_assertions_passed':True,'results':str(out)},indent=2))

if __name__ == '__main__':
    main()
