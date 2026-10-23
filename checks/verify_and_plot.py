"""Reproducible numerical checks and the manuscript's explanatory figures.

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
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from plot_growth_example import plot_growth_example

ROOT = Path(__file__).resolve().parents[1]

def constants(s):
    c = 1.0 / zeta(s, 1)
    A = gamma(1 - 1/s) * c**(1/s)
    B = (2**(1/s)-1) * A
    D = c / ((s-1)*(2-s)) if 1 < s < 2 else None
    return c, A, B, D

def third_coefficient(s):
    """E_s for internal vertices (negative sign for the leaf mean)."""
    if not 1 < s < 2:
        raise ValueError('The third power correction requires 1 < s < 2.')
    return (1/zeta(s,1))**(2/s)*gamma(2-2/s)/(2-s)

def coefficient_integral(s):
    """Independent quadrature of E_s; split at 1 and invert the tail."""
    low=quad(lambda y: -np.expm1(-y)/y if y else 1.0,0,1,
             weight='alg',wvar=(1-2/s,0),epsabs=1e-11)[0]
    high=quad(lambda x: -np.expm1(-1/x) if x else 1.0,0,1,
              weight='alg',wvar=(2/s-2,0),epsabs=1e-11)[0]
    return (1/zeta(s,1))**(2/s)*(low+high)/s

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
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'savefig.bbox':'tight'})
    s=np.linspace(1.27,1.98,500)
    star=1+1/np.sqrt(2);s0=(1+np.sqrt(3))/2
    fig,ax=plt.subplots(figsize=(3.45,2.45))
    ax.plot(s,2-s,color='#b75c23',lw=1.7,label=r'First bias: $2-s$')
    ax.plot(s,1/(2*s),color='#176b89',lw=2,label=r'Gaussian scale: $1/(2s)$')
    ax.plot(s,2/s-1,color='#82518c',lw=1.6,ls='-.',label=r'Second bias: $2/s-1$')
    ax.plot(s,np.maximum(0,1+1/s-s),color='#555555',lw=1.3,ls='--',
            label=r'New error: $\max(0,1+1/s-s)$')
    for x,label in ((s0,r'$s_0$'),(1.5,r'$3/2$'),(star,r'$s_\star$')):
        ax.plot([x,x],[0,1/(2*x)],color='#aaaaaa',ls=':',lw=.8)
        ax.scatter([x],[1/(2*x)],color='#202020',s=12,zorder=4)
        ax.text(x,-.055,label,fontsize=8,ha='center')
    ax.set(xlim=(1.27,1.98),ylim=(-.075,.8),xlabel=r'Power $s$',ylabel=r'Exponent of $n$')
    ax.legend(frameon=False,loc='upper right',fontsize=7)
    ax.grid(axis='y',alpha=.15)
    fig.tight_layout(pad=.4)
    (ROOT/'figures').mkdir(exist_ok=True)
    fig.savefig(ROOT/'figures/centering_scales.pdf')
    fig.savefig(ROOT/'figures/centering_scales.png',dpi=220)
    plt.close(fig)

def main():
    rng = np.random.default_rng(20261020)
    report = {'seed':20261020, 'enumeration':[], 'coupling_trials':0,
              'mean_expansion':[], 'third_term':[], 'boundary_s2':[],
              'coefficient_quadrature':[], 'notes':[
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
    for s in (1.3,1.45,1.6,1.8):
        c,A,B,D=constants(s)
        E=third_coefficient(s)
        integral=coefficient_integral(s)
        assert abs(integral-E)<1e-10*max(1,E)
        report['coefficient_quadrature'].append(dict(s=s,closed_form=E,integral=integral))
        for n in (10000,100000,1000000):
            mu,error=mean_series(n,s)
            residual=mu-A*n**(1/s)+D*n**(2-s)
            report['third_term'].append(dict(n=n,s=s,internal_mean=mu,
                Q_n=residual/n**(2/s-1),E_s=E,
                error_after_third_term=residual-E*n**(2/s-1),
                proved_error_scale=1+n**(1+1/s-s),series_truncation_bound=error))
    c,A,B,_=constants(2)
    for n in (10000,100000,1000000):
        mu,error=mean_series(n,2)
        report['boundary_s2'].append(dict(n=n,internal_mean=mu,
            residual=mu-A*np.sqrt(n)+c/2*np.log(n),series_truncation_bound=error))
    plot_growth_example()
    plot_scales()
    out=ROOT/'checks'/'verification_results.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'enumerated_cases':len(report['enumeration']),
                      'coupling_trials':report['coupling_trials'],
                      'large_n_cases':len(report['mean_expansion']),
                      'third_term_cases':len(report['third_term']),
                      'boundary_cases':len(report['boundary_s2']),
                      'coefficient_integrals':len(report['coefficient_quadrature']),
                      'all_assertions_passed':True,'results':str(out)},indent=2))

if __name__ == '__main__':
    main()
