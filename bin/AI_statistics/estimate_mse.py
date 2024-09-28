from scipy.special import logit, expit
from scipy import stats as st
import numpy as np


def estimate_w_null(x, n, B):
    b = np.log(B)
    delta = b * (n - 2 * x)
    return expit(-delta)


def expectation_vectorized(func, n, B, p):
    b = np.log(B)
    p1 = expit(logit(p) + b)
    x = np.arange(0, n + 1)
    
    pmf_values = st.binom.pmf(x[:, None], n, p1[None, :])
    func_values = func(x[:, None], n, B )
    
    expectations = np.sum(pmf_values * func_values, axis=0)
    return expectations
        

def effect_size_estimate(x, n, B, w=None):
    if w is None:
        w = estimate_w_null(x, n, B)
    x, n, B = map(np.asarray, [x, n, B])

    b = np.log(B)
    p = x / n
    return w * (logit(p) - b) + (1 - w) * (logit(p) + b)
    

def es_var_vectorized(n, B, p, w=None):
    kwargs = dict(n=n, B=B, p=p, w=w)
    exp = expectation_vectorized(effect_size_estimate, **kwargs)
    def x2(*args, **kwargs):
        return effect_size_estimate(*args, **kwargs) ** 2
    return expectation_vectorized(x2, **kwargs) - exp ** 2


def calc_variance(n, B, n_points=101):
    ess = np.linspace(0, 1, n_points)
    yvals = expectation_vectorized(effect_size_estimate, n, B, ess)
    vars = es_var_vectorized(n, B, ess)
    return (yvals - ess) ** 2 + vars
