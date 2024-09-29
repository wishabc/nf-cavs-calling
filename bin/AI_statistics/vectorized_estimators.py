from scipy.special import logit, expit
from scipy import stats as st
import numpy as np

# Vectorized functions for effect size estimates, expectation, variance and MSE


def estimate_w_null(x, n, B):
    b = np.log(B)
    delta = b * (n - 2 * x)
    return expit(-delta)


def mode1_expectation_vectorized(func, n, B, p):
    b = np.log(B)
    p1 = expit(logit(p) + b)
    x = np.arange(0, n + 1)
    
    pmf_values = st.binom.pmf(x[:, None], n, p1[None, :])
    func_values = func(x[:, None], n, B )
    
    expectations = np.sum(pmf_values * func_values, axis=0)
    return expectations


def es_estimate_vectorized(x, n, B, w=None):
    if w is None:
        w = estimate_w_null(x, n, B)
    x, n, B = map(np.asarray, [x, n, B])

    b = np.log(B)
    p = x / n
    return w * (logit(p) - b) + (1 - w) * (logit(p) + b)
    

def es_variance_vectorized(n, B, p, w=None):
    assert np.all(p >= 0) and np.all(p <= 1), "p must be in [0, 1]"
    kwargs = dict(n=n, B=B, p=p, w=w)
    exp = mode1_expectation_vectorized(es_estimate_vectorized, **kwargs)
    def es_estimate_squared(*args, **kwargs):
        return es_estimate_vectorized(*args, **kwargs) ** 2
    return mode1_expectation_vectorized(es_estimate_squared, **kwargs) - exp ** 2


def calc_variance(n, B, n_points=101):
    ess = np.linspace(0, 1, n_points)
    yvals = mode1_expectation_vectorized(es_estimate_vectorized, n, B, ess)
    vars = es_variance_vectorized(n, B, ess)
    return (yvals - ess) ** 2 + vars

def logit_es(p, d=1/128):
    return np.log2(p + d) - np.log2(1 - p + d)

def inv_logit_es(es, d=1/128):
    s = 2 ** es
    return (1 + (s - 1) * (d + 1)) / (s + 1)
