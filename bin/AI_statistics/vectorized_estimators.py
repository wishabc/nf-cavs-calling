from scipy.special import logit, expit, logsumexp
from scipy import stats as st
import numpy as np
from typing import Tuple
from scipy import interpolate



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


def aggregate_effect_size(es, weights):
    """
    Aggregate effect sizes using weights

    Args:
    - es: effect sizes shape (n_samples, ...)
    - weights: weights shape (n_samples,)
    """
    return np.average(es, weights=weights, axis=0)

def logit_es(p, d=1/128):
    return np.log2(p + d) - np.log2(1 - p + d)


def inv_logit_es(es, d=1/128):
    s = 2 ** es
    return (1 + (s - 1) * (d + 1)) / (s + 1)

## P-values

def calc_bimodal_pvalues(dist1: st.rv_discrete, dist2: st.rv_discrete, x: np.ndarray, w: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    log_p_right = logsumexp(
        [dist1.logsf(x - 1), dist2.logsf(x - 1)], 
        b=[w, 1 - w]
    )
    log_p_left = logsumexp(
        [dist1.logcdf(x), dist2.logcdf(x)], 
        b=[w, 1 - w]
    )
    log_p_both = log_pval_both(log_p_right, log_p_left)
    return log_p_right, log_p_left, log_p_both


def log_pval_both(log_p_right, log_p_left):
    return np.min(np.stack([log_p_right, log_p_left]), axis=0) - np.log(2)


def stouffer_combine_log_pvals(log_pvals, weights):
    """
    A vectorized version of Stouffer's method for combining p-values
    """
    if weights is None:
        weights = np.ones_like(log_pvals)
    pvals = np.exp(log_pvals)
    # return st.combine_pvalues(pvals, weights=weights, method='stouffer')[1]
    return st.norm.logsf(weights.dot(st.norm.isf(pvals)) / np.linalg.norm(weights))


# implementation of Storey method for FDR estimation
def qvalue(pvals, bootstrap=False):
    m, pvals = len(pvals), np.asarray(pvals)
    ind = np.argsort(pvals)
    rev_ind = np.argsort(ind)
    pvals = pvals[ind]
    # Estimate proportion of features that are truly null.
    kappa = np.arange(0.05, 0.96, 0.01)
    pik = np.array([sum(pvals > k) / (m * (1-k)) for k in kappa])

    if bootstrap:
        minpi0 = np.quantile(pik, 0.1)
        W = np.array([(pvals >= l).sum() for l in kappa])
        mse = (W / (np.square(m^2) * np.square(1 - kappa))) * (1 - (W / m)) + np.square((pik - minpi0))
        
        if np.any(np.isnan(mse)) or np.any(np.isinf(mse)):
            # Case 1: mse contains NaN or Inf
            mask = np.isnan(mse)  # This will return a boolean mask where True indicates NaN positions

        else:
            # Case 2: mse contains only finite values
            mask = mse == mse.min()
        pi0 = pik[mask][0]
    else:
        cs = interpolate.UnivariateSpline(kappa, pik, k=3, s=None, ext=0)
        pi0 = float(cs(1.))
    
    pi0 = min(pi0, 1)
    # Compute the q-values.
    qvals = np.zeros(len(pvals))
    qvals[-1] = pi0 * pvals[-1]
    for i in np.arange(m - 2, -1, -1):
        qvals[i] = min(pi0 * m * pvals[i]/float(i+1), qvals[i+1])
    qvals = qvals[rev_ind]
    return qvals