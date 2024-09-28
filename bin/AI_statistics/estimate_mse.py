from scipy.special import logit, expit
from scipy import stats as st
import pandas as pd
import sys
import numpy as np
from tqdm import tqdm

tqdm.pandas()

def effect_size_estimate(x, n, B):
    x, n, B = map(np.asarray, [x, n, B])

    b = np.log(B)
    delta = b * (n - 2 * x)
    p = x / n
    return np.where(
        p == 0.0, 0.0,
        np.where(
            p == 1.0, 1.0,
            expit(expit(-delta) * (logit(p) - b) + expit(delta) * (logit(p) + b))
        )
    )

def expectation_vectorized(func, n, B, p):
    b = np.log(B)
    p1 = expit(logit(p) + b)
    x = np.arange(0, n + 1)
    
    pmf_values = st.binom.pmf(x[:, None], n, p1[None, :])
    func_values = func(x[:, None], n, B )
    
    expectations = np.sum(pmf_values * func_values, axis=0)
    return expectations

def es_var_vectorized(n, B, p):
    exp = expectation_vectorized(effect_size_estimate, n=n, B=B, p=p)
    def x2(x, n, B):
        return effect_size_estimate(x, n, B ) ** 2
    return expectation_vectorized(x2, n=n, B=B, p=p) - exp ** 2

def calc_variance(n, B, n_points=101):
    ess = np.linspace(0, 1, n_points)
    yvals = expectation_vectorized(effect_size_estimate, n, B, ess)
    vars = es_var_vectorized(n, B, ess)
    return (yvals - ess) ** 2 + vars



# if __name__ == '__main__':
#     df = pd.read_table(sys.argv[1])
#     initial_df_len = len(df.index)
    
#     data = df[['BAD', 'coverage']].drop_duplicates()
#     data['inverse_mse'] = data.progress_apply(lambda row: 1 / calc_mse(row['coverage'], row['BAD']).mean(), axis=1)
#     df = df.merge(data, on=['BAD', 'coverage'])
#     assert initial_df_len == len(df.index)

#     df.to_csv(sys.argv[1], index=False, sep='\t')