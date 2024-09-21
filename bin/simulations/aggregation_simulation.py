import numpy as np
import scipy.stats as st
from simulations import cached_method, BinomialSamplingModel, BinomialScoringModel, BinomialModel

"""
Functions and classes to estimate detection power of aggregation.
Not used as part of the pipeline
"""

def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False
    

def aggregate_pvals(pvals, weights):
    """
    A vectorized version of Stouffer's method for combining p-values
    """
    if weights is None:
        weights = np.ones_like(pvals)
    # return st.combine_pvalues(pvals, weights=weights, method='stouffer')[1]
    return st.norm.logsf(weights.dot(st.norm.isf(pvals)) / np.linalg.norm(weights))


class AggregatedBinomialModel:
    __child_model__ = BinomialModel
    def __init__(self, ns, effect=0, Bs=None, indivs=None):
        self.e = effect
        self.ns = np.array(ns)
        self.indivs = self._validate_list_argument(indivs)
        self.Bs = self._validate_list_argument(Bs)
        self.weights = self.ns  # np.sqrt(self.ns)
        self._random_state_mod = 2 ** 32
        self.models = [self.__child_model__(n, effect, B) for n, B in zip(self.ns, self.Bs)]

    @cached_method
    def get_effect_model(self, effect=0):
        if effect == self.e:
            return self
        return self.__class__(self.ns, effect, Bs=self.Bs)
    
    def _validate_list_argument(self, arg):
        if arg is None:
            arg = np.ones(len(self.ns))
        elif is_iterable(arg):
            arg = np.array(arg)
            assert len(arg) == len(self.ns)
        else:
            arg = np.ones(len(self.ns)) * arg
        return arg

    @classmethod
    def from_model(cls, other: 'AggregatedBinomialModel'):
        return cls(other.ns, other.e, other.Bs, other.indivs)

class AggregatedBinomialSamplingModel(AggregatedBinomialModel):
    __child_model__ = BinomialSamplingModel

    @cached_method
    def get_samples(self, n_itter, random_state=0, bad_phasing_mode=None):
        """
        Get n_itter samples from each model
        returns:
        - samples: np.array of shape (n_models, n_itter), dtype=int
        - phase: np.array of shape (n_models, n_itter), dtype=int
        - bad_phasing_mode: None, 1 or 2, if not None, all models are forced to use the specified BAD phasing mode
        """
        assert bad_phasing_mode in [1, 2, None]
        samples = np.empty((len(self.models), n_itter), dtype=int)
        phases = np.empty((len(self.models), n_itter), dtype=int)
        for i, (model, indiv) in enumerate(zip(self.models, self.indivs)):
            sample, phase = model.get_samples(
                    n_itter,
                    random_state=(random_state + n_itter * 10 * i) % self._random_state_mod,
                    bad_phasing_mode=bad_phasing_mode,
                    phase_random_state=(random_state + hash(indiv)) % self._random_state_mod,
                )
            samples[i, :] = sample
            phases[i, :] = phase

        return samples, phases


class AggregatedBinomialScoringModel(AggregatedBinomialModel):
    __child_model__ = BinomialScoringModel

    def aggregate_pvals(self, pvals, weights):
        """
        A vectorized version of Stouffer's method for combining p-values
        """
        if weights is None:
            weights = np.ones_like(pvals)
        # return st.combine_pvalues(pvals, weights=self.weights, method='stouffer')
        return st.norm.logsf(self.weights.dot(st.norm.isf(pvals)) / np.linalg.norm(self.weights))

    @cached_method
    def aggregated_log_p_values(self, samples):
        p_right, p_left, p_both = map(
            np.stack,
            zip(*[model.calc_pvalues(sample) for sample, model in zip(samples, self.models)])
        )
        agg_log_p_right = aggregate_pvals(p_right, self.weights) 
        agg_log_p_left = aggregate_pvals(p_left, self.weights)
        #         agg_log_p = aggregate_pvals(p_both, self.weights)
        agg_log_p = np.log(2) + np.min([agg_log_p_left, agg_log_p_right], axis=0)
        side = np.where(agg_log_p_right < agg_log_p_left, 1, -1)
        return agg_log_p, side


class AggregatedSamplingPowerEstimator:
    def __init__(self, null_model: AggregatedBinomialSamplingModel, scoring_model: AggregatedBinomialScoringModel, n_itter=10000, random_state=None, bad_phasing_mode=None):
        assert null_model.ns == scoring_model.ns, 'Models have different ns'
        if random_state is None:
            random_state = np.random.randint(0, 100000)
        self.random_state = random_state
        self.n_itter = n_itter
        self.null_model = null_model
        self.scoring_model = scoring_model
        assert bad_phasing_mode in [1, 2, None]
        self.bad_phasing_mode = bad_phasing_mode

    @cached_method
    def sensitivity(self, effect, signif_tr, correct_indices=False):
        if effect == self.null_model.e:
            print('(!) Comparing same effect size models')
        effect_model = self.null_model.get_effect_model(effect)
        samples, _ = effect_model.get_samples(self.n_itter, random_state=self.random_state, bad_phasing_mode=self.bad_phasing_mode)
        log_pvals, side = self.scoring_model.aggregated_log_p_values(samples)
        if correct_indices:
            correct_indices = side == (1 if effect - self.null_model.e > 0 else -1)
        else:
            correct_indices = np.ones_like(log_pvals, dtype=bool)
        return np.sum((log_pvals <= np.log(signif_tr)) & correct_indices) / self.n_itter

    @cached_method
    def specificity(self, signif_tr):
        samples, _ = self.null_model.get_samples(self.n_itter, random_state=self.random_state, bad_phasing_mode=self.bad_phasing_mode)
        log_pvals, _ = self.scoring_model.aggregated_log_p_values(samples)
        return 1 - np.sum(log_pvals <= np.log(signif_tr)) / self.n_itter

    @cached_method
    def correct_side_sensitivity(self, effect, signif_tr):
        return self.sensitivity(effect=effect, signif_tr=signif_tr, correct_indices=True)
