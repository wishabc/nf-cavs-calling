from scipy.special import logsumexp, expit, logit
from base_models import cached_method, BimodalScoringModel, BimodalEffectModel, cached_property
from aggregation_model import AggregatedBimodalSamplingModel, AggregatedBimodalScoringModel
import numpy as np


class CachedScoringModel:
    """
    A model with pre-calculated p-values for all possible x values
    """

    def __init__(self, model: BimodalScoringModel):
        self.model = model
        if not isinstance(model, BimodalScoringModel):
            raise TypeError(f"{model.__class__.__name__} must inherit from BimodalScoringModel")
        self.all_observations = model.all_observations
    
    def _calc_all_pvalues(self):
        self._p_right, self._p_left, self._p_both = self.model.calc_pvalues(self.all_observations)
    
    @cached_property(init_method='_calc_all_pvalues')
    def p_right(self):
        return self._p_right

    @cached_property(init_method='_calc_all_pvalues')
    def p_left(self):
        return self._p_left

    @cached_property(init_method='_calc_all_pvalues')
    def p_both(self):
        return self._p_both
    
    def p_value(self, side='both'):
        assert side in ['right', 'left', 'both']
        if side == 'right':
            return self.p_right
        elif side == 'left':
            return self.p_left
        elif side == 'both':
            return self.p_both
        
    @cached_property
    def effect_size_estimates(self):
        return self.model.effect_size_estimate(self.all_observations)
    
    def get_signif_indices(self, signif_tr, side='both'):
        return self.p_value(side=side) <= signif_tr


class ExactPowerEstimator:
    def __init__(self, null_model: BimodalEffectModel, scoring_model: CachedScoringModel, bad_phasing_mode=None):
        self.null_model = null_model
        self.scoring_model = scoring_model
        assert bad_phasing_mode in [None, 1, 2]
        assert scoring_model.model.compatible_with(null_model)
        self.bad_phasing_mode = bad_phasing_mode

    @cached_method
    def sensitivity(self, effect, signif_tr, side='both', correct_indices=None):
        """
        Calculate the sensitivity of the scoring model to a given effect size

        Args:
            - effect: effect size
            - signif_tr: p-value significance threshold
            - side: 'right', 'left' or 'both'
            - correct_indices: indices of self.all_x with the correct side of the effect. If None, all indices are used
        """
        if effect == self.null_model.e:
            print('(!) Effect and null models are the same (!)')
        effect_model = self.null_model.get_effect_model(effect)
        signif_indices = self.scoring_model.get_signif_indices(signif_tr, side=side)
        if correct_indices is None:
            correct_indices = np.ones_like(signif_indices, dtype=bool)
        return self.sum_probability_for_mode(effect_model, signif_indices & correct_indices)

    @cached_method
    def specificity(self, signif_tr, side='both'):
        signif_indices = self.scoring_model.get_signif_indices(signif_tr, side=side)
        return 1 - self.sum_probability_for_mode(self.null_model, signif_indices)

    def correct_side_sensitivity(self, effect, signif_tr, side='both'):
        correct_side_indices = (self.scoring_model.p_right < self.scoring_model.p_left) == (
                    (effect - self.null_model.e) > 0)
        return self.sensitivity(effect, signif_tr, side=side, correct_indices=correct_side_indices)

    def sum_probability_for_mode(self, model: BimodalEffectModel, indicators: np.ndarray):
        """
        Calculate the sum of probabilities for the given mode

        For true distribution, knowing the correct BAD phasing mode.
        By default positive effect corresponds to preference towards copied allele (mode=1)  
        """
        log_pmf_for_mode = model.get_log_pmf_for_mode(self.bad_phasing_mode)
        return np.exp(logsumexp(log_pmf_for_mode[indicators])) if indicators.sum() != 0 else 0
    
    @cached_method
    def get_effect_log_pmf(self, effect):
        effect_model = self.null_model.get_effect_model(effect)
        return effect_model.get_log_pmf_for_mode(self.bad_phasing_mode)
    
    @cached_method
    def fraction_effect_size_stats(self, effect):
        log_pmf = self.get_effect_log_pmf(effect)
        p_estimates = expit(self.scoring_model.effect_size_estimates * np.log(2))
        p_expectation = np.exp(logsumexp(
                            log_pmf,
                            b=p_estimates
                        ))
        p_variance = np.exp(logsumexp(
                            log_pmf,
                            b=(p_estimates - p_expectation) ** 2
                        ))
        return p_expectation, p_variance


class AggregatedSamplingPowerEstimator:
    def __init__(self, null_model: AggregatedBimodalSamplingModel, scoring_model: AggregatedBimodalScoringModel, n_itter=10000, random_state=None, bad_phasing_mode=None):
        scoring_model.compatible_with(null_model)
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
