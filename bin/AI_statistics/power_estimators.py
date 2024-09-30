from scipy.special import logsumexp
from base_models import cached_method, ScoringModel, cached_property, EffectModel, SamplingModel
import numpy as np


class CachedScoringModel:
    """
    A model with pre-calculated p-values for all possible x values
    """

    def __init__(self, model: ScoringModel):
        self.model = model
        if not isinstance(model, ScoringModel):
            raise TypeError(f"{model.__class__.__name__} must inherit from ScoringModel")
        self.all_observations = model.all_observations
    
    def _calc_all_pvalues(self):
        self._log_p_right, self._log_p_left, self._log_p_both = self.model.calc_log_pvalues(self.all_observations)
    
    @cached_property(init_method='_calc_all_pvalues')
    def log_p_right(self):
        return self._log_p_right

    @cached_property(init_method='_calc_all_pvalues')
    def log_p_left(self):
        return self._log_p_left

    @cached_property(init_method='_calc_all_pvalues')
    def log_p_both(self):
        return self._log_p_both
    
    def log_p_value(self, side='both'):
        assert side in ['right', 'left', 'both']
        if side == 'right':
            return self.log_p_right
        elif side == 'left':
            return self.log_p_left
        elif side == 'both':
            return self.log_p_both
        
    @cached_property
    def effect_size_estimates(self):
        return self.model.calc_effect_size(self.all_observations)
    
    @cached_property
    def fraction_effect_size_estimates(self):
        return self.model.calc_effect_size(self.all_observations, return_frac=True)
    
    def get_signif_indices(self, signif_tr, side='both'):
        return self.log_p_value(side=side) <= np.log(signif_tr)


class ExactPowerEstimator:
    def __init__(self, null_model: EffectModel, scoring_model: CachedScoringModel, bad_phasing_mode=None):
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
    def false_positive_rate(self, signif_tr, side='both'):
        signif_indices = self.scoring_model.get_signif_indices(signif_tr, side=side)
        return self.sum_probability_for_mode(self.null_model, signif_indices)

    def correct_side_sensitivity(self, effect, signif_tr, side='both'):
        correct_side_indices = (self.scoring_model.log_p_right < self.scoring_model.log_p_left) == ((effect - self.null_model.e) > 0)
        return self.sensitivity(effect, signif_tr, side=side, correct_indices=correct_side_indices)

    def sum_probability_for_mode(self, model: EffectModel, indicators: np.ndarray):
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
        p_estimates = self.scoring_model.fraction_effect_size_estimates
        p_expectation = np.exp(
            logsumexp(
                log_pmf,
                b=p_estimates
            )
        )
        p_variance = np.exp(
            logsumexp(
                log_pmf,
                b=(p_estimates - p_expectation) ** 2
            )
        )
        return p_expectation, p_variance


class SamplingPowerEstimator:
    def __init__(self, null_model: SamplingModel, scoring_model: ScoringModel, n_itter=10000, random_state=None, bad_phasing_mode=None):
        assert scoring_model.compatible_with(null_model), 'Scoring model cannot score the null model, wrong Ns'
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
        log_pval_right, log_pval_left, log_pval_both = self.scoring_model.calc_log_pvalues(samples)
        side = np.where(log_pval_right < log_pval_left, 1, -1)
        if correct_indices:
            correct_indices = side == (1 if effect - self.null_model.e > 0 else -1)
        else:
            correct_indices = np.ones_like(log_pval_both, dtype=bool)
        return np.sum((log_pval_both <= np.log(signif_tr)) & correct_indices) / self.n_itter

    @cached_method
    def false_positive_rate(self, signif_tr):
        samples, _ = self.null_model.get_samples(self.n_itter, random_state=self.random_state, bad_phasing_mode=self.bad_phasing_mode)
        log_pvals = self.scoring_model.calc_log_pvalues(samples).both
        return np.sum(log_pvals <= np.log(signif_tr)) / self.n_itter

    @cached_method
    def correct_side_sensitivity(self, effect, signif_tr):
        return self.sensitivity(effect=effect, signif_tr=signif_tr, correct_indices=True)
    
    def fraction_effect_size_stats(self, effect):
        effect_model = self.null_model.get_effect_model(effect)
        samples, _ = effect_model.get_samples(self.n_itter, random_state=self.random_state, bad_phasing_mode=self.bad_phasing_mode)
        p_estimates = self.scoring_model.calc_effect_size(samples, return_frac=True)
        p_expectation = np.mean(p_estimates)
        p_variance = np.var(p_estimates)
        return p_expectation, p_variance
