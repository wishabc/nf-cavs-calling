import numpy as np
from base_models import cached_method, BimodalEffectModel, BimodalSamplingModel, BimodalScoringModel, SamplingModel, ScoringModel
from collections.abc import Sequence
from vectorized_estimators import stouffer_combine_log_pvals, aggregate_effect_size, log_pval_both


def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False
    

def _validate_list_argument(len_models, arg):
    if arg is None:
        arg = np.ones(len_models)
    elif is_iterable(arg):
        arg = np.array(arg)
        assert len(arg) == len_models
    else:
        arg = np.ones(len_models) * arg
    return arg


class AggregatedBimodalModel:
    __child_model__ = BimodalEffectModel

    def __init__(self, models: Sequence[BimodalEffectModel], indivs=None, weights=None):
        assert all(isinstance(model, BimodalEffectModel) for model in models)
        self.models = [self.__child_model__.from_model(model) for model in models]
        len_models = len(self.models)
        self.indivs = _validate_list_argument(len_models, indivs)
        self.weights = _validate_list_argument(len_models, weights)
        if weights is None:
            self.weights *= np.sqrt(np.array([model.n for model in self.models]))

        self._random_state_mod = 2 ** 32
        self.e = aggregate_effect_size([model.e for model in self.models], weights=self.weights)

    @cached_method
    def get_effect_model(self, effect):
        """
        Get a model with a different effect size
        """
        if effect == self.e:
            return self
        return self.__class__.from_model(self, effect=effect)

    @classmethod
    def from_model(cls, other: 'AggregatedBimodalModel', effect=None):
        new_models = [model.from_model(model, effect=effect) for model in other.models]
        return cls(new_models, other.indivs, other.weights)

    def compatible_with(self, other: 'AggregatedBimodalModel'):
        return np.all(x.compatible_with(y) for x, y in zip(self.models, other.models))

    @property
    def all_observations(self):
        raise NotImplementedError


class AggregatedBimodalSamplingModel(SamplingModel, AggregatedBimodalModel):
    __child_model__ = BimodalSamplingModel

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


class AggregatedBimodalScoringModel(ScoringModel, AggregatedBimodalModel):
    __child_model__ = BimodalScoringModel

    def calc_log_pvalues(self, samples):
        self.check_samples(samples)
        log_p_right, log_p_left, _ = map(
            np.stack,
            zip(*[model.calc_pvalues(sample) for sample, model in zip(samples, self.models)])
        )
        agg_log_p_right = stouffer_combine_log_pvals(log_p_right, self.weights) 
        agg_log_p_left = stouffer_combine_log_pvals(log_p_left, self.weights)
        agg_log_p = log_pval_both(agg_log_p_left, agg_log_p_right)
        return agg_log_p_right, agg_log_p_left, agg_log_p

    def effect_size_estimate(self, samples):
        self.check_samples(samples)
        return aggregate_effect_size([model.effect_size_estimate(sample) for sample, model in zip(samples, self.models)], weights=self.weights)
    
    def check_samples(self, samples):
        samples = np.asarray(samples)
        assert samples.ndim == 2, f'(!) samples should be 2D array, got: {samples.ndim}'
        assert samples.shape[0] == len(self.models), f'(!) number of models does not match number of samples'