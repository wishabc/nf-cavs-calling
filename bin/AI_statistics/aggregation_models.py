import numpy as np
from base_models import cached_method, EffectModel, BimodalBaseModel, BimodalSamplingModel, BimodalScoringModel, SamplingModel, ScoringModel, Pvalues
from collections.abc import Sequence
from vectorized_estimators import stouffer_combine_log_pvals, aggregate_effect_size, log_pval_both, generate_cartesian_product


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


class AggregatedBimodalModel(EffectModel):
    __child_model__ = BimodalBaseModel

    def __init__(self, models: Sequence[BimodalBaseModel], indivs=None, weights=None,
                 size_constraint=1e7):
        assert all(isinstance(model, BimodalBaseModel) for model in models)
        self.models = [self.__child_model__.from_model(model) for model in models]
        len_models = len(self.models)
        self.indivs = _validate_list_argument(len_models, indivs)
        self.weights = _validate_list_argument(len_models, weights)
        if weights is None:
            self.weights *= np.sqrt(np.array([model.n for model in self.models]))

        self.size_constraint = size_constraint

        self._random_state_mod = 2 ** 32
        self.e = aggregate_effect_size([model.e for model in self.models], weights=self.weights)

    @classmethod
    def from_model(cls, other: 'AggregatedBimodalModel', effect=None):
        new_models = [model.from_model(model, effect=effect) for model in other.models]
        return cls(new_models, other.indivs, other.weights)

    def compatible_with(self, other: 'AggregatedBimodalModel'):
        return np.all(x.compatible_with(y) for x, y in zip(self.models, other.models))

    def check_sizes(self):
        cartesian_size = np.prod(np.array([model.n for model in self.models]))
        if cartesian_size > self.size_constraint:
            raise ValueError(f"(!) Cartesian product of {cartesian_size} elements is too large. Increase size_constraint to proceed.")

    @property
    def all_observations(self):
        self.check_sizes()
        individual_observations = [model.all_observations for model in self.models]
        return generate_cartesian_product(individual_observations)
    
    def get_log_pmf_for_mode(self, bad_phasing_mode):
        self.check_sizes()
        individual_log_pmfs = [model.get_log_pmf_for_mode(bad_phasing_mode) for model in self.models]
        return generate_cartesian_product(individual_log_pmfs)


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

    def calc_log_pvalues(self, samples) -> Pvalues:
        self.check_samples(samples)
        log_p_right, log_p_left, _ = map(
            np.stack,
            zip(*[model.calc_pvalues(sample) for sample, model in zip(samples, self.models)])
        )
        agg_log_p_right = stouffer_combine_log_pvals(log_p_right, self.weights) 
        agg_log_p_left = stouffer_combine_log_pvals(log_p_left, self.weights)
        agg_log_p = log_pval_both(agg_log_p_left, agg_log_p_right)
        return Pvalues([agg_log_p_right, agg_log_p_left, agg_log_p])

    def get_effect_size_frac(self, samples):
        self.check_samples(samples)
        return aggregate_effect_size([model.effect_size_estimate(sample, return_frac=True) for sample, model in zip(samples, self.models)], weights=self.weights)
    
    def check_samples(self, samples):
        samples = np.asarray(samples)
        assert samples.ndim == 2, f'(!) samples should be 2D array, got: {samples.ndim}'
        assert samples.shape[0] == len(self.models), f'(!) number of models does not match number of samples'
