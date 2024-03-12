import numpy as np
import scipy.stats as st
from scipy.special import expit, logsumexp
from functools import wraps


def is_iterable(obj):
    try:
        iter(obj)
        return True
    except TypeError:
        return False


def cached_property(func=None, *, init_method=None):
    """
    A decorator for cached properties. The property is calculated once and then stored in the instance.
    If the property is accessed again, the stored value is returned.
    If an init_method is provided, it is called to initialize all properties.
    """
    if func is None:
        # If the decorator is used without arguments,
        # it acts as a factory for the actual decorator,
        # passing the init_method to the actual decorator
        return lambda f: cached_property(f, init_method=init_method)

    cache_attr_name = f"_{func.__name__}"

    @property
    def wrapper(self):
        if not hasattr(self, cache_attr_name):
            # If an init_method is provided, call it to initialize all properties
            if init_method is not None:
                getattr(self, init_method)()
            else:
                setattr(self, cache_attr_name, func(self))
        return getattr(self, cache_attr_name)

    @wrapper.setter
    def wrapper(self, value):
        setattr(self, cache_attr_name, value)

    return wrapper


def hash_array(array):
    return hash(array.tobytes())


def cached_method(func):
    """
    A decorator for cached methods.
    The result of the method is calculated once for each set of arguments and then stored in the instance.
    If the method is called again with the same arguments, the stored value is returned.
    """
    cache_attr_name = f"_{func.__name__}_cache"

    @wraps(func)
    def wrapper(self, *args, **kwargs):
        if not hasattr(self, cache_attr_name):
            setattr(self, cache_attr_name, {})
        cache = getattr(self, cache_attr_name)

        hashed_args = tuple(
            hash_array(arg) if isinstance(arg, np.ndarray) else arg
            for arg in args
        )
        key = (hashed_args, tuple(sorted(kwargs.items())))
        if key not in cache:
            cache[key] = func(self, *args, **kwargs)
        return cache[key]

    return wrapper


def aggregate_pvals(pvals, weights):
    """
    A vectorized version of Stouffer's method for combining p-values
    """
    if weights is None:
        weights = np.ones_like(pvals)
    # return st.combine_pvalues(pvals, weights=weights, method='stouffer')[1]
    return st.norm.logsf(weights.dot(st.norm.isf(pvals)) / np.linalg.norm(weights))


class BinomialAllelicImbalanceModel:
    """
    A model to generate allelic imbalance data
    """

    def __init__(self, n, effect=0, B=1):
        """
        Args:
        - n: number of reads
        - effect: effect size
        - B: background allelic dosage
        """
        self.n = n
        self.e = effect
        self.p = expit(self.e * np.log(2) + np.log(B))
        self.p2 = expit(self.e * np.log(2) - np.log(B))
        self.B = B

        self.dist = st.binom(n, self.p)
        self.dist2 = st.binom(n, self.p2)

        self._x = np.arange(self.n + 1)

    @cached_method
    def get_samples(self, n_itter, random_state, phased=True):
        """
        Get n_itter samples from the model
        returns:
        - samples: np.array of shape (n_itter,), dtype=int, number of reads
        - phase: np.array of shape (n_itter,), dtype=bool, copy number phase (according to the major allele)
        """
        if phased:
            return self.dist.rvs(size=n_itter, random_state=random_state), np.ones(n_itter, dtype=int)
        else:
            n_phased = st.binom(n_itter, 0.5).rvs(random_state=random_state)
            samples = np.concatenate([self.dist.rvs(size=n_phased, random_state=random_state),
                                      self.dist2.rvs(size=n_itter - n_phased, random_state=random_state)])
            phase = np.concatenate([np.ones(n_phased, dtype=int), np.zeros(n_itter - n_phased, dtype=int)])
            order = np.argsort(np.random.rand(n_itter))
            return samples[order], phase[order]

    @cached_method
    def get_effect_model(self, effect=0):
        if effect == self.e:
            return self
        return BinomialAllelicImbalanceModel(self.n, effect, B=self.B)

    @property
    def true_log_pmf(self):
        return self.dist.logpmf(self._x)

    @property
    def bimodal_pmf(self):
        return 0.5 * (self.dist.pmf(self._x) + self.dist2.pmf(self._x))

    def estimate_w(self, x):
        if self.e == 0:
            b = np.log(self.B)
            delta = b * (self.n - 2 * x)
            return expit(-delta)
        return expit(self.dist.logpmf(x) - self.dist2.logpmf(x))

    def sum_probability(self, indicators):
        return np.exp(logsumexp(self.true_log_pmf[indicators])) if indicators.sum() != 0 else 0

    @property
    def x(self):
        return self._x


class BinomialAllelicImbalanceScoringModel(BinomialAllelicImbalanceModel):
    """
    A model to score allelic imbalance data
    """

    def __init__(self, other, **scoring_params):
        super().__init__(
            other.n,
            scoring_params.get('e', other.e),
            scoring_params.get('B', other.B)
        )

    def calc_pvalues(self, x):
        w = self.estimate_w(x)
        p_right = w * self.dist.sf(x - 1) + (1 - w) * self.dist2.sf(x - 1)
        p_left = w * self.dist.cdf(x) + (1 - w) * self.dist2.cdf(x)
        p_both = 2 * np.minimum(p_right, p_left)
        # p_both_alt = 2 * (
        #     w * np.minimum(self.dist.sf(x - 1), self.dist.cdf(x)) +
        #     (1 - w) * np.minimum(self.dist2.sf(x - 1), self.dist2.cdf(x))
        # )
        return p_right, p_left, p_both

    def _calc_all_pvalues(self):
        self._p_right, self._p_left, self._p_both = self.calc_pvalues(self.x)

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

    def predicted_imbalanced_idx(self, signif_tr, side='both'):
        return self.p_value(side=side) <= signif_tr


class ExactPowerEstimator:
    def __init__(self, null_model, **scoring_params):
        self.null_model = null_model
        self.scoring_model = BinomialAllelicImbalanceScoringModel(null_model, **scoring_params)

    @cached_method
    def get_signif_indices(self, signif_tr, side='both'):
        return self.scoring_model.predicted_imbalanced_idx(signif_tr, side=side)

    @cached_method
    def sensitivity(self, effect, signif_tr, side='both'):
        model = self.null_model.get_effect_model(effect)
        return model.sum_probability(self.get_signif_indices(signif_tr, side=side))

    @cached_method
    def specificity(self, signif_tr, side='both'):
        return 1 - self.null_model.sum_probability(self.get_signif_indices(signif_tr, side=side))

    @cached_method
    def correct_side_sensitivity(self, effect, signif_tr, side='both'):
        model = self.null_model.get_effect_model(effect)
        imbalanced_indices = self.get_signif_indices(signif_tr, side=side)
        if effect == self.null_model.e:
            print('(!) Comparing same effect size models')
            return model.sum_probability(imbalanced_indices)
        correct_side_indices = (self.scoring_model.p_right < self.scoring_model.p_left) == (
                    (effect - self.null_model.e) > 0)
        return model.sum_probability(imbalanced_indices & correct_side_indices)


class AggregatedBinomialModel:
    def __init__(self, ns, effect=0, Bs=None):
        self.ns = np.array(ns)
        self.e = effect
        if Bs is None:
            Bs = np.ones(len(ns))
        elif is_iterable(Bs):
            Bs = np.array(Bs)
            assert len(Bs) == len(ns)
        else:
            Bs = np.ones(len(ns)) * Bs
        self.Bs = Bs
        self.models = [BinomialAllelicImbalanceModel(n, effect, B) for n, B in zip(ns, Bs)]
        self.weights = self.ns  # np.sqrt(self.ns)

    @cached_method
    def get_effect_model(self, effect=0):
        if effect == self.e:
            return self
        return AggregatedBinomialModel(self.ns, effect, Bs=self.Bs)

    @cached_method
    def get_samples(self, n_itter, random_state=0, phased=True):
        """
        Get n_itter samples from each model
        returns:
        - samples: np.array of shape (n_models, n_itter), dtype=int
        - phase: np.array of shape (n_models, n_itter), dtype=bool
        """
        samples, phase = zip(
            *[model.get_samples(n_itter, random_state=random_state + n_itter * 10 * i, phased=phased) for i, model in
              enumerate(self.models)])
        return np.stack(samples), np.stack(phase)


class AggregatedBinomialScoringModel(AggregatedBinomialModel):
    def __init__(self, other, **scoring_params):
        super().__init__(
            other.ns,
            scoring_params.get('e', other.e),
            scoring_params.get('Bs', other.Bs)
        )
        self.models = [BinomialAllelicImbalanceScoringModel(model) for model in self.models]

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
        return agg_log_p_right, agg_log_p_left, agg_log_p, side


class AggregatedSamplingPowerEstimator:
    def __init__(self, null_model, n_itter=10000, random_state=None, phased=True, **scoring_params):
        if random_state is None:
            random_state = np.random.randint(0, 100000)
        self.random_state = random_state
        self.n_itter = n_itter
        self.null_model = null_model
        self.scoring_model = AggregatedBinomialScoringModel(null_model, **scoring_params)
        self.phased = phased

    @cached_method
    def sensitivity(self, effect, signif_tr):
        model = self.null_model.get_effect_model(effect)
        samples, _ = model.get_samples(self.n_itter, random_state=self.random_state, phased=self.phased)
        _, _, log_pvals, _ = self.scoring_model.aggregated_log_p_values(samples)
        return np.sum(log_pvals <= np.log(signif_tr)) / self.n_itter

    @cached_method
    def specificity(self, signif_tr):
        samples, _ = self.null_model.get_samples(self.n_itter, random_state=self.random_state, phased=self.phased)
        _, _, log_pvals, _ = self.scoring_model.aggregated_log_p_values(samples)
        return 1 - np.sum(log_pvals <= np.log(signif_tr)) / self.n_itter

    @cached_method
    def correct_side_sensitivity(self, effect, signif_tr):
        model = self.null_model.get_effect_model(effect)
        samples, _ = model.get_samples(self.n_itter, random_state=self.random_state, phased=self.phased)
        log_p_right, log_p_left, log_pvals, side = self.scoring_model.aggregated_log_p_values(samples)
        if effect == self.null_model.e:
            print('(!) Comparing same effect size models')
        else:
            correct_side_indices = side == int((effect - self.null_model.e) > 0) * 2 - 1
            log_pvals[~correct_side_indices] = 0
        return np.sum(log_pvals <= np.log(signif_tr)) / self.n_itter
