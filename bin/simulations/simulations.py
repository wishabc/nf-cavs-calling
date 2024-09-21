import numpy as np
import scipy.stats as st
from scipy.special import expit, logsumexp
from functools import wraps


"""
Functions and classes to simulate allelic imbalance data and estimate detection power.
Not used as part of the pipeline
"""


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


class BinomialModel:
    """
    Base class to store parameters of a binomial model
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
        self.p1 = expit(self.e * np.log(2) + np.log(B))
        self.p2 = expit(self.e * np.log(2) - np.log(B))
        self.B = B

        self.dist1 = st.binom(n, self.p1)
        self.dist2 = st.binom(n, self.p2)

        self.all_x = np.arange(self.n + 1)
    
    @classmethod
    def from_model(cls, other: 'BinomialModel'):
        return cls(other.n, other.e, other.B)
    
    @cached_method
    def get_effect_model(self, effect=0):
        """
        Get a model with a different effect size
        """
        if effect == self.e:
            return self
        return self.__class__(self.n, effect, B=self.B)


class BinomialSamplingModel(BinomialModel):
    """
    A model to simulate allelic imbalance data
    """

    @cached_method
    def get_samples(self, n_itter, random_state, bad_phasing_mode=None, phase_random_state=None):
        """
        Get n_itter samples from the model
        
        Args:
        - n_itter: number of samples
        - random_state: random state
        - bad_phasing_mode: None, 1 or 2, if not None, the model is forced to use the specified BAD phasing mode

        returns:
        - samples: np.array of shape (n_itter,), dtype=int, number of reads
        - phase: np.array of shape (n_itter,), dtype=int, copy number phase (according to the major allele)
        """
        assert bad_phasing_mode in [1, 2, None]
        if bad_phasing_mode is not None:
            dist = self.dist1 if bad_phasing_mode == 1 else self.dist2
            return dist.rvs(size=n_itter, random_state=random_state), bad_phasing_mode * np.ones(n_itter, dtype=int)
        else:
            if phase_random_state is None:
                phase_random_state = random_state
            n_phased = st.binom(n_itter, 0.5).rvs(random_state=phase_random_state)
            samples = np.concatenate(
                [
                    self.dist1.rvs(
                        size=n_phased,
                        random_state=random_state
                    ),
                    self.dist2.rvs(
                        size=n_itter - n_phased,
                        random_state=random_state
                    )
                ]
            )
            phase = np.ones(n_itter, dtype=int)
            phase[n_phased:] = 2
            rng = np.random.RandomState(phase_random_state)
            order = np.argsort(rng.rand(n_itter)) 
            return samples[order], phase[order]



class BinomialScoringModel(BinomialModel):
    """
    A model to score allelic imbalance data
    """
    def calc_pvalues(self, x):
        w = self.estimate_w(x)
        p_right = w * self.dist1.sf(x - 1) + (1 - w) * self.dist2.sf(x - 1)
        p_left = w * self.dist1.cdf(x) + (1 - w) * self.dist2.cdf(x)
        p_both = 2 * np.minimum(p_right, p_left)
        # p_both_alt = 2 * (
        #     w * np.minimum(self.dist.sf(x - 1), self.dist.cdf(x)) +
        #     (1 - w) * np.minimum(self.dist2.sf(x - 1), self.dist2.cdf(x))
        # )
        return p_right, p_left, p_both

    def estimate_w(self, x):
        if self.e == 0:
            b = np.log(self.B)
            delta = b * (self.n - 2 * x)
            return expit(-delta)
        return expit(self.dist1.logpmf(x) - self.dist2.logpmf(x))


class CachedBinomialScoringModel(BinomialScoringModel):
    """
    A model with pre-calculated p-values for all possible x values
    """ 
    
    def _calc_all_pvalues(self):
        self._p_right, self._p_left, self._p_both = self.calc_pvalues(self.all_x)
    
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
    
    def get_signif_indices(self, signif_tr, side='both'):
        return self.p_value(side=side) <= signif_tr


class ExactPowerEstimator:
    def __init__(self, null_model: BinomialModel, scoring_model: CachedBinomialScoringModel, bad_phasing_mode=1):
        self.null_model = null_model
        self.scoring_model = scoring_model
        assert bad_phasing_mode in [None, 1, 2]
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


    # For true distribution, knowing the correct BAD phasing mode.
    # By default positive effect corresponds to preference towards copied allele (mode=1)    
    def sum_probability_for_mode(self, model, indicators):
        if self.bad_phasing_mode is None:
            log_pmf_for_mode = logsumexp([model.dist1.logpmf(model.all_x), model.dist2.logpmf(model.all_x)], axis=0) - np.log(2)
        else:
            dist = model.dist1 if self.bad_phasing_mode == 1 else model.dist2
            log_pmf_for_mode = dist.logpmf(model.all_x)
        return np.exp(logsumexp(log_pmf_for_mode[indicators])) if indicators.sum() != 0 else 0

