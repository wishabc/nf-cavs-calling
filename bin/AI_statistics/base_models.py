import numpy as np
from functools import wraps
from scipy.special import logsumexp
import scipy.stats as st


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


def hash_element(element):
    if isinstance(element, np.ndarray):
        return hash(element.tobytes())
    return hash(element)


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
            hash_element(arg)
            for arg in args
        )
        hashed_kwargs = tuple(
            (key, hash_element(value))
            for key, value in kwargs.items()
        )
        key = (hashed_args, hashed_kwargs)
        if key not in cache:
            cache[key] = func(self, *args, **kwargs)
        return cache[key]

    return wrapper


class BimodalEffectModel:
    """
    Base class to store parameters of a model
    """
    def __init__(self, n, effect, B, **kwargs):
        self.e = effect
        self.n = n
        self.B = B

    
    @property
    def all_observations(self):
        return np.arange(self.n + 1)
    
    @property
    def dist1(self) -> st.rv_discrete:
        raise NotImplementedError
    
    @property
    def dist2(self) -> st.rv_discrete:
        raise NotImplementedError
    
    @classmethod
    def from_model(cls, other: 'BimodalEffectModel', effect=None) -> 'BimodalEffectModel':
        raise NotImplementedError
    
    @cached_method
    def get_effect_model(self, effect):
        """
        Get a model with a different effect size
        """
        if effect == self.e:
            return self
        return self.__class__.from_model(self, effect=effect)
    
    def get_log_pmf_for_mode(self, bad_phasing_mode):
        assert bad_phasing_mode in [1, 2, None]
        if bad_phasing_mode is None:
            log_pmf_for_mode = logsumexp([self.dist1.logpmf(self.all_observations), self.dist2.logpmf(self.all_observations)], axis=0) - np.log(2)
        else:
            dist = self.dist1 if bad_phasing_mode == 1 else self.dist2
            log_pmf_for_mode = dist.logpmf(self.all_observations)
        return log_pmf_for_mode
    
    def compatible_with(self, other: 'BimodalEffectModel'):
        return self.n == other.n


class BimodalSamplingModel(BimodalEffectModel):
    """
    A model to simulate allelic imbalance data
    """
    def get_samples(self, **kwargs):
        raise NotImplementedError


class BimodalScoringModel(BimodalEffectModel):
    """
    Base class that contains calc_pvalues method
    """
    def calc_pvalues(self, observations):
        raise NotImplementedError

    def effect_size_estimate(self, observations):
        raise NotImplementedError
