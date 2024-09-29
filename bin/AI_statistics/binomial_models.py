import scipy.stats as st
from scipy.special import expit
import numpy as np
from base_models import BimodalEffectModel, cached_method
from estimate_mse import es_estimate_vectorized, estimate_w_null, es_variance_vectorized


class BinomialModel(BimodalEffectModel):
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
        super.__init__(n, effect, B)
        self.p1 = expit(self.e * np.log(2) + np.log(B))
        self.p2 = expit(self.e * np.log(2) - np.log(B))

        self.dist1 = st.binom(n, self.p1)
        self.dist2 = st.binom(n, self.p2)
    
    @classmethod
    def from_model(cls, other: 'BinomialModel', effect=None):
        if effect is None:
            effect = other.e
        return cls(other.n, effect, other.B)
    

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
            return estimate_w_null(x, self.n, self.B)
        return expit(self.dist1.logpmf(x) - self.dist2.logpmf(x))

    def effect_size_estimate(self, x):
        w = self.estimate_w(x)
        return es_estimate_vectorized(x, self.n, self.B, w)