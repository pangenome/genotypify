"beta mixture with spikes model"
from functools import partial
from typing import NamedTuple

import jax
import jax.numpy as jnp
import numpy as np
from jax import jit, tree_map, vmap
from jax.scipy.special import betaln, gammaln, logsumexp, xlog1py, xlogy


def id_print(x, **kw):
    return x


def _logbinom(n, k):
    return gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1)


def _wf_trans(s, N, a, b):
    # X` = Y / N where:
    #
    #     Y ~ Binomial(N, p') | 0 < Y < N ,
    #     p' = p + p(1-p)(s/2),
    #     p ~ Beta(a,b)
    #
    # E(Y | 0 < Y < N) = N(p' - p'^n)
    # EX' = Ep'(1-p'^{N-1})
    a, b = id_print((a, b), what="ab")
    EX = (a * (2 + 2 * a + b * (2 + s))) / (2.0 * (a + b) * (1 + a + b))
    # EX = 0.5 * (
    #     (a * (2 + 2 * a + b * (2 + s))) / (a + b) / (1 + a + b)
    #     - (
    #         (2 * (a + b + N) + b * N * s)
    #         * jnp.exp(
    #             gammaln(a + b) + gammaln(a + N) - gammaln(a) - gammaln(1 + a + b + N)
    #         )
    #     )
    # )
    # var(X') = E var(X'|X) + var E(X' | X) = E p'(1-p')/N + var(x + x(1-x)*(s/2))

    # E(X'|X) = p'(1-p'^(N-1))
    # E p'(1-p') / N
    Evar = (
        (a * b * (4 - a * (-2 + s) + b * (2 + s)))
        / (2.0 * (a + b) * (1 + a + b) * (2 + a + b))
        / N
    )
    varE = (
        a
        * b
        * (
            4 * (1 + a + b) * (2 + a + b) * (3 + a + b)
            - 4 * (a - b) * (1 + a + b) * (3 + a + b) * s
            + (a + a**3 - a**2 * (-2 + b) + b * (1 + b) ** 2 - a * b * (2 + b))
            * s**2
        )
    ) / (4.0 * (a + b) ** 2 * (1 + a + b) ** 2 * (2 + a + b) * (3 + a + b))
    Evar, varE = id_print((Evar, varE), what="Evar/VarE")
    var = Evar + varE
    # EX = E(p')
    # var = E(p'(1-p')/N) + var(p) + (s/2)^2 var(p(1-p)) + s * cov(p, p(1-p))
    u = EX * (1 - EX) / var - 1.0
    u = id_print(u, what="u")
    a1 = u * EX
    b1 = u * (1 - EX)
    a1, b1 = id_print((a1, b1), what="a1b1")
    return a1, b1


class BetaMixture(NamedTuple):
    """Mixture of beta pdfs:

    M = len(c) - 1
    p(x) = sum_{i=0}^{M} c[i] x^(a[i] - 1) (1-x)^(b[i] - 1) / beta(a[i], b[i])
    """

    a: np.ndarray
    b: np.ndarray
    log_c: np.ndarray

    @classmethod
    def uniform(cls, M) -> "BetaMixture":
        return cls.interpolate(lambda x: 1.0, M)

    @classmethod
    def interpolate(cls, f, M, norm=False) -> "BetaMixture":
        # bernstein polynomial basis:
        # sum_{i=0}^(M) f(i/M) binom(M,i) x^i (1-x)^(M-i)
        #   = sum_{i=1}^(M+1) f((i-1)/M) binom(M,i-1) x^(i-1) (1-x)^(M-i+2-1)
        #   = sum_{i=1}^(M+1) f((i-1)/M) binom(M,i-1) * beta(i,M-i+2) x^(i-1) (1-x)^(M-i+2-1) / beta(i, M-i+2)
        #   = sum_{i=1}^(M+1) f((i-1)/M) (1+M)^-1 x^(i-1) (1-x)^(M-i+2-1) / beta(i, M-i+2)
        i = jnp.arange(1, M + 2, dtype=float)
        c = vmap(f)((i - 1) / M) / (1 + M)
        c0 = jnp.isclose(c, 0.0)  # work around nans in gradients
        c_safe = jnp.where(c0, 1.0, c)
        log_c = jnp.where(c0, -jnp.inf, jnp.log(c_safe))
        if norm:
            log_c -= logsumexp(log_c)
        return cls(a=i, b=M - i + 2, log_c=log_c)

    @property
    def M(self):
        return len(self.log_c) - 1

    @property
    def moments(self):
        "return a matrix A such that A @ self.c = [EX, EX^2] where X ~ self"
        a = self.a
        b = self.b
        EX = a / (a + b)  # EX
        EX2 = a * (a + 1) / (a + b) / (1 + a + b)  # EX2
        return jnp.array([EX, EX2 - EX**2])

    def __call__(self, x):
        return np.exp(
            self.log_c
            + xlogy(self.a - 1, x)
            + xlog1py(self.b - 1, -x)
            - betaln(self.a, self.b)
        ).sum()

    def plot(self, K=100, ax=None) -> None:
        if ax is None:
            import matplotlib.pyplot as plt

            ax = plt.gca()

        x = np.linspace(0.0, 1.0, K)
        y = np.vectorize(self)(x)
        ax.plot(x, y)


class SpikedBeta(NamedTuple):
    log_p0: float
    log_p1: float
    f_x: BetaMixture

    def sample_component(self, rng):
        sub1, sub2, sub3 = jax.random.split(rng, 4)
        i = jax.random.categorical(sub1, self.f_x.log_c)
        p = jax.random.beta(sub2, self.f_x.a[i], self.f_x.b[i])
        j = jax.random.categorical(
            sub3, jnp.array([self.log_p0, self.log_p1, self.log_r])
        )
        return jnp.array([0.0, 1.0, p])[j]

    @property
    def log_r(self):
        return jnp.log1p(-jnp.exp(self.log_p0 + self.log_p1))

    @property
    def M(self):
        return self.f_x.M

    def plot(self):
        import matplotlib.pyplot as plt

        self.f_x.plot()
        plt.bar(0.0, self.p0, 0.05, alpha=0.2, color="tab:red")
        plt.bar(1.0, self.p1, 0.05, alpha=0.2, color="tab:red")


def transition(f: SpikedBeta, s: float, Ne: float, n: int, d: int):
    """Given a prior distribution on population allele frequency, compute posterior after
    observing data at dt generations in the future"""
    # var(X') = E var(X'|X) + var E(X' | X)
    #         ~= var(X'|X = EX) + var E(X' | X)
    a = f.f_x.a
    b = f.f_x.b
    log_c = f.f_x.log_c
    # s = id_print(s, what="s")
    # probability mass for fixation. p(beta=0) = p0 + \int_0^1 f(x) (1-x)^n, and similarly for p1
    log_p0 = f.log_p0 + logsumexp(
        log_c + gammaln(a + b) + gammaln(b + Ne) - gammaln(b) - gammaln(a + b + Ne)
    )
    log_p1 = f.log_p1 + logsumexp(
        log_c + gammaln(a + b) + gammaln(a + Ne) - gammaln(a) - gammaln(a + b + Ne)
    )
    # new parameters of each mixture component after w-f sampling
    a1, b1 = _wf_trans(s, Ne, a, b)
    # now model binomial sampling
    # p(d_k|d_{k-1},...,d_1) = \int p(d_k|x_k) p(x_k | d_k-1,..,d1)
    # probability of data arising from each mixing component
    # a1, b1 = id_print((a1, b1), what="wf_trans")
    ret = _binom_sampling(n, d, a1, b1, log_c, log_p0, log_p1)
    ret = id_print(ret, what="ret")
    return ret


def _binom_sampling(n, d, a, b, log_c, log_p0, log_p1):
    log_r = jnp.log1p(-jnp.exp(log_p0 + log_p1))
    log_r, _, _ = id_print((log_r, log_p0, log_p1), what="r/p0/p1")
    a1 = a + d
    b1 = b + n - d
    # probability of the data given each mixing component
    log_p_mix = betaln(a1, b1) - betaln(a, b) + _logbinom(n, d)
    log_p_mix = id_print(log_p_mix, what="p_mix")
    log_p01 = log_c + log_p_mix
    lp0 = log_p0 + jnp.log(d == 0)
    lp1 = log_p1 + jnp.log(d == n)
    ll = logsumexp(jnp.concatenate([lp0[None], log_r + log_p01, lp1[None]]))
    # p(p=1|d) = p(d|p=1)*p01 / pd
    log_p0 = lp0 - ll
    log_p1 = lp1 - ll
    # update mixing coeffs
    # p(c | data) = p(data | c) * p(c) / p(data)
    log_c1 = log_r + log_p01
    log_c1 -= logsumexp(log_c1)
    # posterior after binomial sampling --
    beta = SpikedBeta(log_p0, log_p1, BetaMixture(a1, b1, log_c1))
    # beta = id_print(beta, what="beta")
    return beta, ll


# @partial(jit, static_argnums=3)
def forward(s, Ne, obs, prior):
    """
    Run the forward algorithm for the BMwS model.

    Args:
        s: selection coefficient at each time point (T - 1,)
        Ne:  diploid effective population size at each time point (T - 1,)
        obs: (sample size, # derived alleles) observed at each time point (T, 2)
        times: number of generations before present when each observation was sampled, sorted descending (T,)

    Returns:
        Tuple (betas, lls). betas are the filtering distributions, and lls are the conditional likelihoods.
    """
    n, d = obs.T

    def _f(fX, d):
        beta, ll = transition(fX, **d)
        return beta, (beta, ll)

    # uniform
    f_x = prior  # BetaMixture.uniform(M)

    beta0, ll0 = _binom_sampling(
        n[-1], d[-1], f_x.a, f_x.b, f_x.log_c, -jnp.inf, -jnp.inf
    )

    if False:
        lls = []
        betas = []
        f_x = beta0
        for i in range(1, 1 + len(Ne)):
            f_x, (_, ll) = _f(
                f_x, {"Ne": Ne[-i], "n": n[-i - 1], "d": d[-i - 1], "s": s[-i]}
            )
            betas.append(f_x)
            lls.append(ll)
    else:
        _, (betas, lls) = jax.lax.scan(
            _f, beta0, {"Ne": Ne, "n": n[:-1], "d": d[:-1], "s": s}, reverse=True
        )

    return (betas, beta0), jnp.concatenate([jnp.array(lls), ll0[None]])


def loglik(s, Ne, obs, prior):
    betas, lls = forward(s, Ne, obs, prior)
    return lls.sum()


def _construct_prior(prior: int | BetaMixture) -> BetaMixture:
    if isinstance(prior, int):
        M = prior
        prior = BetaMixture.uniform(M)
    return prior


@partial(jit, static_argnums=(3, 4))
def sample_paths(
    s: np.ndarray,
    Ne: np.ndarray,
    obs: np.ndarray,
    k: int,
    seed: int = 1,
    prior: int | BetaMixture = BetaMixture.uniform(100),
):
    """
    Sample allele frequency paths from posterior distribution.

    Args:
        s: selection coefficient at each time point (T - 1,)
        Ne:  diploid effective population size at each time point (T - 1,)
        obs: (sample size, # derived alleles) observed at each time point (T, 2)
        k: number of paths to sample
        seed: seed for random number generator
        prior: prior on initial allele frequency

    Returns:
        Array of shape (k, T), containing k samples from the allele frequency posterior.

    Notes:
        - obs[0] denotes the most recent observation; obs[-1] is the most ancient.
        - s and Ne control the Wright-Fisher transitions that occur *between* each time points.
          Therefore, there they have one less entry than the number of observations.
    """
    rng = jax.random.PRNGKey(seed)

    prior = _construct_prior(prior)

    def _sample_spikebeta(beta: SpikedBeta, rng, cond):
        # -1, M represent the special fixed states 0/1
        log_p0p = beta.log_p0 + jnp.log(cond != -1)
        log_p1p = beta.log_p1 + jnp.log(cond != -2)
        log_p = jnp.concatenate(
            [log_p0p[None], log_p1p[None], beta.log_r + beta.f_x.log_c]
        )
        sub1, sub2 = jax.random.split(rng)
        s = jax.random.categorical(sub1, log_p) - 2
        x = jnp.where(
            s < 0,
            jnp.take(jnp.array([0.0, 1.0]), 2 + s),
            jax.random.beta(sub2, beta.f_x.a[s], beta.f_x.b[s]),
        )
        return (s, x)

    (betas, beta_n), _ = forward(s, Ne, obs, prior)

    beta0 = tree_map(lambda a: a[0], betas)
    beta1n = tree_map(
        lambda a, b: jnp.concatenate([a[1:], b[None]]), betas, beta_n
    )
    betas = tree_map(lambda a, b: jnp.concatenate([a, b[None]]), betas, beta_n)

    def _f(tup, beta):
        rng, s1 = tup
        rng, sub1, sub2 = jax.random.split(rng, 3)
        s, x = _sample_spikebeta(beta, sub1, s1)
        return (rng, s), x

    def _g(rng, _):
        rng, sub = jax.random.split(rng)
        s0, x0 = _sample_spikebeta(beta0, sub, 0)
        _, xs = jax.lax.scan(_f, (sub, s0), beta1n)
        return rng, jnp.concatenate([x0[None], xs])

    _, ret = jax.lax.scan(_g, rng, None, length=k)
    return ret, betas
