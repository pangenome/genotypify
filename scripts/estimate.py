import logging
from functools import partial

import jax.numpy as jnp
import numpy as np
import scipy.optimize
import scipy.stats
from jax import jit, lax, value_and_grad
from jax.example_libraries.optimizers import adagrad

import bmws.betamix
from bmws.betamix import BetaMixture

logger = logging.getLogger(__name__)


def _obj(s, Ne, obs, prior, lam):
    ll = bmws.betamix.loglik(s, Ne, obs, prior)
    return -ll + lam * (jnp.diff(s) ** 2).sum()


obj = jit(value_and_grad(_obj))


@partial(jit, static_argnums=(3, 4))
def empirical_bayes(s, obs, Ne, M, num_steps=100, learning_rate=1.0) -> BetaMixture:
    "maximize marginal likelihood w/r/t prior hyperparameters"
    opt_init, opt_update, get_params = adagrad(learning_rate)
    params = jnp.zeros(2)
    opt_state = opt_init(params)

    def beta_pdf(x, a, b):
        from jax.scipy.special import betaln, xlog1py, xlogy

        # assumes a, b >= 1.
        x01 = jnp.isclose(x, 0.0) | jnp.isclose(x, 1.0)
        safe_x = jnp.where(x01, 0.5, x)
        return jnp.where(
            x01,
            0.0,
            jnp.exp(xlogy(a - 1, safe_x) + xlog1py(b - 1, -safe_x) - betaln(a, b)),
        )

    def interp(a, b):
        return BetaMixture.interpolate(lambda x: beta_pdf(x, a, b), M, norm=True)

    def loss_fn(log_ab):
        a, b = 1.0 + jnp.exp(log_ab)
        # a, b = id_print((a, b), what="ab")
        prior = interp(a, b)
        return _obj(s, Ne, obs, prior, 0.0)

    def step(i, opt_state):
        value, grads = value_and_grad(loss_fn)(get_params(opt_state))
        # value, grads = id_print((value, grads), what="vg")
        opt_state = opt_update(i, grads, opt_state)
        return opt_state

    opt_state = lax.fori_loop(0, num_steps, step, opt_state)
    # for i in range(num_steps):
    #     opt_state = step(i, opt_state)
    a_star, b_star = 1.0 + jnp.exp(get_params(opt_state))
    return interp(a_star, b_star)


@partial(jit, static_argnums=5)
def jittable_estimate(obs, Ne, lam, prior, learning_rate=0.1, num_steps=100):
    opt_init, opt_update, get_params = adagrad(learning_rate)
    params = jnp.zeros(len(Ne))
    opt_state = opt_init(params)

    @value_and_grad
    def loss_fn(s):
        return _obj(s, Ne, obs, prior, lam)

    def step(i, opt_state):
        value, grads = loss_fn(get_params(opt_state))
        opt_state = opt_update(i, grads, opt_state)
        return opt_state

    opt_state = lax.fori_loop(0, num_steps, step, opt_state)

    return get_params(opt_state)


def estimate_em(
    obs,
    Ne,
    em_iterations: int = 3,
    lam: float = 1.0,
    solver_options: dict = {},
):
    M = 100
    s = np.zeros(len(obs) - 1)
    for i in range(em_iterations):
        prior = empirical_bayes(s, obs, Ne, M)
        s = estimate(obs, Ne, lam=lam, prior=prior, solver_options=solver_options)

    return s, prior


def estimate(
    obs,
    Ne,
    lam: float = 1.0,
    prior: int | BetaMixture = 100,
    solver_options: dict = {},
):
    args = (Ne, obs, prior)
    x0 = np.zeros(len(obs) - 1)
    kwargs = {}
    optimizer = "L-BFGS-B"
    args += (lam,)
    kwargs["bounds"] = [[-0.2, 0.2]] * len(x0)
    i = 0

    def shim(x, *args):
        nonlocal i
        ret = tuple(np.array(y, dtype=np.float64) for y in obj(x, *args))
        logger.debug("i=%d f=%f |df|=%f", i, ret[0], np.linalg.norm(ret[1]))
        i += 1
        return ret

    # from jax.experimental import optimizers
    #
    # opt_init, opt_update, get_params = optimizers.adagrad(
    #     solver_options.get("learning_rate", 1.0)
    # )
    # opt_state = opt_init(x0)
    #
    # def step(step, opt_state):
    #     value, grads = obj(get_params(opt_state), *args)
    #     opt_state = opt_update(step, grads, opt_state)
    #     return value, opt_state, grads
    #
    # for i in range(1000):
    #     value, opt_state, grads = step(i, opt_state)
    #     logger.debug("f=%f |grad|=%f", value, np.linalg.norm(grads))
    #
    # logger.debug("value=%f grad=%s", value, grads)
    #
    # return get_params(opt_state)

    res = scipy.optimize.minimize(
        shim,
        x0,
        jac=True,
        args=args,
        method=optimizer,
        options=solver_options,
        **kwargs
    )
    logger.debug("Optimization result: %s", res)
    return res.x


def _prep_data(data):
    times, Ne, obs = zip(
        *sorted([(ob.t, ob.Ne, (ob.sample_size, ob.num_derived)) for ob in data])[::-1]
    )
    if len(times) != len(set(times)):
        raise ValueError("times should be distinct")
    obs = np.array(obs)
    assert np.all(obs[:, 1] <= obs[:, 0])
    assert np.all(obs[:, 1] >= 0)
    return np.array(Ne), obs, times
