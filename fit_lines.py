# -*- coding: utf-8 -*-
"""
Created on Tue Jul 29 15:57:06 2014

@author: michael
"""

import numpy as np
import emcee as mc
import matplotlib.pyplot as plt


def gendata(n=20, xmin=0., xmax=10., m=1., b=0., xsigma=1.0, ysigma=1.0):
    """
    Generate synthetic data set suitable for linear regression tests.
    Noise is uncorrelated Gaussian in both variables.

    Arguments:
    n = number of data points,
    xmin, xmax = minimum and maximum x value (inclusive)
    m = true slope
    b = true intercept
    xsigma = std.dev. of Gaussian noise in x variable
    ysigma = std.dev. of Gaussian noise in y variable

    Returns (xs, ys) where ys = b + m * xs + noise.
    """
    xs = np.linspace(xmin, xmax, n)  # "true" x values
    ys = b + m*xs + ysigma * np.random.randn(n)  # ys correspond to true xs
    xs += xsigma * np.random.randn(n)  # now add noise to xs
    return (xs, ys)


def lstsqfit(xs, ys):
    """
    Given xs and ys return the (slope, intercept) of a linear least squares
    fit. xs and ys must have the same shape.
    """
    assert xs.shape == ys.shape
    A = np.vstack([xs, np.ones_like(xs)]).T
    m, b = np.linalg.lstsq(A, ys)[0]
    return (m, b)


def lnlike(theta, xs, ys):
    """
    Log-likelihood of data model.

    Arguments:
    theta = [m, b, sigmax, sigmay, xt[0], xt[1], ...] = parameters
    xs = observed x values
    ys = observed y values

    Returns:
    Unnormalized (base-e) log likelihood of the observed data given model
    parameters
    """
    from numpy import dot

    assert len(theta) == len(xs) + 4
    m = theta[0]
    b = theta[1]
    sigmax = theta[2]
    sigmay = theta[3]
    xt = np.array(theta[4:])

    sigmapart = -np.log(sigmax) - np.log(sigmay)
    dxss = dot(xs, xs)
    dxst = dot(xs, xt)
    dxtt = dot(xt, xt)
    dyss = dot(ys, ys)
    dxtys = dot(xt, ys)
    sumx = np.sum(xt)
    sumy = np.sum(ys)
    xpart = -(dxss - 2. * dxst + dxtt)/(2. * sigmax ** 2.)
    ypart = -(dyss - 2. * m * dxtys + m * m * dxtt
              - 2. * b * sumy + 2. * b * m * sumx
              + len(xs) * b * b) / (2. * sigmay ** 2.)
    return sigmapart + xpart + ypart


def lnprior(theta):
    """
    Log-prior distribution

    Arguments:
    theta = [m, b, sigmax, sigmay, xt[0], xt[1], ...] = parameters

    Returns:
    Unnormalized (base-e) log prior of the observed data given model parameters
    """
    assert len(theta) == len(xs) + 4
    m = theta[0]
    b = theta[1]
    sigmax = theta[2]
    sigmay = theta[3]
    xt = theta[4:]
    # Note that the non-informative prior in m is pi(m) = 0.5 * (1+m**2)^(-3/2)
    # The priors for sigmax and sigmay are uniform over a range
    if (-100 < b < 100 and 0.08 < sigmax < 0.12 and 0 < sigmay < 2 and
            min(xt) > -100 and max(xt) < 100):
        return -1.5 * np.log(1. + m ** 2.)
    return -np.inf


def lnprob(theta, xs, ys):
    """
    Log-probability

    Arguments:
    theta = [m, b, sigmax, sigmay, xt[0], xt[1], ...] = parameters
    xs = observed x values
    ys = observed y values

    Returns:
    Unnormalized (base-e) log posterior probability of the observed data given
    model parameters and prior distribution.
    """
    lp = lnprior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lnlike(theta, xs, ys) + lp


#%% Main program - (spyder cell magic)
if __name__ == "__main__":
    import time

    np.random.seed(42)  # reproducibility
    (xsigma, ysigma) = (0.1, 1.2)  # true values
    xs, ys = gendata(n=10, xsigma=xsigma, ysigma=ysigma)
    (fitm, fitb) = lstsqfit(xs, ys)

    print "Data"
    print "    x\t    y"
    print "-------------"
    for i in xrange(len(xs)):
        print "% 4.2f\t% 5.1f" % (xs[i], ys[i])
    print

    # initial guess for MCMC parameters
    t = [1., 0., xsigma, ysigma]
    t.extend(xs)

    ndim = len(t)
    nwalkers = 3 * ndim
    burnin = 200
    nsteps = 10000
    # starting points for walkers in a small ball around our initial guess
    pos = [t + 1e-4 * np.random.randn(ndim) for i in range(nwalkers)]

    # MCMC sampler
    print "Running MCMC chain with %d walkers and %d steps with %d burnin." % \
        (nwalkers, nsteps, burnin)
    start = time.clock()
    sampler = mc.EnsembleSampler(nwalkers, ndim, lnprob, args=(xs, ys))
    sampler.run_mcmc(pos, nsteps)
    # remove burn-in samples and flatten chain
    samples = sampler.chain[:, burnin:, :].reshape((-1, ndim))
    stop = time.clock()
    print "Finished. Took %f seconds." % (stop - start)

#%% Plotting - (spyder cell magic)
    plt.errorbar(xs, ys, yerr=ysigma, xerr=xsigma, fmt='k.')
    plt.plot(xs, fitm*xs + fitb, 'b-')
    lines = samples[np.random.random_integers(0, len(samples) - 1, 50), :2]
    for (m, b) in lines:
        plt.plot(xs, m * xs + b, 'k-', alpha=0.1)
    plt.axis([min(xs)-1, max(xs)+1, min(ys)-1, max(ys)+1])
    (mmin, bmin) = samples[samples[:, 0] == np.min(samples[:, 0])][0, :2]
    (mmax, bmax) = samples[samples[:, 0] == np.max(samples[:, 0])][0, :2]
    plt.plot(xs, mmin * xs + bmin, 'g:', xs, mmax * xs + bmax, 'g:')
    plt.plot(xs, 1.0 * xs + 0., 'r-')
    plt.show()

    import triangle
    fig = triangle.corner(samples[:, :4],
                          labels=["m", "b", "$\sigma_x$", "$\sigma_y$"],
                          truths=[1., 0., xsigma, ysigma],
                          quantiles=[0.159, 0.5, 0.841])  # 1-sigma quantiles
    fig.show()
    fig2 = triangle.corner(samples[:, :2],
                           labels=["m", "b"],
                           truths=[1., 0.],
                           quantiles=[0.159, 0.5, 0.841])  # 1-sigma quantiles
    fig2.show()
