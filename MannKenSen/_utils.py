"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""

from __future__ import division
import datetime
import numpy as np
from scipy.stats import norm, rankdata
from collections import namedtuple


# Helper functions from the original pymannkendall.py
def __preprocessing(x):
    x = np.asarray(x)

    # Convert datetime objects to numeric timestamps if necessary
    if np.issubdtype(x.dtype, np.datetime64):
        x = x.astype('datetime64[s]').astype(float)
    elif x.dtype == 'O' and len(x) > 0:
        if isinstance(x[0], datetime.datetime):
            x = np.array([val.timestamp() for val in x])

    x = x.astype(float)

    if x.ndim == 2:
        (n, c) = x.shape
        if c == 1:
            x = x.flatten()
    return x, (1 if x.ndim == 1 else x.shape[1])

def __missing_values_analysis(x, method='skip'):
    if method.lower() == 'skip':
        if x.ndim == 1:
            x = x[~np.isnan(x)]
        else:
            x = x[~np.isnan(x).any(axis=1)]
    return x, len(x)

def __mk_score(x, n):
    s = 0
    demo = np.ones(n)
    for k in range(n - 1):
        s += np.sum(demo[k + 1:n][x[k + 1:n] > x[k]]) - np.sum(demo[k + 1:n][x[k + 1:n] < x[k]])
    return s

def __variance_s(x, n):
    unique_x, tp = np.unique(x, return_counts=True)
    if n == len(unique_x):
        return (n * (n - 1) * (2 * n + 5)) / 18
    return (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) * (2 * tp + 5))) / 18

def __z_score(s, var_s):
    if var_s == 0:
        return 0
    if s > 0:
        return (s - 1) / np.sqrt(var_s)
    return (s + 1) / np.sqrt(var_s) if s < 0 else 0

def __p_value(z, alpha):
    p = 2 * (1 - norm.cdf(abs(z)))
    h = abs(z) > norm.ppf(1 - alpha / 2)
    trend = 'decreasing' if z < 0 and h else 'increasing' if z > 0 and h else 'no trend'
    return p, h, trend

def __mk_probability(p, s):
    """
    Computes the Mann-Kendall probability.
    """
    C = 1 - p / 2
    Cd = C if s <= 0 else p / 2
    return C, Cd

def __sens_estimator_unequal_spacing(x, t):
    """
    Computes Sen's slope for unequally spaced data using a vectorized approach.
    """
    n = len(x)
    if n < 2:
        return np.array([])

    # Create all pairs of indices
    i, j = np.triu_indices(n, k=1)

    # Calculate differences
    x_diff = x[j] - x[i]
    t_diff = t[j] - t[i]

    # Avoid division by zero
    valid_mask = t_diff != 0

    return x_diff[valid_mask] / t_diff[valid_mask]

def __confidence_intervals(slopes, var_s, alpha):
    """
    Computes the confidence intervals for Sen's slope.
    """
    n = len(slopes)
    if n == 0 or var_s == 0:
        return np.nan, np.nan

    # For a two-sided confidence interval
    Z = norm.ppf(1 - alpha / 2)

    # Ranks of the lower and upper confidence limits (1-based)
    C = Z * np.sqrt(var_s)
    M1 = (n - C) / 2
    M2 = (n + C) / 2

    sorted_slopes = np.sort(slopes)

    # Interpolate to find the values at the fractional ranks
    ranks = np.arange(1, n + 1)
    lower_ci = np.interp(M1, ranks, sorted_slopes)
    upper_ci = np.interp(M2, ranks, sorted_slopes)

    return lower_ci, upper_ci
