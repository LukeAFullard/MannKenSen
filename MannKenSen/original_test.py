"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
from ._utils import __preprocessing, __missing_values_analysis, __mk_score, __variance_s, __z_score, __p_value, sens_slope

def original_test(x_old, t_old, alpha=0.05):
    """
    Mann-Kendall test for unequally spaced time series.
    Input:
        x_old: a vector of data
        t_old: a vector of timestamps corresponding to x_old
        alpha: significance level (default 0.05)
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept
    """
    res = namedtuple('Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept'])
    x, _ = __preprocessing(x_old)

    x_mk, n = __missing_values_analysis(x)

    if n < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan)

    s = __mk_score(x_mk, n)
    var_s = __variance_s(x_mk, n)
    Tau = s / (0.5 * n * (n - 1))

    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)

    slope, intercept = sens_slope(x_old, t_old)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept)
