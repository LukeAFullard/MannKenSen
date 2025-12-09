"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
from ._utils import __preprocessing, __missing_values_analysis, __mk_score, __variance_s, __z_score, __p_value, __sens_estimator_unequal_spacing, __confidence_intervals, __mk_probability

def original_test(x_old, t_old, alpha=0.05):
    """
    Mann-Kendall test for unequally spaced time series.
    Input:
        x_old: a vector of data
        t_old: a vector of timestamps corresponding to x_old
        alpha: significance level (default 0.05)
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd
    """
    res = namedtuple('Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])
    x, _ = __preprocessing(x_old)
    t, _ = __preprocessing(t_old)

    # Handle missing values for slope calculation
    mask = ~np.isnan(x) & ~np.isnan(t)
    x_filtered = x[mask]
    t_filtered = t[mask]

    if len(x_filtered) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    s = __mk_score(x_filtered, len(x_filtered))
    var_s = __variance_s(x_filtered, len(x_filtered))
    Tau = s / (0.5 * len(x_filtered) * (len(x_filtered) - 1))

    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)
    C, Cd = __mk_probability(p, s)

    slopes = __sens_estimator_unequal_spacing(x_filtered, t_filtered)
    slope = np.nanmedian(slopes) if len(slopes) > 0 else np.nan

    if np.isnan(slope):
        intercept = np.nan
    else:
        intercept = np.nanmedian(x_filtered) - np.nanmedian(t_filtered) * slope

    lower_ci, upper_ci = __confidence_intervals(slopes, var_s, alpha)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)
