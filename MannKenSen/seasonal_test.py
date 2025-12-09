"""
This script provides a modified version of the Seasonal Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
from ._utils import (__preprocessing, __missing_values_analysis, __mk_score,
                   __variance_s, __z_score, __p_value,
                   __sens_estimator_unequal_spacing, __confidence_intervals)

def seasonal_test(x_old, t_old, period=12, alpha=0.05):
    """
    Seasonal Mann-Kendall test for unequally spaced time series.
    Input:
        x_old: a vector of data
        t_old: a vector of timestamps
        period: seasonal cycle (default 12)
        alpha: significance level (default 0.05)
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci
    """
    res = namedtuple('Seasonal_Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci'])

    x, _ = __preprocessing(x_old)
    t, _ = __preprocessing(t_old)

    # Filter out NaN values from the start
    mask = ~np.isnan(x) & ~np.isnan(t)
    x = x[mask]
    t = t[mask]

    if len(x) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan)

    # Determine seasons for each data point
    seasons = np.floor(t - 1) % period

    s = 0
    var_s = 0
    denom = 0
    all_slopes = []

    for i in range(int(period)):
        season_mask = seasons == i
        season_x = x[season_mask]
        season_t = t[season_mask]
        n = len(season_x)

        if n > 1:
            s += __mk_score(season_x, n)
            var_s += __variance_s(season_x, n)
            denom += 0.5 * n * (n - 1)
            all_slopes.extend(__sens_estimator_unequal_spacing(season_x, season_t))

    Tau = s / denom if denom != 0 else 0
    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)

    if not all_slopes:
        slope, intercept, lower_ci, upper_ci = np.nan, np.nan, np.nan, np.nan
    else:
        slopes_array = np.asarray(all_slopes)
        slope = np.nanmedian(slopes_array)
        intercept = np.nanmedian(x) - np.nanmedian(t) * slope
        lower_ci, upper_ci = __confidence_intervals(slopes_array, var_s, alpha)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci)
