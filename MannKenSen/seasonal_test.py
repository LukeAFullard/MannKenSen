"""
This script provides a modified version of the Seasonal Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
from ._utils import __preprocessing, __missing_values_analysis, __mk_score, __variance_s, __z_score, __p_value, __sens_estimator_unequal_spacing

def seasonal_sens_slope(x_old, t_old, period=12):
    res = namedtuple('Seasonal_Sens_Slope_Test', ['slope', 'intercept'])
    x, _ = __preprocessing(x_old)
    t, _ = __preprocessing(t_old)

    # Pad and reshape data and time vectors
    if x.ndim == 1:
        if np.mod(len(x), period) != 0:
            pad_len = period - np.mod(len(x), period)
            x = np.pad(x, (0, pad_len), 'constant', constant_values=(np.nan,))
            t = np.pad(t, (0, pad_len), 'constant', constant_values=(np.nan,))

        x = x.reshape(-1, period)
        t = t.reshape(-1, period)

    d = []
    for i in range(period):
        season_x = x[:, i]
        season_t = t[:, i]

        mask = ~np.isnan(season_x) & ~np.isnan(season_t)
        season_x, season_t = season_x[mask], season_t[mask]

        if len(season_x) > 1:
            d.extend(__sens_estimator_unequal_spacing(season_x, season_t))

    if not d:
        return res(np.nan, np.nan)

    slope = np.nanmedian(np.asarray(d))

    mask = ~np.isnan(x_old) & ~np.isnan(t_old)
    valid_x = x_old[mask]
    valid_t = t_old[mask]
    intercept = np.nanmedian(valid_x) - np.nanmedian(valid_t) * slope

    return res(slope, intercept)

def multivariate_test(x_old, t_old, alpha=0.05, period=12):
    res = namedtuple('Multivariate_Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept'])
    x, c = __preprocessing(x_old)

    s = 0
    var_s = 0
    denom = 0

    for i in range(c):
        season_data, n = __missing_values_analysis(x[:, i] if c > 1 else x)
        if n > 1:
            s += __mk_score(season_data, n)
            var_s += __variance_s(season_data, n)
            denom += 0.5 * n * (n - 1)

    Tau = s / denom if denom != 0 else 0
    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)

    slope, intercept = seasonal_sens_slope(x_old.flatten(), t_old.flatten(), period=period)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept)

def seasonal_test(x_old, t_old, period=12, alpha=0.05):
    """
    Seasonal Mann-Kendall test for unequally spaced time series.
    Input:
        x_old: a vector of data
        t_old: a vector of timestamps
        period: seasonal cycle (default 12)
        alpha: significance level (default 0.05)
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept
    """
    res = namedtuple('Seasonal_Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept'])
    x, _ = __preprocessing(x_old)
    t, _ = __preprocessing(t_old)

    if x.ndim == 1:
        if np.mod(len(x), period) != 0:
            pad_len = period - np.mod(len(x), period)
            x = np.pad(x, (0, pad_len), 'constant', constant_values=(np.nan,))
            t = np.pad(t, (0, pad_len), 'constant', constant_values=(np.nan,))
        x = x.reshape(-1, period)

    trend, h, p, z, Tau, s, var_s, slope, intercept = multivariate_test(x, t, alpha=alpha, period=period)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept)
