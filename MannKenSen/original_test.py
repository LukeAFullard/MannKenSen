"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
import pandas as pd
from pandas import DataFrame
from ._utils import __preprocessing, __mk_score, __variance_s, __z_score, __p_value, __sens_estimator_unequal_spacing, __confidence_intervals, __mk_probability

def original_test(x, t, alpha=0.05, hicensor=False):
    """
    Mann-Kendall test for unequally spaced time series.
    Input:
        x: a vector of data, or a DataFrame from prepare_censored_data.
        t: a vector of timestamps corresponding to x.
        alpha: significance level (default 0.05).
        hicensor (bool): If True, applies the high-censor rule, where all
                         values below the highest left-censor limit are
                         treated as censored at that limit.
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd
    """
    res = namedtuple('Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])

    # Input validation and preparation
    if isinstance(x, DataFrame) and all(col in x.columns for col in ['value', 'censored', 'cen_type']):
        data = x.copy()
    elif hasattr(x, '__iter__') and any(isinstance(i, str) for i in x):
        raise TypeError("Input data `x` contains strings. Please pre-process it with `prepare_censored_data` first.")
    else:
        x_proc, _ = __preprocessing(x)
        data = pd.DataFrame({
            'value': x_proc,
            'censored': np.zeros(len(x_proc), dtype=bool),
            'cen_type': np.full(len(x_proc), 'not', dtype=object)
        })

    t_proc, _ = __preprocessing(t)
    data['t'] = t_proc

    # Handle missing values
    mask = ~np.isnan(data['value']) & ~np.isnan(data['t'])
    data_filtered = data[mask].copy()

    # Apply HiCensor rule if requested
    if hicensor and 'lt' in data_filtered['cen_type'].values:
        max_lt_censor = data_filtered.loc[data_filtered['cen_type'] == 'lt', 'value'].max()
        hi_censor_mask = data_filtered['value'] < max_lt_censor
        data_filtered.loc[hi_censor_mask, 'censored'] = True
        data_filtered.loc[hi_censor_mask, 'cen_type'] = 'lt'
        data_filtered.loc[hi_censor_mask, 'value'] = max_lt_censor

    x_filtered = data_filtered['value'].to_numpy()
    t_filtered = data_filtered['t'].to_numpy()
    censored_filtered = data_filtered['censored'].to_numpy()
    cen_type_filtered = data_filtered['cen_type'].to_numpy()


    if len(x_filtered) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    if np.any(censored_filtered):
        s, var_s, D = _mk_score_and_var_censored(x_filtered, t_filtered, censored_filtered, cen_type_filtered)
        Tau = s / D if D > 0 else 0
    else:
        s = __mk_score(x_filtered, len(x_filtered))
        var_s = __variance_s(x_filtered, len(x_filtered))
        Tau = s / (0.5 * len(x_filtered) * (len(x_filtered) - 1))

    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)
    C, Cd = __mk_probability(p, s)

    if np.any(censored_filtered):
        slopes = _sens_estimator_censored(x_filtered, t_filtered, cen_type_filtered)
    else:
        slopes = __sens_estimator_unequal_spacing(x_filtered, t_filtered)

    slope = np.nanmedian(slopes) if len(slopes) > 0 else np.nan

    if np.isnan(slope):
        intercept = np.nan
    else:
        intercept = np.nanmedian(x_filtered) - np.nanmedian(t_filtered) * slope

    lower_ci, upper_ci = __confidence_intervals(slopes, var_s, alpha)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)
