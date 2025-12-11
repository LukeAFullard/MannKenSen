"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
from collections import namedtuple
from ._utils import (__mk_score, __variance_s, __z_score,
                   __p_value, __sens_estimator_unequal_spacing,
                   __confidence_intervals, __mk_probability,
                   _mk_score_and_var_censored, _sens_estimator_censored,
                   _prepare_data)
from .plotting import plot_trend


def original_test(x, t, alpha=0.05, hicensor=False, plot_path=None, lt_mult=0.5, gt_mult=1.1, sens_slope_method='lwp', tau_method='b'):
    """
    Mann-Kendall test for unequally spaced time series.
    Input:
        x: a vector of data, or a DataFrame from prepare_censored_data.
        t: a vector of timestamps corresponding to x.
        alpha: significance level (default 0.05).
        hicensor (bool): If True, applies the high-censor rule, where all
                         values below the highest left-censor limit are
                         treated as censored at that limit.
        plot_path (str, optional): If provided, a plot of the trend analysis
                                   is saved to this file path.
        lt_mult (float): The multiplier for left-censored data (default 0.5).
        gt_mult (float): The multiplier for right-censored data (default 1.1).
        sens_slope_method (str): The method to use for handling ambiguous slopes
                                 in censored data. See `_sens_estimator_censored`
                                 for details.
        tau_method (str): The method for calculating Kendall's Tau ('a' or 'b').
                          Default is 'b', which accounts for ties.
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd

    A namedtuple containing the following fields:
        - trend: The trend of the data ('increasing', 'decreasing', or 'no trend').
        - h: A boolean indicating whether the trend is significant.
        - p: The p-value of the test.
        - z: The Z-statistic.
        - Tau: Kendall's Tau, a measure of correlation between the data and time.
               Ranges from -1 (perfectly decreasing) to +1 (perfectly increasing).
        - s: The Mann-Kendall score.
        - var_s: The variance of `s`.
        - slope: The Sen's slope.
        - intercept: The intercept of the trend line.
        - lower_ci: The lower confidence interval of the slope.
        - upper_ci: The upper confidence interval of the slope.
        - C: The confidence of the trend direction.
        - Cd: The confidence that the trend is decreasing.
    """
    res = namedtuple('Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])

    data_filtered, _ = _prepare_data(x, t, hicensor)

    # Check for tied timestamps
    if len(data_filtered['t']) != len(np.unique(data_filtered['t'])):
        import warnings
        warnings.warn("Tied timestamps detected in the time vector `t`. Corresponding data points will be excluded from the Sen's slope calculation.", UserWarning)

    x_filtered = data_filtered['value'].to_numpy()
    t_filtered = data_filtered['t'].to_numpy()
    censored_filtered = data_filtered['censored'].to_numpy()
    cen_type_filtered = data_filtered['cen_type'].to_numpy()


    if len(x_filtered) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    s, var_s, D = _mk_score_and_var_censored(
        x_filtered, t_filtered, censored_filtered, cen_type_filtered,
        tau_method=tau_method
    )
    Tau = s / D if D > 0 else 0

    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)
    C, Cd = __mk_probability(p, s)

    if np.any(censored_filtered):
        slopes = _sens_estimator_censored(
            x_filtered, t_filtered, cen_type_filtered,
            lt_mult=lt_mult, gt_mult=gt_mult, method=sens_slope_method
        )
    else:
        slopes = __sens_estimator_unequal_spacing(x_filtered, t_filtered)

    slope = np.nanmedian(slopes) if len(slopes) > 0 else np.nan

    if np.isnan(slope):
        intercept = np.nan
    else:
        intercept = np.nanmedian(x_filtered) - np.nanmedian(t_filtered) * slope

    lower_ci, upper_ci = __confidence_intervals(slopes, var_s, alpha)

    results = res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)

    if plot_path:
        plot_trend(data_filtered, results, plot_path, alpha)

    return results
