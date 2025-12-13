"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
import pandas as pd
import warnings
from collections import namedtuple
from ._stats import (_z_score, _p_value, _sens_estimator_unequal_spacing,
                     _confidence_intervals, _mk_probability,
                     _mk_score_and_var_censored, _sens_estimator_censored)
from ._helpers import (_prepare_data, _aggregate_by_group)
from .plotting import plot_trend


def original_test(x, t, alpha=0.05, hicensor=False, plot_path=None, lt_mult=0.5, gt_mult=1.1, sens_slope_method='nan', tau_method='b', agg_method='none', min_size=10):
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
        sens_slope_method (str): The method for handling ambiguous slopes in censored data.
            - 'nan' (default): Sets ambiguous slopes (e.g., between two left-censored
                               values) to `np.nan`, effectively removing them from the
                               median slope calculation. This is a statistically neutral
                               approach.
            - 'lwp': Sets ambiguous slopes to 0, mimicking the LWP-TRENDS R script.
                     This may bias the slope towards zero.
        tau_method (str): The method for calculating Kendall's Tau ('a' or 'b').
                          Default is 'b', which accounts for ties.
        agg_method (str): The method for aggregating data at tied timestamps.
            - 'none' (default): No aggregation. A warning is issued if ties are present.
            - 'median': Use the median of values and times. For censored data, this
                        is a simple heuristic.
            - 'robust_median': A more statistically robust median for censored data.
            - 'middle': Use the observation closest to the middle of the time period.
        min_size (int): Minimum sample size. Warnings issued if n < min_size.
                       Set to None to disable check.
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

    Statistical Assumptions:
    ----------------------
    The Mann-Kendall test and Sen's slope estimator are non-parametric methods
    and do not require data to be normally distributed. However, they rely on
    the following assumptions:

    1.  **Independence**: The data points are serially independent. The presence
        of autocorrelation (serial correlation) can violate this assumption and
        affect the significance of the results.
    2.  **Monotonic Trend**: The trend is assumed to be monotonic, meaning it is
        consistently in one direction (either increasing or decreasing) over
        the time period. The test is not suitable for detecting non-monotonic
        (e.g., cyclical) trends.
    3.  **Homogeneity of Variance**: While not strictly required, the test is
        most powerful when the variance of the data is homogeneous over time.
    4.  **Continuous Data**: The test is designed for continuous data. Although it
        can handle ties, a large number of ties can reduce its power.
    5.  **Normal Approximation**: For sample sizes typically > 10, the test
        statistic `S` is assumed to be approximately normally distributed, which
        is used to calculate the Z-score and p-value. This approximation may be
        less accurate for very small sample sizes with many ties.
    6.  **Linear Trend (for Sen's Slope)**: Sen's slope provides a linear estimate
        of the trend, and the confidence intervals are based on this assumption.
    """
    res = namedtuple('Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])

    # --- Input Validation ---
    valid_sens_slope_methods = ['nan', 'lwp']
    if sens_slope_method not in valid_sens_slope_methods:
        raise ValueError(f"Invalid `sens_slope_method`. Must be one of {valid_sens_slope_methods}.")

    valid_tau_methods = ['a', 'b']
    if tau_method not in valid_tau_methods:
        raise ValueError(f"Invalid `tau_method`. Must be one of {valid_tau_methods}.")

    valid_agg_methods = ['none', 'median', 'robust_median', 'middle']
    if agg_method not in valid_agg_methods:
        raise ValueError(f"Invalid `agg_method`. Must be one of {valid_agg_methods}.")

    data_filtered, is_datetime = _prepare_data(x, t, hicensor)

    n = len(data_filtered)

    # Sample size validation
    if n < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan,
                  np.nan, np.nan, np.nan, np.nan)

    if min_size is not None and n < min_size:
        import warnings
        warnings.warn(
            f"Sample size (n={n}) is below recommended minimum (n={min_size}). "
            f"Results may be unreliable. Consider using more data or setting "
            f"min_size=None to suppress this warning.",
            UserWarning
        )

    # Handle tied timestamps
    if len(data_filtered['t']) != len(np.unique(data_filtered['t'])):
        if agg_method == 'none':
            warnings.warn(
                "Tied timestamps detected in the time vector `t`. "
                "The Sen's slope calculation may be affected. "
                "Consider using an aggregation method via the `agg_method` parameter.",
                UserWarning
            )
        else:
            agg_data_list = [
                _aggregate_by_group(group, agg_method, is_datetime)
                for _, group in data_filtered.groupby('t')
            ]
            data_filtered = pd.concat(agg_data_list, ignore_index=True)

    x_filtered = data_filtered['value'].to_numpy()
    t_filtered = data_filtered['t'].to_numpy()
    censored_filtered = data_filtered['censored'].to_numpy()
    cen_type_filtered = data_filtered['cen_type'].to_numpy()


    if len(x_filtered) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    s, var_s, D, Tau = _mk_score_and_var_censored(
        x_filtered, t_filtered, censored_filtered, cen_type_filtered,
        tau_method=tau_method
    )

    z = _z_score(s, var_s)
    p, h, trend = _p_value(z, alpha)
    C, Cd = _mk_probability(p, s)

    if np.any(censored_filtered):
        slopes = _sens_estimator_censored(
            x_filtered, t_filtered, cen_type_filtered,
            lt_mult=lt_mult, gt_mult=gt_mult, method=sens_slope_method
        )
    else:
        slopes = _sens_estimator_unequal_spacing(x_filtered, t_filtered)

    slope = np.nanmedian(slopes) if len(slopes) > 0 else np.nan

    if np.isnan(slope):
        intercept = np.nan
    else:
        intercept = np.nanmedian(x_filtered) - np.nanmedian(t_filtered) * slope

    lower_ci, upper_ci = _confidence_intervals(slopes, var_s, alpha)

    results = res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)

    if plot_path:
        plot_trend(data_filtered, results, plot_path, alpha)

    return results
