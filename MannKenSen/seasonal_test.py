"""
This script provides a modified version of the Seasonal Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
import pandas as pd
from pandas import DataFrame
from collections import namedtuple
import warnings
from ._stats import (_z_score, _p_value,
                   _sens_estimator_unequal_spacing, _confidence_intervals,
                   _mk_probability, _mk_score_and_var_censored,
                   _sens_estimator_censored)
from ._datetime import (_get_season_func, _get_cycle_identifier, _get_time_ranks)
from ._helpers import (_prepare_data, _aggregate_by_group)
from .plotting import plot_trend


def seasonal_test(x, t, period=12, alpha=0.05, agg_method='none', season_type='month', hicensor=False, plot_path=None, lt_mult=0.5, gt_mult=1.1, sens_slope_method='nan', tau_method='b', time_method='absolute', min_size_per_season=5):
    """
    Seasonal Mann-Kendall test for unequally spaced time series.
    Input:
        x: a vector of data, or a DataFrame from prepare_censored_data.
        t: a vector of timestamps.
        period: seasonal cycle (default 12).
        alpha: significance level (default 0.05).
        hicensor (bool): If True, applies the high-censor rule, where all
                         values below the highest left-censor limit are
                         treated as censored at that limit.
        plot_path (str, optional): If provided, saves a plot of the trend
                                   analysis to this file path.
        agg_method: method for aggregating multiple data points within a season-year.
                    'none' (default): performs analysis on all data points.
                    'median': (LWP method) uses the median of values and times.
                              For censored data, this is a simple heuristic.
                    'robust_median': uses a more statistically robust median for
                                     censored data. Note: The logic to determine
                                     if the result is censored is a heuristic
                                     from the LWP-TRENDS R script and may not be
                                     universally robust.
                    'middle': uses the observation closest to the middle of the
                              time period.
        season_type: For datetime inputs, specifies the type of seasonality.
                     'year', 'month', 'day_of_week', 'quarter', 'hour', 'week_of_year',
                     'day_of_year', 'minute', 'second'.
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
        time_method (str): The method for handling timestamps in the seasonal test.
                           'absolute' (default): Uses the precise numeric timestamps.
                                                 This is statistically robust for
                                                 unequally spaced data.
                           'rank': Uses cycle-based ranks (1, 2, 3,...) for time,
                                   matching the LWP-TRENDS R script's methodology.
                                   This may be useful for result replication.
        min_size_per_season (int): Minimum observations per season.
                                   Warnings issued if any season < this.
    Output:
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
    The Seasonal Mann-Kendall test extends the standard test by accounting for
    seasonality. It relies on the following assumptions:

    1.  **Independent Seasons**: The trend is analyzed for each season
        independently, and the results are then combined. This assumes that the
        data from different seasons are independent.
    2.  **Serial Independence within Seasons**: The data points within each
        season are assumed to be serially independent.
    3.  **Monotonic Trend per Season**: The test assumes a monotonic trend
        within each season, but the direction and magnitude of the trend can
        vary between seasons.
    4.  **Consistent Seasonal Definition**: The definition of seasons (e.g.,
        'month', 'quarter') must be appropriate for the data and consistent
        throughout the time series.
    5.  **Homogeneity of Trend**: The combined test statistic assumes that the
        trends in each season are homogeneous (i.e., in the same direction).
        If some seasons have increasing trends while others have decreasing
        trends, the test may fail to detect a significant overall trend.
    """
    res = namedtuple('Seasonal_Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])

    # --- Input Validation ---
    valid_agg_methods = ['none', 'median', 'robust_median', 'middle']
    if agg_method not in valid_agg_methods:
        raise ValueError(f"Invalid `agg_method`. Must be one of {valid_agg_methods}.")

    valid_sens_slope_methods = ['nan', 'lwp']
    if sens_slope_method not in valid_sens_slope_methods:
        raise ValueError(f"Invalid `sens_slope_method`. Must be one of {valid_sens_slope_methods}.")

    valid_tau_methods = ['a', 'b']
    if tau_method not in valid_tau_methods:
        raise ValueError(f"Invalid `tau_method`. Must be one of {valid_tau_methods}.")

    valid_time_methods = ['absolute', 'rank']
    if time_method not in valid_time_methods:
        raise ValueError(f"Invalid `time_method`. Must be one of {valid_time_methods}.")

    data_filtered, is_datetime = _prepare_data(x, t, hicensor)

    if is_datetime:
        season_func = _get_season_func(season_type, period)

    if len(data_filtered) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    # --- Aggregation Logic ---
    if agg_method != 'none':
        if is_datetime:
            t_pd = pd.to_datetime(data_filtered['t_original'])
            cycles = _get_cycle_identifier(t_pd, season_type)
            seasons_agg = season_func(t_pd) if season_type != 'year' else np.ones(len(t_pd))
        else:
            t_numeric_agg = data_filtered['t'].to_numpy()
            t_normalized = t_numeric_agg - t_numeric_agg[0]
            cycles = np.floor(t_normalized / period)
            seasons_agg = np.floor(t_normalized % period)

        data_filtered['cycle'] = cycles
        data_filtered['season_agg'] = seasons_agg

        agg_data_list = [
            _aggregate_by_group(group, agg_method, is_datetime)
            for _, group in data_filtered.groupby(['cycle', 'season_agg'])
        ]
        data_filtered = pd.concat(agg_data_list, ignore_index=True)


    # --- Trend Analysis ---
    if is_datetime and season_type != 'year':
        t_pd = pd.to_datetime(data_filtered['t_original'])
        seasons = season_func(t_pd)
        cycles = _get_cycle_identifier(t_pd, season_type)
        season_range = np.unique(seasons)
    elif not is_datetime:
        t_normalized = data_filtered['t'] - data_filtered['t'].min()
        seasons = (np.floor(t_normalized) % period).astype(int)
        cycles = np.floor(t_normalized / period)
        season_range = range(int(period))
    else: # is_datetime and season_type == 'year'
        seasons = np.ones(len(data_filtered))
        cycles = _get_cycle_identifier(pd.to_datetime(data_filtered['t_original']), season_type) if is_datetime else np.zeros(len(data_filtered))
        season_range = [1]


    data_filtered['season'] = seasons
    data_filtered['cycle'] = cycles

    # Sample size validation per season
    if min_size_per_season is not None:
        season_counts = data_filtered.groupby('season').size()
        min_season_n = season_counts.min()

        if min_season_n < min_size_per_season:
            import warnings
            warnings.warn(
                f"Minimum season size (n={min_season_n}) is below recommended "
                f"minimum (n={min_size_per_season}). Results may be unreliable.",
                UserWarning
            )

    s, var_s, denom = 0, 0, 0
    all_slopes = []
    tau_weighted_sum = 0
    denom_sum = 0


    for i in season_range:
        season_mask = data_filtered['season'] == i
        season_data = data_filtered[season_mask]
        season_x = season_data['value'].to_numpy()
        season_t_raw = season_data['t'].to_numpy()
        season_censored = season_data['censored'].to_numpy()
        season_cen_type = season_data['cen_type'].to_numpy()
        n = len(season_x)

        if n > 1:
            # DEVELOPER NOTE:
            # The 'absolute' time_method (default) uses true numeric timestamps,
            # which is statistically robust for unequally spaced data.
            # The 'rank' method uses cycle-based ranks (1, 2, 3,...) for time,
            # matching the LWP-TRENDS R script's methodology.
            if time_method == 'rank':
                season_cycles = season_data['cycle'].to_numpy()
                season_t = _get_time_ranks(season_t_raw, season_cycles)
            else: # 'absolute'
                season_t = season_t_raw

            s_season, var_s_season, d_season, tau_season = _mk_score_and_var_censored(
                season_x, season_t, season_censored, season_cen_type,
                tau_method=tau_method
            )
            s += s_season
            var_s += var_s_season

            # For weighted average of Tau
            if d_season > 0:
                tau_weighted_sum += tau_season * d_season
                denom_sum += d_season


            if np.any(season_censored):
                all_slopes.extend(_sens_estimator_censored(
                    season_x, season_t, season_cen_type,
                    lt_mult=lt_mult, gt_mult=gt_mult, method=sens_slope_method
                ))
            else:
                all_slopes.extend(_sens_estimator_unequal_spacing(season_x, season_t))

    Tau = tau_weighted_sum / denom_sum if denom_sum > 0 else 0
    z = _z_score(s, var_s)
    p, h, trend = _p_value(z, alpha)
    C, Cd = _mk_probability(p, s)

    if not all_slopes:
        slope, intercept, lower_ci, upper_ci = np.nan, np.nan, np.nan, np.nan
    else:
        slope = np.nanmedian(np.asarray(all_slopes))
        intercept = np.nanmedian(data_filtered['value']) - np.nanmedian(data_filtered['t']) * slope
        lower_ci, upper_ci = _confidence_intervals(np.asarray(all_slopes), var_s, alpha)

    results = res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)

    if plot_path:
        plot_trend(data_filtered, results, plot_path, alpha)

    return results
