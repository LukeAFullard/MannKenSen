"""
This script provides a modified version of the Seasonal Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""
from collections import namedtuple
import numpy as np
import pandas as pd
from ._utils import (__preprocessing, __missing_values_analysis, __mk_score,
                   __variance_s, __z_score, __p_value,
                   __sens_estimator_unequal_spacing, __confidence_intervals,
                   __mk_probability, _get_season_func, _is_datetime_like)

def seasonal_test(x_old, t_old, period=12, alpha=0.05, agg_method='none', season_type='month'):
    """
    Seasonal Mann-Kendall test for unequally spaced time series.
    Input:
        x_old: a vector of data
        t_old: a vector of timestamps
        period: seasonal cycle (default 12)
        alpha: significance level (default 0.05)
        agg_method: method for aggregating multiple data points within a season-year.
                    'none' (default): performs analysis on all data points.
                    'median': uses the median of values and times for each season-year.
                    'middle': uses the observation closest to the middle of the time period.
        season_type: For datetime inputs, specifies the type of seasonality.
                     'year', 'month', 'day_of_week', 'quarter', 'hour', 'week_of_year',
                     'day_of_year', 'minute', 'second'.
    Output:
        trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd
    """
    res = namedtuple('Seasonal_Mann_Kendall_Test', ['trend', 'h', 'p', 'z', 'Tau', 's', 'var_s', 'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'])

    x_raw = np.asarray(x_old)
    t_raw = np.asarray(t_old)

    is_datetime = _is_datetime_like(t_raw)

    if is_datetime:
        season_func = _get_season_func(season_type, period)

    mask = ~np.isnan(x_raw)
    x, t = x_raw[mask], t_raw[mask]

    if len(x) < 2:
        return res('no trend', False, np.nan, 0, 0, 0, 0, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

    # --- Aggregation Logic ---
    if agg_method != 'none':
        if is_datetime:
            t_pd = pd.to_datetime(t)
            years_agg = t_pd.year
            seasons_agg = season_func(t_pd) if season_type != 'year' else np.ones(len(t_pd))
        else:
            t_numeric_agg, _ = __preprocessing(t)
            years_agg = np.floor((t_numeric_agg - 1) / period)
            seasons_agg = np.floor(t_numeric_agg - 1) % period

        unique_year_seasons = np.unique(np.column_stack((years_agg, seasons_agg)), axis=0)
        agg_x, agg_t = [], []

        for year, season in unique_year_seasons:
            group_mask = (years_agg == year) & (seasons_agg == season)
            group_x, group_t = x[group_mask], t[group_mask]

            if len(group_x) > 1:
                if agg_method == 'median':
                    agg_x.append(np.median(group_x))
                    agg_t.append(pd.to_datetime(group_t).to_series().median() if is_datetime else np.median(group_t))
                elif agg_method == 'middle':
                    t_numeric_group, _ = __preprocessing(group_t)
                    closest_idx = np.argmin(np.abs(t_numeric_group - np.mean(t_numeric_group)))
                    agg_x.append(group_x[closest_idx])
                    agg_t.append(group_t[closest_idx])
            else:
                agg_x.append(group_x[0])
                agg_t.append(group_t[0])

        x, t = np.array(agg_x), np.array(agg_t)

    # --- Trend Analysis ---
    t_numeric, _ = __preprocessing(t)
    if is_datetime and season_type != 'year':
        seasons = season_func(pd.to_datetime(t))
        season_range = np.unique(seasons)
    elif not is_datetime:
        seasons = (np.floor(t_numeric - 1) % period).astype(int)
        season_range = range(int(period))
    else: # is_datetime and season_type == 'year'
        seasons = np.ones(len(t))
        season_range = [1]


    s, var_s, denom = 0, 0, 0
    all_slopes = []

    for i in season_range:
        season_mask = seasons == i
        season_x, season_t = x[season_mask], t_numeric[season_mask]
        n = len(season_x)

        if n > 1:
            s += __mk_score(season_x, n)
            var_s += __variance_s(season_x, n)
            denom += 0.5 * n * (n - 1)
            all_slopes.extend(__sens_estimator_unequal_spacing(season_x, season_t))

    Tau = s / denom if denom != 0 else 0
    z = __z_score(s, var_s)
    p, h, trend = __p_value(z, alpha)
    C, Cd = __mk_probability(p, s)

    if not all_slopes:
        slope, intercept, lower_ci, upper_ci = np.nan, np.nan, np.nan, np.nan
    else:
        slope = np.nanmedian(np.asarray(all_slopes))
        intercept = np.nanmedian(x) - np.nanmedian(t_numeric) * slope
        lower_ci, upper_ci = __confidence_intervals(np.asarray(all_slopes), var_s, alpha)

    return res(trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci, C, Cd)
