"""
This script provides a function to test for seasonality in a time series
using the Kruskal-Wallis H-test.
"""
from collections import namedtuple
import numpy as np
import pandas as pd
from scipy.stats import kruskal
from ._utils import __preprocessing, _get_season_func, _is_datetime_like

def seasonality_test(x_old, t_old, period=12, alpha=0.05, season_type='month'):
    """
    Performs a Kruskal-Wallis H-test to determine if there is a statistically
    significant difference between the distributions of seasons in a time series.

    Input:
        x_old: a vector of data
        t_old: a vector of timestamps
        period: seasonal cycle (default 12)
        alpha: significance level (default 0.05)
        season_type: For datetime inputs, specifies the type of seasonality.
                     'month', 'day_of_week', etc.
    Output:
        A named tuple with the following fields:
        - h_statistic: The Kruskal-Wallis H-statistic.
        - p_value: The p-value of the test.
        - is_seasonal: A boolean indicating if seasonality was detected.
    """
    res = namedtuple('Seasonality_Test', ['h_statistic', 'p_value', 'is_seasonal'])

    x_raw = np.asarray(x_old)
    t_raw = np.asarray(t_old)

    is_datetime = _is_datetime_like(t_raw)

    if is_datetime:
        season_func = _get_season_func(season_type, period)

    mask = ~np.isnan(x_raw)
    x, t = x_raw[mask], t_raw[mask]

    if len(x) < 2:
        return res(np.nan, np.nan, False)

    if is_datetime:
        seasons = season_func(pd.to_datetime(t))
    else:
        t_numeric, _ = __preprocessing(t)
        seasons = (np.floor(t_numeric - 1) % period).astype(int)

    # Kruskal-Wallis H-test requires at least two groups
    unique_seasons = np.unique(seasons)
    if len(unique_seasons) < 2:
        return res(np.nan, np.nan, False)

    seasonal_data = [x[seasons == s] for s in unique_seasons]

    # The test requires at least 5 observations in each group
    if any(len(group) < 5 for group in seasonal_data):
        return res(np.nan, np.nan, False)

    h_statistic, p_value = kruskal(*seasonal_data)

    is_seasonal = p_value < alpha

    return res(h_statistic, p_value, is_seasonal)
