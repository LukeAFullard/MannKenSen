"""
This script provides a function to automatically detect seasonality and
perform the appropriate Mann-Kendall trend test.
"""
from collections import namedtuple
import numpy as np

from .seasonality_test import seasonality_test
from .seasonal_test import seasonal_test
from .original_test import original_test
from ._utils import _is_datetime_like


def automated_seasonal_test(x_old, t_old, alpha=0.05, seasonality_alpha=0.05):
    """
    Automatically detects the presence and type of seasonality in the data and
    performs the most appropriate Mann-Kendall trend test.

    It first checks if the time series is datetime-like. If not, it falls
    back to the non-seasonal `original_test`.

    If the series is datetime-like, it tests for several types of
    seasonality (monthly, bimonthly, quarterly, biannual) using the
    Kruskal-Wallis test.

    - If significant seasonality is detected (p < seasonality_alpha), it selects the
      seasonality with the strongest effect (highest Kruskal-Wallis statistic)
      and runs the `seasonal_test`.
    - If no significant seasonality is found, it runs the non-seasonal
      `original_test`.

    Input:
        x_old: A vector of data.
        t_old: A vector of timestamps corresponding to x_old. Must be
               datetime-like for seasonality detection.
        alpha: Significance level for the trend test (default 0.05).
        seasonality_alpha: Significance level for the seasonality test (default 0.05).

    Output:
        A named tuple containing the test results, including the type of
        seasonality detected ('none' if no seasonality was found or data
        was not datetime-like).
    """
    AutomatedMannKendallTest = namedtuple(
        'AutomatedMannKendallTest',
        [
            'seasonality', 'trend', 'h', 'p', 'z', 'Tau', 's', 'var_s',
            'slope', 'intercept', 'lower_ci', 'upper_ci', 'C', 'Cd'
        ]
    )

    t_raw = np.asarray(t_old)
    if not _is_datetime_like(t_raw):
        # Fallback to original test if data is not datetime-like
        res = original_test(x_old, t_old, alpha=alpha)
        return AutomatedMannKendallTest('none (numeric data)', *res)

    seasonalities_to_check = {
        'month': 12,
        'bimonth': 6,
        'quarter': 4,
        'biannual': 2,
    }

    results = []
    for season_type, period in seasonalities_to_check.items():
        try:
            test_res = seasonality_test(x_old, t_old, period=period, season_type=season_type)
            if test_res.p_value is not None and not np.isnan(test_res.p_value):
                results.append({
                    'season_type': season_type,
                    'period': period,
                    'stat': test_res.h_statistic,
                    'p_value': test_res.p_value,
                })
        except (ValueError, IndexError):
            # Not enough data for this seasonal aggregation
            continue

    significant_seasons = [res for res in results if res['p_value'] < seasonality_alpha]

    if not significant_seasons:
        res = original_test(x_old, t_old, alpha=alpha)
        return AutomatedMannKendallTest('none', *res)

    best_season = max(significant_seasons, key=lambda x: x['stat'])
    res = seasonal_test(
        x_old,
        t_old,
        period=best_season['period'],
        alpha=alpha,
        season_type=best_season['season_type']
    )
    return AutomatedMannKendallTest(best_season['season_type'], *res)
