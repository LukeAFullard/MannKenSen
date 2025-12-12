import pytest
import numpy as np
import pandas as pd
import os
from MannKenSen.seasonality_test import seasonality_test
from MannKenSen.plotting import plot_seasonal_distribution
from MannKenSen import original_test, seasonal_test

@pytest.fixture
def seasonal_data():
    # 5 years of monthly data to ensure at least 5 points per season
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=60, freq='ME'))
    x = 10 * np.sin(2 * np.pi * t.month / 12) + 50 + np.random.rand(60)
    return x, t

@pytest.fixture
def non_seasonal_data():
    # 5 years of monthly data
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=60, freq='ME'))
    x = np.linspace(0, 10, 60) + np.random.rand(60)
    return x, t

def test_seasonality_test_detects_seasonality(seasonal_data):
    x, t = seasonal_data
    result = seasonality_test(x, t)
    assert result.is_seasonal
    assert result.p_value < 0.05

def test_seasonality_test_rejects_non_seasonal(non_seasonal_data):
    x, t = non_seasonal_data
    result = seasonality_test(x, t)
    assert not result.is_seasonal
    assert result.p_value > 0.05

def test_plot_seasonal_distribution(seasonal_data):
    x, t = seasonal_data
    save_path = "test_plot.png"

    # Ensure the file does not exist before the test
    if os.path.exists(save_path):
        os.remove(save_path)

    returned_path = plot_seasonal_distribution(x, t, save_path=save_path)

    assert returned_path == save_path
    assert os.path.exists(save_path)

    # Clean up the created file
    os.remove(save_path)

def test_trend_plotting():
    """
    Tests that the plotting functionality in original_test and
    seasonal_test creates a file.
    """
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=20, freq='YE'))
    x = np.arange(20)

    # Test original_test plotting
    original_plot_path = "original_test_plot.png"
    if os.path.exists(original_plot_path):
        os.remove(original_plot_path)

    original_test(x, t, plot_path=original_plot_path)
    assert os.path.exists(original_plot_path)
    os.remove(original_plot_path)

    # Test seasonal_test plotting
    seasonal_plot_path = "seasonal_test_plot.png"
    if os.path.exists(seasonal_plot_path):
        os.remove(seasonal_plot_path)

    seasonal_test(x, t, plot_path=seasonal_plot_path)
    assert os.path.exists(seasonal_plot_path)
    os.remove(seasonal_plot_path)

# New tests for biweekly seasonality
def test_biweekly_seasonality():
    """
    Tests the biweekly seasonality functionality.
    """
    # Create a dataset with a known biweekly pattern
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=104, freq='W'))
    x = np.array([i % 2 for i in range(104)]) # Alternating values every week, creating a biweekly pattern

    result = seasonal_test(x, t, season_type='biweekly', period=26)
    assert result.trend == 'no trend'

def test_biweekly_seasonality_53_week_year():
    """
    Tests the biweekly seasonality with a 53-week year.
    """
    t = pd.to_datetime(pd.date_range(start='2015-01-01', end='2015-12-31', freq='D'))
    x = np.arange(len(t))

    result = seasonal_test(x, t, season_type='biweekly', period=27)
    assert result.trend == 'increasing'

# Edge case tests based on code review feedback
def test_day_of_year_seasonality_leap_year():
    """
    Tests 'day_of_year' seasonality, ensuring it handles leap years correctly.
    The 'period' for day_of_year is dynamic and should not be specified.
    """
    # Data spans a leap year (2020) and a common year (2021)
    t = pd.to_datetime(pd.date_range(start='2020-01-01', end='2021-12-31', freq='D'))
    x = np.arange(len(t))

    # The function should dynamically determine the number of seasons
    result = seasonal_test(x, t, season_type='day_of_year')
    assert result.trend == 'increasing'

def test_week_of_year_seasonality_53_week_year():
    """
    Tests 'week_of_year' seasonality for a year with 53 weeks.
    """
    # 2015 was a 53-week year
    t = pd.to_datetime(pd.date_range(start='2015-01-01', end='2015-12-31', freq='D'))
    x = np.arange(len(t))

    result = seasonal_test(x, t, season_type='week_of_year', period=53)
    assert result.trend == 'increasing'


# Parameterized test for robust season types
@pytest.mark.parametrize("season_type, period, freq, n_periods", [
    ('year', 1, 'YE', 20),
    ('month', 12, 'ME', 60),
    ('day_of_week', 7, 'D', 365 * 2),
    ('quarter', 4, 'QE', 40),
    ('hour', 24, 'h', 168 * 2),
    ('biweekly', 26, 'W', 104 * 2),
    ('minute', 60, 'min', 1440 * 2),
    ('second', 60, 's', 3600 * 2),
])
def test_general_season_types(season_type, period, freq, n_periods):
    """
    Tests the more straightforward season types in the season_map.
    Edge cases like 'day_of_year' and 'week_of_year' are tested separately.
    """
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=n_periods, freq=freq))
    x = np.arange(len(t))

    result = seasonal_test(x, t, season_type=season_type, period=period)
    assert result.trend == 'increasing'

def test_seasonal_test_aggregation_methods():
    """
    Test the 'median' and 'middle' aggregation methods in seasonal_test.
    """
    # Create a dataset with multiple observations per season-year
    dates = pd.to_datetime(['2020-01-10', '2020-01-20', '2021-01-15', '2022-01-05', '2022-01-25',
                            '2020-02-10', '2020-02-20', '2021-02-15', '2022-02-05', '2022-02-25',
                           '2023-01-10', '2023-01-20', '2024-01-15', '2025-01-05', '2025-01-25'])
    np.random.seed(0)
    values = np.arange(15) + np.random.normal(0, 2, 15) # Clear increasing trend with more noise

    # Test with no aggregation
    result_none = seasonal_test(x=values, t=dates, period=12, season_type='month')
    assert result_none.trend == 'increasing'

    # Test with 'median' aggregation
    # The data for Jan 2020 will be aggregated to a single point with value 1.5
    # The data for Jan 2022 will be aggregated to a single point with value 4.5
    result_median = seasonal_test(x=values, t=dates, period=12, season_type='month', agg_method='median')
    assert result_median.trend == 'increasing'

    # Test with 'middle' aggregation
    # The data for Jan 2020 will be aggregated to the point on 2020-01-20 (value 2)
    # The data for Jan 2022 will be aggregated to the point on 2022-01-25 (value 5)
    result_middle = seasonal_test(x=values, t=dates, period=12, season_type='month', agg_method='middle')
    assert result_middle.trend == 'increasing'

    # The scores and slopes should be different for each aggregation method
    assert result_none.s != result_median.s
    assert result_median.s != result_middle.s
    assert result_none.slope != result_median.slope
    assert result_median.slope != result_middle.slope

def test_seasonality_test_insufficient_data():
    """Test seasonality_test with insufficient data."""
    x = [1, 2, 3]
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=3, freq='ME'))
    result = seasonality_test(x, t)
    assert np.isnan(result.h_statistic)
    assert np.isnan(result.p_value)
    assert not result.is_seasonal

def test_plot_seasonal_distribution_insufficient_data():
    """Test plot_seasonal_distribution with insufficient data."""
    x = [1]
    t = [pd.to_datetime('2020-01-01')]
    result = plot_seasonal_distribution(x, t, save_path='test.png')
    assert result is None

def test_seasonality_test_insufficient_unique_values():
    """
    Test seasonality_test returns no trend if a season has enough points
    but not enough unique values.
    """
    # Create a dataset where one season has 5 identical points
    t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=60, freq='ME'))
    # Convert to numpy array to make it mutable
    x = np.array(10 * np.sin(2 * np.pi * t.month / 12) + 50 + np.random.rand(60))

    # Corrupt the data for January (month=1) to have only one unique value
    x[t.month == 1] = 42

    result = seasonality_test(x, t)
    assert not result.is_seasonal
    assert np.isnan(result.h_statistic)
    assert np.isnan(result.p_value)

def test_seasonal_time_method():
    """
    Validates that the `time_method` parameter in `seasonal_test` produces
    different results for 'absolute' and 'rank' methods.
    """
    # Create a dataset with a known trend and some irregularity in spacing
    # to highlight the difference between the time methods.
    dates = pd.to_datetime([
        '2020-01-10', '2020-02-15', '2020-03-01', # Early in the month
        '2021-01-25', '2021-02-20', '2021-03-28', # Late in the month
        '2022-01-15', '2022-02-15', '2022-03-15', # Middle of the month
        '2023-01-05', '2023-02-10', '2023-03-20',
    ])
    values = np.arange(len(dates)) + np.random.normal(0, 1, len(dates))

    # --- Absolute method (default) ---
    # This method should consider the precise timestamps.
    result_abs = seasonal_test(
        x=values,
        t=dates,
        season_type='month',
        time_method='absolute'
    )
    assert result_abs.trend == 'increasing'

    # --- Rank method (LWP-TRENDS style) ---
    # This method should rank by year, ignoring the intra-month spacing.
    result_rank = seasonal_test(
        x=values,
        t=dates,
        season_type='month',
        time_method='rank'
    )
    assert result_rank.trend == 'increasing'

    # --- Verification ---
    # The statistical results (slope and variance) should be different
    # because the underlying time vectors are treated differently. The S-statistic
    # is expected to be the same for this test case because the chronological
    # order of the data points within each season does not change.
    assert result_abs.slope != result_rank.slope, "Sen's slope should differ between time methods."
    assert result_abs.var_s == result_rank.var_s, "Variance should be the same because it is rank-based."
    assert result_abs.s == result_rank.s, "S-statistic should be the same for this test case."

    # Test that an invalid method raises an error
    with pytest.raises(ValueError):
        seasonal_test(x=values, t=dates, season_type='month', time_method='invalid_method')
