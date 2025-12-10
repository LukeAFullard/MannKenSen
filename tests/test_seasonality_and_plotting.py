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
