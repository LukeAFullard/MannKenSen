import pytest
import numpy as np
import pandas as pd
from MannKenSen.original_test import original_test
from MannKenSen.seasonal_test import seasonal_test

@pytest.fixture
def unequal_linear_data():
    t = np.array([1, 2, 4, 7, 8, 10])
    x = 2 * t + 5 + np.random.rand(len(t))
    return t, x

@pytest.fixture
def unequal_seasonal_data():
    t = np.array([1, 2, 3, 13, 14, 15, 25, 26, 27])
    x = 10 * np.sin(2 * np.pi * t / 12) + np.random.rand(len(t))
    return t, x

def test_original_test_unequal_spacing(unequal_linear_data):
    t, x = unequal_linear_data
    result = original_test(x, t)

    assert result.trend == 'increasing'
    assert result.h
    assert result.slope == pytest.approx(2.0, abs=0.5)
    assert result.lower_ci <= result.slope <= result.upper_ci
    assert result.C > 0.95
    assert result.Cd < 0.05


def test_seasonal_test_day_of_week_seasonality():
    # 10 weeks of daily data
    t = pd.to_datetime(pd.date_range(start='2023-01-01', periods=70, freq='D'))
    # Constant value on non-trend days
    x = np.full(70, 5.0)
    # Clear, strong increasing trend on Mondays (dayofweek=0)
    mondays = t.dayofweek == 0
    x[mondays] = np.arange(10, 110, 10)  # [10, 20, ..., 100]

    result = seasonal_test(x, t, period=7, season_type='day_of_week')

    assert result.trend == 'increasing'
    assert result.h
    assert result.s == 45  # s for n=10 is 45, s for other constant days is 0


def test_seasonal_test_parameter_validation():
    t = pd.to_datetime(pd.date_range(start='2023-01-01', periods=10, freq='D'))
    x = np.arange(10)

    # Test that an incorrect period for a given season_type raises an error
    with pytest.raises(ValueError, match="For season_type='day_of_week', period must be 7."):
        seasonal_test(x, t, period=12, season_type='day_of_week')

    # Test that an unknown season_type raises an error
    with pytest.raises(ValueError, match="Unknown season_type: 'invalid_season'."):
        seasonal_test(x, t, period=7, season_type='invalid_season')


def test_seasonal_test_datetime_aggregation():
    # 5 years of data with two observations in Jan
    t = np.array([
        '2023-01-10', '2023-01-20',
        '2024-01-12', '2024-01-22',
        '2025-01-15', '2025-01-25',
        '2026-01-11', '2026-01-21',
        '2027-01-14', '2027-01-24',
    ], dtype='datetime64')
    x = np.array([10, 100, 20, 200, 30, 300, 40, 400, 50, 500])

    result = seasonal_test(x, t, period=12, agg_method='median')

    assert result.trend == 'increasing'
    assert result.h
    assert result.s == 10 # For x_agg = [55, 110, 165, 220, 275]
    assert result.C > 0.95
    assert result.Cd < 0.05


@pytest.fixture
def seasonal_data_with_duplicates():
    # Two observations in the first 'month' of each of 5 'years'
    t = np.array([1.1, 1.9, 13.1, 13.9, 25.1, 25.9, 37.1, 37.9, 49.1, 49.9])
    # Clear increasing trend year-over-year
    x = np.array([10, 100, 20, 200, 30, 300, 40, 400, 50, 500])
    return t, x


def test_seasonal_test_aggregation(seasonal_data_with_duplicates):
    t, x = seasonal_data_with_duplicates

    # Test 'median' aggregation
    result_med = seasonal_test(x, t, period=12, agg_method='median')
    assert result_med.trend == 'increasing'
    assert result_med.s == 10

    # Test 'middle' aggregation
    result_mid = seasonal_test(x, t, period=12, agg_method='middle')
    assert result_mid.trend == 'increasing'
    assert result_mid.s == 10

    # Test 'none' (default) aggregation
    result_none = seasonal_test(x, t, period=12, agg_method='none')
    assert result_none.trend == 'increasing'
    assert result_none.s == 25


def test_original_test_with_datetime64():
    t = np.array(['2023-01-01', '2023-01-02', '2023-01-04', '2023-01-07', '2023-01-08', '2023-01-10'], dtype='datetime64')
    x = 2 * np.arange(len(t)) + 5 + np.random.rand(len(t))
    result = original_test(x, t)

    assert result.trend == 'increasing'
    assert result.h
    assert result.C > 0.95
    assert result.Cd < 0.05


def test_original_test_with_datetime_objects():
    from datetime import datetime
    t = [datetime(2023, 1, 1), datetime(2023, 1, 2), datetime(2023, 1, 4), datetime(2023, 1, 7), datetime(2023, 1, 8), datetime(2023, 1, 10)]
    x = 2 * np.arange(len(t)) + 5 + np.random.rand(len(t))
    result = original_test(x, t)

    assert result.trend == 'increasing'
    assert result.h
    assert result.C > 0.95
    assert result.Cd < 0.05


def test_seasonal_test_unequal_spacing(unequal_seasonal_data):
    t, x = unequal_seasonal_data
    result = seasonal_test(x, t, period=12)
    assert result.trend == 'no trend'
    assert not result.h
    assert result.C < 0.95


def test_original_test_no_trend():
    np.random.seed(0)
    t = np.linspace(0, 10, 20)
    x = np.full_like(t, 5.0)
    x += np.random.normal(0, 0.1, size=x.shape)

    result = original_test(x, t)

    assert result.trend == 'no trend'
    assert not result.h
    assert result.slope == pytest.approx(0.0, abs=0.1)
    assert result.lower_ci <= result.slope <= result.upper_ci


def test_seasonal_test_with_trend():
    t = np.arange(0, 48, 1)
    seasonal_component = np.sin(2 * np.pi * t / 12)
    trend_component = 0.1 * t
    x = seasonal_component + trend_component

    result = seasonal_test(x, t, period=12)

    assert result.trend == 'increasing'
    assert result.h
    assert result.slope == pytest.approx(0.1, abs=0.05)
    assert result.lower_ci <= result.slope <= result.upper_ci
    assert result.C > 0.95
    assert result.Cd < 0.05


def test_original_test_with_nan():
    t = np.array([1, 2, 3, 4, 5, 6])
    x = np.array([1, 2, np.nan, 4, 5, 6])
    result = original_test(x, t)
    assert result.trend == 'increasing'
    assert result.h
    assert result.slope == pytest.approx(1.0, abs=0.1)
    assert result.lower_ci <= result.slope <= result.upper_ci
    assert result.C > 0.95
    assert result.Cd < 0.05


def test_seasonal_test_with_nan():
    t = np.arange(0, 48, 1).astype(float)
    x = 0.1 * t + np.sin(2 * np.pi * t / 12)
    x[5] = np.nan
    x[20] = np.nan

    result = seasonal_test(x, t, period=12)

    assert result.trend == 'increasing'
    assert result.h
    assert result.slope == pytest.approx(0.1, abs=0.05)
    assert result.lower_ci <= result.slope <= result.upper_ci
    assert result.C > 0.95
    assert result.Cd < 0.05
