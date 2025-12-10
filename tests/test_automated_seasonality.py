
import numpy as np
import pandas as pd
import pytest
from MannKenSen import automated_seasonal_test

# Test case 1: Data with clear monthly seasonality
def test_with_monthly_seasonality():
    # Generate synthetic data with a clear monthly pattern and a linear trend
    dates = pd.to_datetime(pd.date_range(start='2000-01-01', periods=120, freq='ME'))
    # Create a seasonal cycle
    seasonal_cycle = np.sin(np.arange(120) * (2 * np.pi / 12)) * 5
    # Create a linear trend
    trend = np.linspace(0, 5, 120)
    # Combine them to get the final data
    data = 10 + seasonal_cycle + trend + np.random.normal(0, 1, 120)

    # Run the automated test
    result = automated_seasonal_test(data, dates)

    # Check that 'month' seasonality was detected
    assert result.seasonality == 'month'
    # Check that the trend is identified as increasing
    assert result.trend == 'increasing'

# Test case 2: Data with no seasonality
def test_with_no_seasonality():
    # Generate synthetic data with only a linear trend
    dates = pd.to_datetime(pd.date_range(start='2000-01-01', periods=120, freq='ME'))
    trend = np.linspace(0, 5, 120)
    data = 10 + trend + np.random.normal(0, 1, 120)

    # Run the automated test
    result = automated_seasonal_test(data, dates)

    # Check that 'none' seasonality was detected
    assert result.seasonality == 'none'
    # Check that the trend is identified as increasing
    assert result.trend == 'increasing'

# Test case 3: Non-datetime data
def test_with_numeric_data():
    # Generate simple numeric data with a linear trend
    t = np.arange(100)
    data = 0.1 * t + np.random.normal(0, 1, 100)

    # Run the automated test
    result = automated_seasonal_test(data, t)

    # Check that seasonality is 'none (numeric data)'
    assert result.seasonality == 'none (numeric data)'
    # Check that the trend is identified as increasing
    assert result.trend == 'increasing'
