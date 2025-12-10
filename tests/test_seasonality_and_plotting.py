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
