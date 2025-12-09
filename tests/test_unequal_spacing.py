import pytest
import numpy as np
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
    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = original_test(x, t)

    assert trend == 'increasing'
    assert h == True
    assert slope == pytest.approx(2.0, abs=0.5)
    assert lower_ci <= slope <= upper_ci

def test_seasonal_test_unequal_spacing(unequal_seasonal_data):
    t, x = unequal_seasonal_data
    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = seasonal_test(x, t, period=12)

    assert trend == 'no trend'
    assert lower_ci <= slope <= upper_ci

def test_original_test_no_trend():
    np.random.seed(0)
    t = np.linspace(0, 10, 20)
    x = np.full_like(t, 5.0)
    x += np.random.normal(0, 0.1, size=x.shape)

    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = original_test(x, t)

    assert trend == 'no trend'
    assert h == False
    assert slope == pytest.approx(0.0, abs=0.1)
    assert lower_ci <= slope <= upper_ci

def test_seasonal_test_with_trend():
    t = np.arange(0, 48, 1)
    seasonal_component = np.sin(2 * np.pi * t / 12)
    trend_component = 0.1 * t
    x = seasonal_component + trend_component

    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = seasonal_test(x, t, period=12)

    assert trend == 'increasing'
    assert h == True
    assert slope == pytest.approx(0.1, abs=0.05)
    assert lower_ci <= slope <= upper_ci

def test_original_test_with_nan():
    t = np.array([1, 2, 3, 4, 5, 6])
    x = np.array([1, 2, np.nan, 4, 5, 6])
    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = original_test(x, t)
    assert trend == 'increasing'
    assert h == True
    assert slope == pytest.approx(1.0, abs=0.1)
    assert lower_ci <= slope <= upper_ci

def test_seasonal_test_with_nan():
    t = np.arange(0, 48, 1).astype(float)
    x = 0.1 * t + np.sin(2 * np.pi * t / 12)
    x[5] = np.nan
    x[20] = np.nan

    trend, h, p, z, Tau, s, var_s, slope, intercept, lower_ci, upper_ci = seasonal_test(x, t, period=12)

    assert trend == 'increasing'
    assert h == True
    assert slope == pytest.approx(0.1, abs=0.05)
    assert lower_ci <= slope <= upper_ci
