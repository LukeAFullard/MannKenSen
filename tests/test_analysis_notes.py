
import numpy as np
import pandas as pd
import pytest
import warnings
from MannKenSen.analysis_notes import get_analysis_note, get_sens_slope_analysis_note
from MannKenSen._utils import _sens_estimator_censored
from MannKenSen.preprocessing import prepare_censored_data
from MannKenSen.original_test import original_test
from MannKenSen.seasonal_test import seasonal_test

def test_get_analysis_note_all_na():
    data = pd.DataFrame({'value': [np.nan, np.nan, np.nan], 'censored': [False, False, False]})
    note = get_analysis_note(data)
    assert note == "Data all NA values"

def test_get_analysis_note_less_than_3_unique_values():
    data = pd.DataFrame({'value': [1, 1, 2, 2, 1], 'censored': [False, False, False, False, False]})
    note = get_analysis_note(data)
    assert note == "< 3 unique values"

def test_get_analysis_note_less_than_5_non_censored_values():
    data = pd.DataFrame({'value': [1, 2, 3, 4, 5], 'censored': [True, True, False, False, False]})
    note = get_analysis_note(data)
    assert note == "< 5 Non-censored values"

def test_get_analysis_note_long_run_of_single_value():
    data = pd.DataFrame({'value': [1]*11 + [2, 3, 4, 5, 6, 7, 8, 9, 10], 'censored': [False]*20})
    note = get_analysis_note(data, post_aggregation=True)
    assert note == "Long run of single value"

def test_get_analysis_note_seasonal_less_than_3_non_na():
    data = pd.DataFrame({
        'value': [1, 2, np.nan, 4, 5, 6],
        'censored': [False, False, False, False, False, False],
        'season': [1, 1, 1, 2, 2, 2]
    })
    note = get_analysis_note(data, is_seasonal=True, post_aggregation=True)
    assert note == "< 3 non-NA values in Season"


def test_get_analysis_note_seasonal_less_than_2_unique_values():
    data = pd.DataFrame({
        'value': [1, 1, 1, 4, 5, 6],
        'censored': [False, False, False, False, False, False],
        'season': [1, 1, 1, 2, 2, 2]
    })
    note = get_analysis_note(data, is_seasonal=True, post_aggregation=True)
    assert note == "< 2 unique values in Season"

def test_get_analysis_note_seasonal_long_run():
    data = pd.DataFrame({
        'value': [1, 1, 1, 1, 2, 5, 6, 7, 8],
        'censored': [False, False, False, False, False, False, False, False, False],
        'season': [1, 1, 1, 1, 1, 2, 2, 2, 2]
    })
    note = get_analysis_note(data, is_seasonal=True, post_aggregation=True)
    assert note == "Long run of single value in a Season"

def test_get_sens_slope_analysis_note_ok():
    slopes = np.array([1, 2, 3, 4, 5])
    t = np.array([1, 2, 3, 4, 5])
    cen_type = np.array(['not', 'not', 'not', 'not', 'not'])
    note = get_sens_slope_analysis_note(slopes, t, cen_type)
    assert note == "ok"

def test_get_sens_slope_analysis_note_influenced_by_left_censored():
    x_raw = ['<1', 5, 10, 12]
    t = np.arange(len(x_raw))
    data = prepare_censored_data(x_raw)
    slopes = _sens_estimator_censored(data['value'].values, t, data['cen_type'].values, method='lwp')
    note = get_sens_slope_analysis_note(slopes, t, data['cen_type'].values)
    assert note == "WARNING: Sen slope influenced by left-censored values"

def test_get_sens_slope_analysis_note_influenced_by_right_censored():
    x_raw = [1, 5, 2, '>10']
    t = np.arange(len(x_raw))
    data = prepare_censored_data(x_raw)
    slopes = _sens_estimator_censored(data['value'].values, t, data['cen_type'].values, method='lwp')
    note = get_sens_slope_analysis_note(slopes, t, data['cen_type'].values)
    assert note == "WARNING: Sen slope influenced by right-censored values"

def test_get_sens_slope_analysis_note_influenced_by_both_censored():
    x_raw = ['<1', 5, 2, '>10']
    t = np.arange(len(x_raw))
    data = prepare_censored_data(x_raw)
    slopes = _sens_estimator_censored(data['value'].values, t, data['cen_type'].values, method='lwp')
    note = get_sens_slope_analysis_note(slopes, t, data['cen_type'].values)
    assert note == "WARNING: Sen slope influenced by left- and right-censored values"

def test_get_sens_slope_analysis_note_based_on_two_censored():
    x_raw = ['<1', '<2', '>3', '>4']
    t = np.arange(len(x_raw))
    data = prepare_censored_data(x_raw)
    slopes = _sens_estimator_censored(data['value'].values, t, data['cen_type'].values, method='lwp')
    note = get_sens_slope_analysis_note(slopes, t, data['cen_type'].values)
    assert note == "WARNING: Sen slope based on two censored values"

def test_get_sens_slope_analysis_note_tied_non_censored():
    x_raw = [1, 2, 1, 2, 1]
    t = np.arange(len(x_raw))
    data = prepare_censored_data(x_raw)
    slopes = _sens_estimator_censored(data['value'].values, t, data['cen_type'].values, method='lwp')
    note = get_sens_slope_analysis_note(slopes, t, data['cen_type'].values)
    assert note == "WARNING: Sen slope based on tied non-censored values"

def test_original_test_min_size_warning():
    """Test that original_test issues a warning for small sample sizes."""
    x = [1, 2, 3]
    t = [1, 2, 3]
    with pytest.warns(UserWarning, match="Sample size .* is below recommended minimum"):
        original_test(x, t, min_size=4)

    # Test that no warning is issued when min_size is None
    with warnings.catch_warnings(record=True) as record:
        original_test(x, t, min_size=None)
        assert len(record) == 0

def test_seasonal_test_min_size_warning():
    """Test that seasonal_test issues a warning for small seasonal sample sizes."""
    # Season 1 has 3 values, Season 2 has 2
    x = [1, 2, 3, 4, 5]
    t = [1, 2, 3, 13, 14] # Using numeric time for simplicity
    with pytest.warns(UserWarning, match="Minimum season size .* is below recommended minimum"):
        seasonal_test(x, t, period=12, min_size_per_season=3)

    # Test that no warning is issued when min_size_per_season is None
    with warnings.catch_warnings(record=True) as record:
        seasonal_test(x, t, period=12, min_size_per_season=None)
        assert len(record) == 0
