
import numpy as np
import pandas as pd
import pytest
from MannKenSen.analysis_notes import get_analysis_note, get_sens_slope_analysis_note
from MannKenSen._utils import _sens_estimator_censored

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
        'value': [1, 1, 1, 1, 2, 6, 7, 8],
        'censored': [False, False, False, False, False, False, False, False],
        'season': [1, 1, 1, 1, 1, 2, 2, 2]
    })
    note = get_analysis_note(data, is_seasonal=True, post_aggregation=True)
    assert note == "Long run of single value in a Season"

def test_get_sens_slope_analysis_note_ok():
    slopes = np.array([1, 2, 3, 4, 5])
    t = np.array([1, 2, 3, 4, 5])
    cen_type = np.array(['not', 'not', 'not', 'not', 'not'])
    note = get_sens_slope_analysis_note(slopes, t, cen_type)
    assert note == "ok"

def test_get_sens_slope_analysis_note_influenced_by_censored():
    x = np.array([10, 2, 3, 4, 5])
    t = np.arange(5)
    cen_type = np.array(['not', 'lt', 'not', 'not', 'not'])
    slopes = _sens_estimator_censored(x, t, cen_type)
    note = get_sens_slope_analysis_note(slopes, t, cen_type)
    assert note == "WARNING: Sen slope influenced by censored values"

def test_get_sens_slope_analysis_note_based_on_two_censored():
    x = np.array([1, 2, 3, 4, 5])
    t = np.arange(5)
    cen_type = np.array(['lt', 'lt', 'lt', 'lt', 'lt'])
    slopes = _sens_estimator_censored(x, t, cen_type)
    note = get_sens_slope_analysis_note(slopes, t, cen_type)
    assert note == "WARNING: Sen slope based on two censored values"

def test_get_sens_slope_analysis_note_tied_non_censored():
    x = np.array([1, 2, 1, 2, 1])
    t = np.arange(5)
    cen_type = np.array(['not', 'not', 'not', 'not', 'not'])
    slopes = _sens_estimator_censored(x, t, cen_type)
    note = get_sens_slope_analysis_note(slopes, t, cen_type)
    assert note == "WARNING: Sen slope based on tied non-censored values"
