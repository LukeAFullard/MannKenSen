import pytest
import numpy as np
import pandas as pd
from MannKenSen import original_test, seasonal_test
from MannKenSen._utils import prepare_censored_data

# Unit tests for the new `prepare_censored_data` function
def test_prepare_censored_data_valid():
    x = [1, '<2', 3, '>4', '5']
    expected_df = pd.DataFrame({
        'value': [1.0, 2.0, 3.0, 4.0, 5.0],
        'censored': [False, True, False, True, False],
        'cen_type': ['not', 'lt', 'not', 'gt', 'not']
    })
    result_df = prepare_censored_data(x)
    pd.testing.assert_frame_equal(result_df, expected_df)

def test_prepare_censored_data_all_numeric():
    x = [1, 2, 3, 4, 5]
    expected_df = pd.DataFrame({
        'value': [1.0, 2.0, 3.0, 4.0, 5.0],
        'censored': [False, False, False, False, False],
        'cen_type': ['not', 'not', 'not', 'not', 'not']
    })
    result_df = prepare_censored_data(x)
    pd.testing.assert_frame_equal(result_df, expected_df)

def test_prepare_censored_data_with_spaces():
    x = [' < 2 ', ' > 4 ']
    expected_df = pd.DataFrame({
        'value': [2.0, 4.0],
        'censored': [True, True],
        'cen_type': ['lt', 'gt']
    })
    result_df = prepare_censored_data(x)
    pd.testing.assert_frame_equal(result_df, expected_df)

def test_prepare_censored_data_invalid_string():
    with pytest.raises(ValueError, match="Could not convert string 'abc' to a float."):
        prepare_censored_data(['abc'])

def test_prepare_censored_data_invalid_censored_format():
    with pytest.raises(ValueError, match="Invalid left-censored value format: '<'"):
        prepare_censored_data(['<'])
    with pytest.raises(ValueError, match="Invalid right-censored value format: '>a'"):
        prepare_censored_data(['>a'])

def test_prepare_censored_data_non_iterable():
    with pytest.raises(TypeError):
        prepare_censored_data(123)
