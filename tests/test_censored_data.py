import pytest
import numpy as np
import pandas as pd
from MannKenSen import original_test, seasonal_test
from MannKenSen._utils import prepare_censored_data

# rpy2 setup to call the original R script for ground truth
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
pandas2ri.activate()
ro.r.source('Example_Files/R/LWPTrends_v2502/LWPTrends_v2502.r')

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

# Integration tests using rpy2 to compare with the original R script
CENSORED_DATA = {
    't': np.arange(2000, 2011),
    'x': [6, 5, '<4', 3, 4, 2, '<2', 5, 2, 1, '<4']
}

SEASONAL_CENSORED_DATA = {
    't': pd.to_datetime(pd.date_range(start='2000-01-01', periods=24, freq='ME')),
    'x': [10, '<8', 12, 11, 15, '>18', 9, 7, 11, 10, 14, 17, 8, '<6', 10, 9, 13, 16, 7, 5, 9, 8, 12, 15]
}


def run_r_test(t, x, hicensor=False):
    """Helper function to run the R script and get results."""
    # Convert data to R DataFrame
    r_df = ro.DataFrame({'myDate': ro.vectors.DateVector([str(d) for d in t]),
                         'Value': ro.StrVector(x)})
    # Pre-process in R
    r_df = ro.r.RemoveAlphaDetect(r_df)
    r_df = ro.r.GetMoreDateInfo(r_df)
    r_df['TimeIncr'] = r_df.rx2('Year') # For non-seasonal test

    # Run the non-seasonal test in R
    if hicensor:
        # The R script's HiCensor is a bit different, we simulate it
        r_df_hicensor = ro.r.ApplyHiCensor(r_df, HiCensor=True)
        result_r = ro.r.MannKendall(r_df_hicensor)
    else:
        result_r = ro.r.MannKendall(r_df)

    # Convert results to a pandas Series for easy comparison
    with localconverter(ro.default_converter + pandas2ri.converter):
        result_pd = ro.conversion.rpy2py(result_r)
    return result_pd.iloc[0]

@pytest.fixture
def python_results():
    """Fixture to run the Python implementation."""
    x_prepared = prepare_censored_data(CENSORED_DATA['x'])
    return original_test(x_prepared, CENSORED_DATA['t'])

@pytest.fixture
def r_results():
    """Fixture to run the R implementation."""
    t_dates = [f"{year}-01-01" for year in CENSORED_DATA['t']]
    return run_r_test(t_dates, [str(v) for v in CENSORED_DATA['x']])

def test_original_test_censored(python_results, r_results):
    """Compare Python and R results for the non-seasonal test."""
    assert python_results.s == pytest.approx(r_results['S'], abs=1e-9)
    assert python_results.var_s == pytest.approx(r_results['VarS'], abs=1e-9)
    assert python_results.p == pytest.approx(r_results['p'], abs=1e-3)

@pytest.fixture
def python_results_hicensor():
    """Fixture to run the Python implementation with HiCensor."""
    x_prepared = prepare_censored_data(CENSORED_DATA['x'])
    return original_test(x_prepared, CENSORED_DATA['t'], hicensor=True)

@pytest.fixture
def r_results_hicensor():
    """Fixture to run the R implementation with HiCensor."""
    t_dates = [f"{year}-01-01" for year in CENSORED_DATA['t']]
    return run_r_test(t_dates, [str(v) for v in CENSORED_DATA['x']], hicensor=True)

def test_original_test_hicensor(python_results_hicensor, r_results_hicensor):
    """Compare Python and R results for the non-seasonal test with HiCensor."""
    assert python_results_hicensor.s == pytest.approx(r_results_hicensor['S'], abs=1e-9)
    assert python_results_hicensor.var_s == pytest.approx(r_results_hicensor['VarS'], abs=1e-9)
    assert python_results_hicensor.p == pytest.approx(r_results_hicensor['p'], abs=1e-3)

def run_r_seasonal_test(t, x, period=12):
    """Helper function to run the R seasonal test and get results."""
    # Convert data to R DataFrame
    r_df = ro.DataFrame({'myDate': ro.vectors.DateVector(t.strftime('%Y-%m-%d')),
                         'Value': ro.StrVector([str(v) for v in x])})

    # Pre-process in R
    r_df = ro.r.RemoveAlphaDetect(r_df)
    r_df = ro.r.GetMoreDateInfo(r_df)
    r_df['TimeIncr'] = r_df.rx2('Month')
    r_df['Season'] = r_df.rx2('Month')

    # Run the seasonal test in R
    result_r = ro.r.SeasonalKendall(r_df)

    # Convert results to a pandas Series
    with localconverter(ro.default_converter + pandas2ri.converter):
        result_pd = ro.conversion.rpy2py(result_r)
    return result_pd.iloc[0]

@pytest.fixture
def python_seasonal_results():
    """Fixture to run the Python seasonal implementation."""
    x_prepared = prepare_censored_data(SEASONAL_CENSORED_DATA['x'])
    return seasonal_test(x_prepared, SEASONAL_CENSORED_DATA['t'], period=12)

@pytest.fixture
def r_seasonal_results():
    """Fixture to run the R seasonal implementation."""
    return run_r_seasonal_test(SEASONAL_CENSORED_DATA['t'], SEASONAL_CENSORED_DATA['x'])

def test_seasonal_test_censored(python_seasonal_results, r_seasonal_results):
    """Compare Python and R results for the seasonal test."""
    assert python_seasonal_results.s == pytest.approx(r_seasonal_results['S'], abs=1e-9)
    assert python_seasonal_results.var_s == pytest.approx(r_seasonal_results['VarS'], abs=1e-9)
    assert python_seasonal_results.p == pytest.approx(r_seasonal_results['p'], abs=1e-3)
