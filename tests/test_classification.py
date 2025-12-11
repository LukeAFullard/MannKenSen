import pytest
from collections import namedtuple
from MannKenSen.classification import classify_trend

# Define a mock result object for testing
TrendResult = namedtuple('TrendResult', ['h', 'C', 'trend'])

def test_classify_no_trend():
    """Test that a non-significant result is classified as 'No Trend'."""
    result = TrendResult(h=False, C=0.99, trend='increasing')
    assert classify_trend(result) == "No Trend"

def test_classify_default_categories():
    """Test the default IPCC-style classification."""
    # Highly Likely Increasing
    result_hl = TrendResult(h=True, C=0.96, trend='increasing')
    assert classify_trend(result_hl) == "Highly Likely Increasing"

    # Very Likely Decreasing
    result_vl = TrendResult(h=True, C=0.92, trend='decreasing')
    assert classify_trend(result_vl) == "Very Likely Decreasing"

    # Likely Increasing
    result_l = TrendResult(h=True, C=0.75, trend='increasing')
    assert classify_trend(result_l) == "Likely Increasing"

    # As Likely as Not
    result_aln = TrendResult(h=True, C=0.50, trend='decreasing')
    assert classify_trend(result_aln) == "As Likely as Not Decreasing"

def test_classify_custom_map():
    """Test classification with a user-provided custom category map."""
    custom_map = {
        0.99: "Virtually Certain",
        0.90: "Extremely Likely",
        0.50: "More Likely Than Not"
    }

    # Virtually Certain Increasing
    result_vc = TrendResult(h=True, C=0.995, trend='increasing')
    assert classify_trend(result_vc, category_map=custom_map) == "Virtually Certain Increasing"

    # Extremely Likely Decreasing
    result_el = TrendResult(h=True, C=0.95, trend='decreasing')
    assert classify_trend(result_el, category_map=custom_map) == "Extremely Likely Decreasing"

    # More Likely Than Not Increasing
    result_ml = TrendResult(h=True, C=0.60, trend='increasing')
    assert classify_trend(result_ml, category_map=custom_map) == "More Likely Than Not Increasing"

def test_classify_edge_cases():
    """Test edge cases, such as confidence values at the exact thresholds."""
    # Exactly on the threshold
    result_exact = TrendResult(h=True, C=0.90, trend='increasing')
    assert classify_trend(result_exact) == "Very Likely Increasing"

    # Just below a threshold
    result_below = TrendResult(h=True, C=0.89, trend='decreasing')
    assert classify_trend(result_below) == "Likely Decreasing"
