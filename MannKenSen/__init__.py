"""
MannKenSen: A Python package for non-parametric trend analysis on unequally spaced time series data.

This package provides implementations of the Mann-Kendall test and Sen's slope
estimator, with additional support for seasonal analysis, seasonality testing,
and plotting utilities.
"""
from .original_test import original_test
from .seasonal_test import seasonal_test
from .seasonality_test import seasonality_test
from .plotting import plot_seasonal_distribution
from .automated_seasonal_test import automated_seasonal_test

__all__ = [
    'original_test',
    'seasonal_test',
    'seasonality_test',
    'plot_seasonal_distribution',
    'automated_seasonal_test',
]
