# pyMannKendall

This project provides a Python implementation of the Mann-Kendall test for trend analysis. The original code was sourced from the `pyMannKendall` package and has been modified to support unequally spaced time series data.

## MannKenSen Package

The `MannKenSen` package is a new addition to this project that provides modified versions of the Mann-Kendall test and Sen's slope estimator to handle unequally spaced time series data.

### `original_test(x, t, alpha=0.05)`

This function performs the Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data.
- `t`: A vector of timestamps corresponding to the data.
- `alpha`: The significance level (default is 0.05).

**Output:**
A named tuple with the following fields:
- `trend`: The trend of the data ('increasing', 'decreasing', or 'no trend').
- `h`: A boolean indicating whether the trend is significant.
- `p`: The p-value of the test.
- `z`: The Z-statistic.
- `Tau`: Kendall's Tau.
- `s`: The Mann-Kendall score.
- `var_s`: The variance of `s`.
- `slope`: The Sen's slope.
- `intercept`: The intercept of the trend line.
- `lower_ci`: The lower confidence interval of the slope.
- `upper_ci`: The upper confidence interval of the slope.

### `seasonal_test(x, t, period=12, alpha=0.05)`

This function performs the seasonal Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data.
- `t`: A vector of timestamps corresponding to the data.
- `period`: The seasonal period (default is 12).
- `alpha`: The significance level (default is 0.05).

**Output:**
A named tuple with the same fields as `original_test`.
