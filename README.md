# MannKenSen

This project provides a Python implementation of the Mann-Kendall test for trend analysis. The original code was sourced from the `pyMannKendall` package and has been modified to support unequally spaced time series data. The `pyMannKendall` code will be deleted once we are production ready.

## Background and Acknowledgements

The statistical methods used in this package, particularly for handling unequally spaced time series and seasonal aggregation, are inspired by the LWP-TRENDS R package developed by Land & Water People in New Zealand. Their robust implementation has served as an excellent reference for this work.

For more information on the original R functions, please see: [LWPTrends_v2502.zip](https://landwaterpeople.co.nz/wp-content/uploads/2025/03/LWPTrends_v2502.zip)

## MannKenSen Package

The `MannKenSen` package provides modified versions of the Mann-Kendall test and Sen's slope estimator to handle unequally spaced time series data.

### `original_test(x, t, alpha=0.05)`

This function performs the Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data.
- `t`: A vector of timestamps corresponding to the data. The function automatically handles numeric vectors, as well as `numpy.datetime64` and Python `datetime` objects by converting them to Unix timestamps.
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
- `C`: The confidence of the trend direction.
- `Cd`: The confidence that the trend is decreasing.


### `seasonal_test(x, t, period=12, alpha=0.05, agg_method='none', season_type='month')`

This function performs the seasonal Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data.
- `t`: A vector of timestamps. Can be a numeric vector or a datetime-like array.
- `period`: The seasonal period. For numeric `t`, this defines the cycle length. Cycles are calculated relative to the start of the time series (e.g., for `t=[1990, 1991, 1992]` and `period=1`, each year is a cycle). For datetime `t`, this must match the expected period of the `season_type`.
- `alpha`: The significance level (default is 0.05).
- `agg_method`: The method for aggregating multiple data points within the same season-year.
  - `'none'` (default): Performs the analysis on all individual data points.
  - `'median'`: Aggregates data using the median of the values and times.
  - `'middle'`: Aggregates data by selecting the observation closest to the middle of the time period.
- `season_type`: For datetime inputs, specifies how to define a season. See the table below for options.

**Output:**
A named tuple with the same fields as `original_test`.

#### `season_type` Options for Datetime Inputs

| `season_type`    | Description                               | Expected `period` | How it Groups Data |
|------------------|-------------------------------------------|-------------------|--------------------|
| `'year'`         | Annual Trend (non-seasonal)               | 1                 | Groups all data into a single season. |
| `'month'`        | Month of the year                         | 12                | Groups data by the calendar month (e.g., all Januaries). |
| `'day_of_week'`  | Day of the week                           | 7                 | Groups data by the day of the week (e.g., all Mondays). |
| `'quarter'`      | Quarter of the year                       | 4                 | Groups data by the calendar quarter (e.g., all Q1s). |
| `'hour'`         | Hour of the day                           | 24                | Groups data by the hour of the day (e.g., all 9 AMs). |
| `'week_of_year'` | ISO week of the year                      | 52 or 53          | Groups data by the ISO week number. |
| `'day_of_year'`  | Day of the year                           | (no validation)   | Groups data by the day of the year (1-366). |
| `'minute'`       | Minute of the hour                        | 60                | Groups data by the minute of the hour. |
| `'second'`       | Second of the minute                      | 60                | Groups data by the second of the minute. |


**Example: Weekly Seasonality**
```python
import numpy as np
import pandas as pd
from MannKenSen import seasonal_test

# Create 4 years of weekly data
t = pd.to_datetime(pd.date_range(start='2020-01-01', periods=208, freq='W'))
x = np.random.rand(208)

# Inject a clear increasing trend in the 10th week of each year
tenth_week_mask = t.isocalendar().week == 10
x[tenth_week_mask] = [10, 20, 30, 40]

# Perform the test for weekly seasonality
result = seasonal_test(x, t, period=53, season_type='week_of_year')
print(result)
```

### `seasonality_test(x, t, period=12, alpha=0.05, season_type='month')`

Performs a Kruskal-Wallis H-test to determine if there is a statistically significant difference between the distributions of data across seasons. This is a common way to test for the presence of seasonality.

**Input:**
- `x`: A vector of data.
- `t`: A vector of timestamps.
- `period`: The seasonal period.
- `alpha`: The significance level to determine if the result is seasonal.
- `season_type`: For datetime inputs, specifies how to define a season.

**Output:**
A named tuple with the fields: `h_statistic`, `p_value`, and `is_seasonal` (a boolean).

**Example: Testing for Seasonality**
```python
# Using the same seasonal data from the previous example
is_seasonal_result = seasonality_test(x, t, period=53, season_type='week_of_year')
print(is_seasonal_result)
# Returns: Seasonality_Test(h_statistic=..., p_value=..., is_seasonal=True)
```

### `plot_seasonal_distribution(x, t, period=12, season_type='month', save_path='seasonal_distribution.png')`

Generates and saves a box plot to visually compare the distribution of values across different seasons. This is a helpful utility for visually inspecting seasonality.

**Input:**
- Accepts the same `x`, `t`, `period`, and `season_type` parameters as `seasonal_test`.
- `save_path`: The file path where the plot image will be saved.

**Output:**
The file path where the plot was saved.

**Example: Plotting Seasonal Distributions**
```python
# Using the same seasonal data from the previous example
plot_path = plot_seasonal_distribution(x, t, period=53, season_type='week_of_year', save_path='my_seasonal_plot.png')
print(f"Plot saved to: {plot_path}")
```
