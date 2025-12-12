# MannKenSen

This project provides a Python implementation of the Mann-Kendall test for trend analysis. The original code was sourced from the `pyMannKendall` package and has been modified to support unequally spaced time series data and regional trend aggregation. The `pyMannKendall` code will be deleted once we are production ready.

## Background and Acknowledgements

The statistical methods used in this package, particularly for handling unequally spaced time series and seasonal aggregation, are inspired by the LWP-TRENDS R package developed by Land & Water People in New Zealand. Their robust implementation has served as an excellent reference for this work.

For more information on the original R functions, please see: [LWPTrends_v2502.zip](https://landwaterpeople.co.nz/wp-content/uploads/2025/03/LWPTrends_v2502.zip)

## Installation

To install the necessary dependencies for this package, run the following command:

```bash
pip install -r requirements.txt
```

## Working with Censored Data

The `MannKenSen` package is designed to work with censored data, where some measurements are only known to be above or below a certain detection limit.

### `prepare_censored_data(x)`

To use this feature, your data must first be pre-processed into a specific `pandas.DataFrame` format. The package provides the `prepare_censored_data` utility for this purpose. It converts an array of numbers and censored strings (e.g., `'<5'`) into the required DataFrame structure.

**Input:**
- `x` (array-like): A 1D array or list containing numeric values and/or censored strings.
  - Left-censored values should be formatted as `'<value'`.
  - Right-censored values should be formatted as `'>value'`.

**Output:**
A `pandas.DataFrame` with the following columns:
- `'value'`: The numeric value of the detection limit.
- `'censored'`: A boolean, `True` if the value was censored.
- `'cen_type'`: The type of censoring (`'lt'`, `'gt'`, or `'not'`).

This pre-processed DataFrame can be passed directly to the `original_test` and `seasonal_test` functions.

**Example: Preparing Censored Data**
```python
import numpy as np
import pandas as pd
from MannKenSen import prepare_censored_data, original_test

# Create a time vector
t = pd.to_datetime(pd.date_range(start='2010-01-01', periods=10, freq='YE'))

# Create a data vector with mixed numeric and censored string values
x_raw = [5, '<4', 3.5, '>6', 6.2, '<4', 3, 2.5, '<2', 2.1]

# Pre-process the censored data
x_prepared = prepare_censored_data(x_raw)

print(x_prepared)
# Expected Output:
#    value  censored cen_type
# 0    5.0     False      not
# 1    4.0      True       lt
# 2    3.5     False      not
# 3    6.0      True       gt
# 4    6.2     False      not
# 5    4.0      True       lt
# 6    3.0     False      not
# 7    2.5     False      not
# 8    2.0      True       lt
# 9    2.1     False      not

# The prepared DataFrame can now be used in the trend test
result = original_test(x=x_prepared, t=t)
print(result)
```

## MannKenSen Package

The `MannKenSen` package provides modified versions of the Mann-Kendall test and Sen's slope estimator to handle unequally spaced time series data.

### `original_test(x, t, alpha=0.05, hicensor=False, plot_path=None, lt_mult=0.5, gt_mult=1.1, sens_slope_method='lwp', tau_method='b', agg_method='none')`

This function performs the Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data, or a pre-processed DataFrame from `prepare_censored_data`.
- `t`: A vector of timestamps corresponding to `x`.
- `alpha`: The significance level (default is 0.05).
- `hicensor` (bool): If `True`, applies the high-censor rule, where all values below the highest left-censor limit are treated as censored at that limit. Default is `False`.
- `plot_path` (str, optional): If a file path is provided, a plot of the trend analysis is saved. Default is `None`.
- `lt_mult` (float): The multiplier for left-censored data in the Sen's slope calculation (default is 0.5).
- `gt_mult` (float): The multiplier for right-censored data in the Sen's slope calculation (default is 1.1).
- `sens_slope_method` (str): The method for handling ambiguous slopes in censored data.
  - `'lwp'` (default): Sets ambiguous slopes to 0, mimicking the LWP-TRENDS R script.
  - `'nan'`: Sets ambiguous slopes to `np.nan`, a more statistically neutral approach.
- `tau_method` (str): The method for calculating Kendall's Tau ('a' or 'b'). Default is `'b'`, which accounts for ties in the data and is the recommended method.
- `agg_method` (str): The method for aggregating data at tied timestamps.
  - `'none'` (default): No aggregation is performed. A warning is issued if ties are present, as this can affect the Sen's slope calculation.
  - `'median'`, `'robust_median'`, `'middle'`: See `seasonal_test` for descriptions. It is recommended to use an aggregation method when tied timestamps are present.

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


### `seasonal_test(x, t, period=12, alpha=0.05, agg_method='none', season_type='month', hicensor=False, plot_path=None, lt_mult=0.5, gt_mult=1.1, sens_slope_method='lwp', tau_method='b')`

This function performs the seasonal Mann-Kendall test on unequally spaced time series data.

**Input:**
- `x`: A vector of data, or a pre-processed DataFrame from `prepare_censored_data`.
- `t`: A vector of timestamps.
- `period`: The seasonal period (default is 12).
- `alpha`: The significance level (default is 0.05).
- `agg_method`: The method for aggregating data within a season-cycle.
  - `'none'` (default): No aggregation is performed.
  - `'median'`: (LWP method) A simple median aggregation.
  - `'robust_median'`: A more robust median for censored data.
  - `'middle'`: Selects the observation closest to the middle of the time period.
- `season_type`: For datetime inputs, specifies how to define a season (default is `'month'`).
- `hicensor` (bool): If `True`, applies the high-censor rule. Default is `False`.
- `plot_path` (str, optional): If provided, saves a plot of the trend analysis.
- `lt_mult` (float): Multiplier for left-censored data (default is 0.5).
- `gt_mult` (float): Multiplier for right-censored data (default is 1.1).
- `sens_slope_method` (str): The method for handling ambiguous slopes in censored data.
  - `'lwp'` (default): Sets ambiguous slopes to 0.
  - `'nan'`: Sets ambiguous slopes to `np.nan`.
- `tau_method` (str): The method for calculating Kendall's Tau ('a' or 'b'). Default is `'b'`, which accounts for ties in the data and is the recommended method.

**Output:**
A named tuple with the same fields as `original_test`.

#### Understanding the `period` Parameter

The `period` parameter is crucial for correct seasonal analysis, and its meaning depends on the type of the time vector `t` you provide.

*   **For Datetime-like Inputs:**
    When `t` is a vector of `datetime` or `numpy.datetime64` objects, the `period` parameter serves as a **validation mechanism** for the chosen `season_type`. For instance, if `season_type` is `'month'`, the package expects `period` to be 12. If a different value is provided, it will raise a `ValueError`. This ensures that the seasonal analysis aligns with the calendar definition of the seasonality. For `season_type`s with variable periods, like `'week_of_year'`, the `period` can be either 52 or 53.

*   **For Numeric Inputs:**
    When `t` is a numeric vector (e.g., years, days since an event), the `period` defines the **length of a full seasonal cycle**. The interpretation of the `period` is entirely dependent on the units of your time vector.
    - If `t` is in **days**, a `period` of `365` or `365.25` would represent a yearly seasonal cycle.
    - If `t` is in **months**, a `period` of `12` would represent a yearly cycle.
    - If `t` is in **years** and you want to test for a biennial cycle, a `period` of `2` would be appropriate.

*   **Caveats and Limitations:**
    - The package **cannot infer the units** of a numeric time vector.
    - It is the **user's responsibility** to provide a `period` that is meaningful and correctly corresponds to the units of their numeric `t` vector.
    - Incorrectly specifying the `period` for numeric data will lead to a mathematically valid but contextually incorrect analysis.

#### `season_type` Options for Datetime Inputs

| `season_type`    | Description                               | Expected `period` | How it Groups Data |
|------------------|-------------------------------------------|-------------------|--------------------|
| `'year'`         | Annual Trend (non-seasonal)               | 1                 | Groups all data into a single season. |
| `'month'`        | Month of the year                         | 12                | Groups data by the calendar month (e.g., all Januaries). |
| `'day_of_week'`  | Day of the week                           | 7                 | Groups data by the day of the week (e.g., all Mondays). |
| `'quarter'`      | Quarter of the year                       | 4                 | Groups data by the calendar quarter (e.g., all Q1s). |
| `'hour'`         | Hour of the day                           | 24                | Groups data by the hour of the day (e.g., all 9 AMs). |
| `'week_of_year'` | ISO week of the year                      | 52 or 53          | Groups data by the ISO week number. |
| `'biweekly'`     | Two-week period of the year               | 26 or 27          | Groups data by two-week periods based on the ISO week number. |
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
result = seasonal_test(x, t, period=52, season_type='week_of_year')
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
is_seasonal_result = seasonality_test(x, t, period=52, season_type='week_of_year')
print(is_seasonal_result)
# Returns: Seasonality_Test(h_statistic=..., p_value=..., is_seasonal=True)
```

### Visualizing Trend Analysis

Both the `original_test` and `seasonal_test` functions include a `plot_path` parameter that allows you to generate and save a visualization of the trend analysis.

**Example: Generating a Trend Plot**
```python
# Using the data from the original_test example
result = original_test(x=x_prepared, t=t, plot_path='trend_analysis.png')
print("Trend plot saved to trend_analysis.png")
```

**Interpreting the Plot:**
- **Data Points:** Non-censored data points are shown as blue circles, while censored data points are marked with red 'x's.
- **Sen's Slope:** The calculated trend is displayed as a dashed black line.
- **Confidence Intervals:** The shaded gray area represents the confidence interval of the slope (e.g., 95% CI for `alpha=0.05`), providing a visual representation of the uncertainty in the trend.
- **Statistics Box:** A text box in the top-left corner provides a summary of the key statistical results: the trend direction, Kendall's Tau, the Sen's slope, and the p-value.

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
plot_path = plot_seasonal_distribution(x, t, period=52, season_type='week_of_year', save_path='my_seasonal_plot.png')
print(f"Plot saved to: {plot_path}")
```

### `analysis_notes` Module

This module provides functions for data quality checks inspired by the LWP-TRENDS R script. These checks generate "Analysis Notes" to warn the user about potential issues with their data that could affect trend analysis results.

#### `get_analysis_note(data, values_col='value', censored_col='censored', is_seasonal=False, post_aggregation=False, season_col='season')`

Performs data quality checks and returns an analysis note.

**Input:**
- `data`: A pandas DataFrame containing the data.
- `values_col`: The name of the column with data values.
- `censored_col`: The name of the boolean column indicating censored data.
- `is_seasonal`: A boolean indicating whether to perform seasonal checks.
- `post_aggregation`: A boolean indicating whether the checks are for after data aggregation.
- `season_col`: The name of the column with season identifiers.

**Output:**
A string containing the analysis note (e.g., "ok", "< 3 unique values").

#### `get_sens_slope_analysis_note(slopes, t, cen_type)`

Provides warnings about censored values used in the derivation of the Sen's Slope.

**Input:**
- `slopes`: An array of calculated slopes.
- `t`: An array of timestamps.
- `cen_type`: An array of censor types ('lt', 'gt', 'not').

**Output:**
A string containing the analysis note (e.g., "ok", "WARNING: Sen slope influenced by censored values").

### `inspect_trend_data(data, trend_period=None, end_year=None, prop_year_tol=0.9, prop_incr_tol=0.9, return_summary=False, custom_increments=None)`

Inspects data availability over a trend period and determines the best time increment for trend analysis.

**Input:**
- `data`: A pandas DataFrame containing at least a time column and a value column.
- `trend_period`: The number of years in the trend period. Defaults to the full range of data.
- `end_year`: The last year of the trend period. Defaults to the last year in the data.
- `prop_year_tol`: The minimum proportion of years in the trend period that must have at least one observation.
- `prop_incr_tol`: The minimum proportion of time increments within the trend period that must have at least one observation.
- `return_summary`: If True, returns a namedtuple containing the modified DataFrame and a summary of data availability.
- `custom_increments`: A dictionary of custom time increments. Keys are the increment names (e.g., 'weekly'), and values are the number of increments in a year. If not provided, a default set of increments is used. Supported increments are: 'annually', 'bi-annually', 'quarterly', 'bi-monthly', 'monthly', 'weekly', 'daily'.

**Output:**
- The filtered DataFrame with a new 'time_increment' column. If no suitable increment is found, this column will be filled with 'none'.
- If `return_summary` is True, a namedtuple `InspectionResult` with `data` and `summary` attributes is returned.

### `regional_test(trend_results, time_series_data, site_col='site', value_col='value', time_col='time', s_col='s', c_col='C')`

Performs a regional trend aggregation analysis on the results of multiple single-site trend tests. This is useful for determining if there is a consistent, statistically significant trend across a network of monitoring sites. The methodology corrects for inter-site correlation, providing a more robust regional assessment.

**Input:**
- `trend_results`: A pandas DataFrame containing the results from running `original_test` or `seasonal_test` on multiple sites. Must contain columns for the site identifier, the Mann-Kendall score `s`, and the confidence `C`.
- `time_series_data`: A pandas DataFrame containing the original time series data for all sites. This is used to calculate the inter-site correlation.
- `site_col`: The name of the site identifier column in both DataFrames.
- `value_col`: The name of the value column in `time_series_data`.
- `time_col`: The name of the time column in `time_series_data`.
- `s_col`: The name of the Mann-Kendall score column in `trend_results`.
- `c_col`: The name of the confidence column in `trend_results`.

**Output:**
A named tuple with the following fields:
- `M`: The total number of sites.
- `TAU`: The aggregate trend strength (proportion of sites trending in the modal direction).
- `VarTAU`: The uncorrected variance of `TAU`.
- `CorrectedVarTAU`: The variance of `TAU` corrected for inter-site correlation.
- `DT`: The aggregate trend direction ('Increasing' or 'Decreasing').
- `CT`: The confidence in the aggregate trend direction.

**Example: Regional Trend Aggregation**
```python
import pandas as pd
import numpy as np
from MannKenSen import original_test, regional_test

# 1. Create synthetic data for three sites
dates = pd.to_datetime(pd.date_range(start='2000-01-01', periods=20, freq='YE'))
sites = ['A', 'B', 'C']
all_ts_data = []

for site in sites:
    np.random.seed(hash(site) % (2**32 - 1))
    noise = np.random.normal(0, 1.0, 20)
    # Site B has a decreasing trend, A and C are increasing
    trend = -0.1 * np.arange(20) if site == 'B' else 0.1 * np.arange(20)
    values = 10 + trend + noise
    df = pd.DataFrame({'time': dates, 'value': values, 'site': site})
    all_ts_data.append(df)

time_series_data = pd.concat(all_ts_data, ignore_index=True)

# 2. Run single-site trend analysis for each site
trend_results = []
for site in sites:
    site_data = time_series_data[time_series_data['site'] == site]
    res = original_test(x=site_data['value'], t=site_data['time'])
    res_dict = res._asdict()
    res_dict['site'] = site
    trend_results.append(res_dict)

trend_results_df = pd.DataFrame(trend_results)

# 3. Run the regional trend aggregation
regional_res = regional_test(trend_results=trend_results_df,
                             time_series_data=time_series_data)

print(regional_res)
# Expected Output (values will vary slightly due to random noise):
# RegionalTrendResult(M=3, TAU=0.666..., VarTAU=..., CorrectedVarTAU=..., DT='Increasing', CT=...)
```

### `classify_trend(result, category_map=None)`

Classifies a trend result into a descriptive, human-readable category based on its statistical significance and confidence. This is useful for interpreting and communicating trend analysis results.

**Input:**
- `result`: The namedtuple returned by `original_test` or `seasonal_test`.
- `category_map` (dict, optional): A dictionary mapping confidence thresholds (float) to descriptive category labels (str). If `None`, a default IPCC-style mapping is used:
  ```python
  {
      0.95: "Highly Likely",
      0.90: "Very Likely",
      0.67: "Likely",
      0.0:  "As Likely as Not"
  }
  ```

**Output:**
A string describing the trend category (e.g., "Highly Likely Increasing", "No Trend").

**Example: Classifying a Trend**
```python
import numpy as np
from MannKenSen import original_test, classify_trend

# Create synthetic data with a clear trend
t = np.arange(20)
x = 0.1 * t + np.random.normal(0, 0.5, 20)

# Perform the trend test
result = original_test(x, t)

# Classify the result
category = classify_trend(result)
print(f"The trend is: {category}")
# Possible output: "The trend is: Very Likely Increasing"

# Example with a custom category map
custom_map = {
    0.99: "Virtually Certain",
    0.90: "Extremely Likely",
}
custom_category = classify_trend(result, category_map=custom_map)
print(f"The custom category is: {custom_category}")
# Possible output: "The custom category is: Extremely Likely Increasing"
```
