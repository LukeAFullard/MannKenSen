# MannKenSen Package Audit

This document outlines a comprehensive audit of the `MannKenSen` Python package. The goal of this audit is to identify any potential issues in the codebase, statistical methodology, and overall structure to ensure its reliability as a scientific package.

## 1. Statistical Robustness

### 1.1 Censored Data Handling

- **`_mk_score_and_var_censored`**:
    - **Issue**: The function uses a tie-breaking method (`delx`, `dely`) that adds a small epsilon, calculated as `np.min(np.diff(unique_xx)) / 1000.0`. This approach is sensitive to the data's scale. For data with very small-magnitude values, this epsilon could be insignificant, failing to break ties. Conversely, for integer-based data, it may work as intended, but it's not a universally robust strategy.
    - **Status**: **Resolved**
    - **Resolution**: The epsilon calculation was made more robust. It is now calculated as half the minimum difference between unique values (`min_diff / 2.0`) when a difference exists, or defaults to `1.0` if all values are identical. This prevents floating-point underflow and ensures a reliable tie-breaking mechanism.
- **`_sens_estimator_censored`**:
    - **Issue**: The default `method='lwp'` sets ambiguous slopes to 0. This is a statistically biased choice, as it actively pulls the median slope towards zero. In a dataset with many ambiguous pairs, this can incorrectly suggest a "no trend" result.
    - **Status**: **Resolved**
    - **Resolution**: The default method was changed from `'lwp'` to `'nan'`. This makes the function's default behavior more statistically neutral by excluding ambiguous slopes from the calculation, rather than biasing the result towards zero. The user-facing functions `original_test` and `seasonal_test` were also updated to use this new, more robust default. The `'lwp'` method remains available as a non-default option.
- **`_aggregate_censored_median`**:
    - **Issue**: The logic for determining if a median is censored (`is_censored = median_val <= max_censored`) is a direct translation of the LWP-TRENDS heuristic. While useful for replication, it's not a statistically derived method and may not be robust, especially for complex data distributions.
    - **Status**: **Resolved**
    - **Resolution**: Added a note to the docstrings of `_aggregate_censored_median` and `seasonal_test` to clarify that the `'robust_median'` aggregation method uses a heuristic for determining censored medians, which is intended to replicate the LWP-TRENDS R script's behavior and may not be universally robust.

### 1.2 Confidence Intervals

- **`__confidence_intervals`**:
    - **Issue**: The function uses `np.interp` to find the confidence intervals. This linear interpolation assumes the slopes are uniformly distributed between ranks, which is not guaranteed. For small or heavily tied datasets, this can lead to inaccuracies. The standard, non-parametric method involves selecting the k-th smallest and (N-k+1)-th smallest slopes from the sorted list of slopes, where k is determined by the Z-statistic and variance.
    - **Recommendation**: Modify the function to use direct indexing of the sorted slopes. This involves calculating the integer ranks `M1` and `M2`, rounding them to the nearest integer, and using those as 0-based indices to select the correct slope values from the sorted array. This avoids the assumption of linearity made by interpolation.

## 2. Code Structure & Readability

- **`_utils.py`**:
    - **Issue**: This file is monolithic, containing over 500 lines of disparate functionality, including core statistics, data preparation, date handling, and aggregation logic. This makes the code hard to navigate, debug, and maintain.
    - **Recommendation**: Break `_utils.py` into smaller, more focused internal modules:
        - `_stats.py`: For core statistical calculations (`_mk_score_and_var_censored`, `_sens_estimator_censored`, `__confidence_intervals`, etc.).
        - `_helpers.py`: For data manipulation and aggregation helpers (`_prepare_data`, `_aggregate_by_group`, etc.).
        - `_datetime.py`: For time- and season-related functions (`_get_season_func`, `_get_cycle_identifier`, etc.).
- **Function Naming**:
    - **Issue**: Several internal functions use a double-underscore prefix (e.g., `__p_value`, `__sens_estimator_unequal_spacing`), which is typically reserved for name mangling in classes, not for private module functions.
    - **Recommendation**: Rename all internal helper functions to use a single leading underscore (e.g., `_p_value`, `_sens_estimator_unequal_spacing`) to follow standard Python conventions for "internal use" functions.
- **Dead Code**:
    - **Issue**: The functions `__mk_score` and `__variance_s` appear to be legacy code from a non-censored implementation. The main workflows use `_mk_score_and_var_censored`, which contains its own, more complex variance calculation. This dead code adds clutter and potential confusion.
    - **Recommendation**: Remove the unused `__mk_score` and `__variance_s` functions to clean up the codebase.

## 3. Error Handling & Validation

- **Input Validation**:
    - **Issue**: User-facing functions like `original_test` and `seasonal_test` do not validate string-based enum parameters like `sens_slope_method`, `tau_method`, or `agg_method`. A user typo (e.g., `agg_method='medain'`) would not raise an immediate error but would cause the logic to default to the `'none'` case, leading to silent failure and incorrect results.
    - **Recommendation**: Add explicit validation checks at the beginning of `original_test` and `seasonal_test` to ensure that these parameters are one of the accepted values. Raise a `ValueError` if an invalid option is provided.
- **`prepare_censored_data`**:
    - **Issue**: The error messages for malformed censored strings are not very descriptive. For example, inputting `'<'` raises `ValueError: Invalid left-censored value format: '<'`, which doesn't explain *why* it's invalid (i.e., missing a number).
    - **Recommendation**: Enhance the error messages to be more informative. For example: `ValueError: Invalid left-censored value format: '<'. Expected a number after the '<' symbol.`

## 4. Testing

- **Coverage Gaps**:
    - **Issue**: The test suite is missing coverage for several important scenarios:
        - The `hicensor=True` functionality is not tested at all.
        - The various `agg_method` options in `seasonal_test` are not explicitly tested.
        - Edge cases in `prepare_censored_data` (e.g., malformed strings like `'<>5'`, strings with extra spaces) are not tested.
        - Behavior with empty inputs or inputs containing only `np.nan`.
    - **Recommendation**: Add specific unit tests to cover these scenarios to ensure the package is robust and reliable.
- **Statistical Validation**:
    - **Issue**: There are no integration tests that validate the statistical output against a known, trusted implementation. While the user has instructed not to use `rpy2` for this, it represents a major gap in ensuring the scientific validity of the package's results.
    - **Recommendation**: As a temporary measure, create static test cases using data and results generated *manually* from the LWP-TRENDS R script. These "snapshot" tests would provide a baseline for correctness, even without a direct R dependency. This should be noted as a limitation.

## 5. Documentation

- **Docstrings**:
    - **Issue**: The docstrings for user-facing functions often refer to internal functions for implementation details. For example, the `sens_slope_method` parameter in `original_test` says "See `_sens_estimator_censored` for details." This is unhelpful for a user who should not need to read internal code.
    - **Recommendation**: Update the docstrings of all public functions (`original_test`, `seasonal_test`, etc.) to be self-contained. The explanation for parameters like `sens_slope_method` should be copied into the user-facing docstring.
- **`README.md`**:
    - **Issue**: The README lacks a "Limitations" section and a clear, prominent example of how to handle censored data, which is a key feature of the package.
    - **Recommendation**: Add a "Limitations" section that explicitly states which features from the LWP-TRENDS R script are *not* implemented (e.g., covariate adjustment). Add a dedicated "Censored Data Analysis" section with a clear code example using `prepare_censored_data`.
