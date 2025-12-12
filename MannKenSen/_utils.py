"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""

from __future__ import division
import datetime
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata




def _is_datetime_like(x):
    """Checks if an array is datetime-like."""
    return np.issubdtype(x.dtype, np.datetime64) or \
           (x.dtype == 'O' and len(x) > 0 and hasattr(x[0], 'year'))

def _get_season_func(season_type, period):
    """
    Returns a function to extract seasonal data based on the season_type,
    and validates the period.
    """
    def get_dt_prop(dt, prop):
        return getattr(dt.dt, prop) if isinstance(dt, pd.Series) else getattr(dt, prop)

    season_map = {
        'year': (1, lambda dt: get_dt_prop(dt, 'year')),
        'month': (12, lambda dt: get_dt_prop(dt, 'month')),
        'day_of_week': (7, lambda dt: get_dt_prop(dt, 'dayofweek')),
        'quarter': (4, lambda dt: get_dt_prop(dt, 'quarter')),
        'hour': (24, lambda dt: get_dt_prop(dt, 'hour')),
        'week_of_year': ([52, 53], lambda dt: get_dt_prop(dt, 'isocalendar')().week),
        'biweekly': ([26, 27], lambda dt: (get_dt_prop(dt, 'isocalendar')().week - 1) // 2),
        'day_of_year': (None, lambda dt: get_dt_prop(dt, 'dayofyear')),
        'minute': (60, lambda dt: get_dt_prop(dt, 'minute')),
        'second': (60, lambda dt: get_dt_prop(dt, 'second')),
    }
    if season_type not in season_map:
        raise ValueError(f"Unknown season_type: '{season_type}'. Must be one of {list(season_map.keys())}")

    expected_period, season_func = season_map[season_type]

    if expected_period is not None:
        if isinstance(expected_period, list):
            if period not in expected_period:
                raise ValueError(f"For season_type='{season_type}', period must be one of {expected_period}.")
        elif period != expected_period:
            raise ValueError(f"For season_type='{season_type}', period must be {expected_period}.")

    return season_func

def _get_cycle_identifier(dt_series, season_type):
    """
    Returns a numeric series that uniquely identifies the larger time cycle
    for each timestamp, used for aggregation.
    """
    dt_accessor = dt_series.dt if isinstance(dt_series, pd.Series) else dt_series

    if season_type in ['month', 'quarter', 'year', 'day_of_year', 'week_of_year', 'biweekly']:
        # The cycle is the year
        return dt_accessor.year.to_numpy()

    elif season_type == 'day_of_week':
        # The cycle is the week, identified by year and week number
        iso_cal = dt_accessor.isocalendar()
        return (iso_cal.year * 100 + iso_cal.week).to_numpy()

    elif season_type in ['hour', 'minute', 'second']:
        # The cycle is the day, identified by the Unix timestamp of the day's start
        # Convert to int64 (nanoseconds) and then to float seconds
        return (dt_accessor.normalize().astype(np.int64) / 10**9)

    else:
        # Default to year if the concept of a cycle is not obvious
        return dt_accessor.year.to_numpy()


def _get_time_ranks(t_values, cycles):
    """Convert timestamps to cycle-based ranks matching R implementation."""
    # Create unique cycle identifiers and sort them to ensure rank order
    unique_cycles = np.unique(cycles)
    ranks = np.zeros_like(t_values, dtype=float)

    # Assign sequential ranks to each cycle
    for i, cycle in enumerate(unique_cycles, start=1):
        mask = cycles == cycle
        ranks[mask] = i

    return ranks


def _rle_lengths(a):
    """
    Calculates the lengths of runs of equal values in an array.
    Equivalent to R's `rle(x)$lengths`.
    """
    if len(a) == 0:
        return np.array([], dtype=int)
    y = a[1:] != a[:-1]
    i = np.append(np.where(y), len(a) - 1)
    return np.diff(np.append(-1, i))


def _mk_score_and_var_censored(x, t, censored, cen_type, tau_method='b'):
    """
    Calculates the Mann-Kendall S statistic and its variance for censored data.
    This is a Python translation of the GetKendal function from the LWP-TRENDS
    R script, which is adapted from the NADA package (Helsel, 2012).
    """
    # 1. Special handling for right-censored ('gt') data
    x_mod = x.copy()
    censored_mod = censored.copy()
    if np.any(cen_type == 'gt'):
        gt_mask = cen_type == 'gt'

        # Calculate a small, data-relative value to break ties
        unique_vals = np.unique(x_mod)
        if len(unique_vals) > 1:
            min_diff = np.min(np.diff(unique_vals))
            # Use a small fraction of the minimum difference, or a default small number
            tie_break_value = min_diff * 0.01 if min_diff > 0 else 1e-9
        else:
            tie_break_value = 1e-9 # Fallback for data with no variance

        # Add the small amount to break ties, treat as uncensored
        max_gt_val = x_mod[gt_mask].max() + tie_break_value
        x_mod[gt_mask] = max_gt_val
        censored_mod[gt_mask] = False

    # 2. Prepare inputs
    xx = x_mod
    cx = censored_mod.astype(bool)
    # Time is treated as uncensored
    yy = rankdata(t, method='ordinal')
    cy = np.zeros_like(yy, dtype=bool)
    n = len(xx)

    if n < 2:
        return 0, 0, 0

    # 3. Calculate delx and dely to break ties
    unique_xx = np.unique(xx)
    delx = np.min(np.diff(unique_xx)) / 1000.0 if len(unique_xx) > 1 else 0.0
    unique_yy = np.unique(yy)
    dely = np.min(np.diff(unique_yy)) / 1000.0 if len(unique_yy) > 1 else 0.0

    # 4. S-statistic calculation using vectorized outer products
    dupx = xx - delx * cx
    dupy = yy - dely * cy

    diffx = dupx[:, np.newaxis] - dupx
    diffy = dupy[:, np.newaxis] - dupy
    signyx = np.sign(diffy * diffx)

    diffcx = cx[:, np.newaxis].astype(int) - cx.astype(int)
    cix = np.sign(diffcx) * np.sign(diffx)
    cix[cix <= 0] = 0
    signyx *= (1 - cix)

    xplus = (cx[:, np.newaxis].astype(int) + cx.astype(int))
    xplus[xplus <= 1] = 0
    xplus[xplus > 1] = 1
    tplus = xplus * np.abs(np.sign(diffx))


    itot = np.sum(np.triu(signyx * (1 - xplus), k=1))
    kenS = itot

    # 5. D (denominator) calculation for Tau
    J = n * (n - 1) / 2.0
    if tau_method == 'a':
        D = J
    else: # Default to Tau-b
        # tt: number of tied pairs in x
        tt = (np.sum(1 - np.abs(np.sign(diffx))) - n) / 2.0
        tt += np.sum(cix) / 2.0
        tt += np.sum(tplus) / 2.0

        # uu: number of tied pairs in y (time)
        diffcy = cy[:, np.newaxis].astype(int) - cy.astype(int)
        ciy = np.sign(diffcy) * np.sign(diffy)
        ciy[ciy <= 0] = 0
        uu = (np.sum(1 - np.abs(np.sign(diffy))) - n) / 2.0
        uu += np.sum(ciy) / 2.0
        yplus = (cy[:, np.newaxis].astype(int) + cy.astype(int))
        yplus[yplus <= 1] = 0
        yplus[yplus > 1] = 1
        uplus = yplus * np.abs(np.sign(diffy))
        uu += np.sum(uplus) / 2.0

        tau_denom = np.sqrt(J - tt) * np.sqrt(J - uu)
        D = tau_denom if tau_denom > 0 else J


    # 6. Variance Calculation (adapted from NADA::cenken)
    varS = n * (n - 1) * (2 * n + 5) / 18.0

    # Add tie correction for x variable (previously in __variance_s)
    unique_x, tp = np.unique(x, return_counts=True)
    if n != len(unique_x):
        varS -= np.sum(tp * (tp - 1) * (2 * tp + 5)) / 18.0

    intg = np.arange(1, n + 1)

    dorder_x = np.argsort(dupx)
    dxx = dupx[dorder_x]
    dcx = cx[dorder_x]

    dorder_y = np.argsort(dupy)
    dyy = dupy[dorder_y]
    dcy = cy[dorder_y]

    # delc correction for ties
    tmpx = dxx - intg * (1 - dcx) * delx
    tmpy = dyy - intg * (1 - dcy) * dely
    rxlng = _rle_lengths(rankdata(tmpx, method='ordinal'))
    rylng = _rle_lengths(rankdata(tmpy, method='ordinal'))

    def var_adj_term(lng):
        lng_vals, lng_counts = np.unique(lng, return_counts=True)
        t1 = np.sum(lng_counts * lng_vals * (lng_vals - 1) * (2 * lng_vals + 5))
        t2 = np.sum(lng_counts * lng_vals * (lng_vals - 1) * (lng_vals - 2))
        t3 = np.sum(lng_counts * lng_vals * (lng_vals - 1))
        return t1, t2, t3

    x1, x2, x3 = var_adj_term(rxlng)
    y1, y2, y3 = var_adj_term(rylng)

    term2 = (x2 * y2) / (9.0 * n * (n - 1) * (n - 2)) if n > 2 else 0
    term3 = (x3 * y3) / (2.0 * n * (n - 1))
    delc = (x1 + y1) / 18.0 - term2 - term3

    # deluc correction for uncensored-censored ties
    x4 = x3
    y4 = y3
    tmpx_uc = intg * dcx - 1
    tmpx_uc[tmpx_uc < 0] = 0
    nrxlng_uc = np.sum(tmpx_uc)
    x1_uc = nrxlng_uc * 2 * 1 * (2 * 2 + 5)
    x2_uc = 0 # (2-2) = 0
    x3_uc = nrxlng_uc * 2 * 1

    tmpy_uc = intg * dcy - 1
    tmpy_uc[tmpy_uc < 0] = 0
    nrylng_uc = np.sum(tmpy_uc)
    y1_uc = nrylng_uc * 2 * 1 * (2 * 2 + 5)
    y2_uc = 0
    y3_uc = nrylng_uc * 2 * 1

    term2_uc = (x2_uc * y2_uc) / (9.0 * n * (n - 1) * (n - 2)) if n > 2 else 0
    term3_uc = (x3_uc * y3_uc) / (2.0 * n * (n - 1))
    deluc = (x1_uc + y1_uc) / 18.0 - term2_uc - term3_uc - (x4 + y4)

    # delu correction for censored-censored ties
    dxx_u = dxx - intg * dcx * delx
    dyy_u = dyy - intg * dcy * dely
    rxlng_u = _rle_lengths(rankdata(dxx_u, method='ordinal'))
    rylng_u = _rle_lengths(rankdata(dyy_u, method='ordinal'))
    x1_u, x2_u, x3_u = var_adj_term(rxlng_u)
    y1_u, y2_u, y3_u = var_adj_term(rylng_u)

    term2_u = (x2_u * y2_u) / (9.0 * n * (n - 1) * (n - 2)) if n > 2 else 0
    term3_u = (x3_u * y3_u) / (2.0 * n * (n - 1))
    delu = (x1_u + y1_u) / 18.0 - term2_u - term3_u

    varS = varS - delc - deluc - delu
    return kenS, varS, D


# Helper functions from the original pymannkendall.py
def __preprocessing(x):
    x = np.asarray(x)

    # Convert datetime objects to numeric timestamps if necessary
    if _is_datetime_like(x):
        x = x.astype('datetime64[s]').astype(float)
    elif x.dtype == 'O' and len(x) > 0:
        if isinstance(x[0], datetime.datetime):
            x = np.array([val.timestamp() for val in x])

    x = x.astype(float)

    if x.ndim == 2:
        (n, c) = x.shape
        if c == 1:
            x = x.flatten()
    return x, (1 if x.ndim == 1 else x.shape[1])

def __missing_values_analysis(x, method='skip'):
    if method.lower() == 'skip':
        if x.ndim == 1:
            x = x[~np.isnan(x)]
        else:
            x = x[~np.isnan(x).any(axis=1)]
    return x, len(x)

def __mk_score(x, n):
    s = 0
    demo = np.ones(n)
    for k in range(n - 1):
        s += np.sum(demo[k + 1:n][x[k + 1:n] > x[k]]) - np.sum(demo[k + 1:n][x[k + 1:n] < x[k]])
    return s

def __variance_s(x, n):
    unique_x, tp = np.unique(x, return_counts=True)
    if n == len(unique_x):
        return (n * (n - 1) * (2 * n + 5)) / 18
    return (n * (n - 1) * (2 * n + 5) - np.sum(tp * (tp - 1) * (2 * tp + 5))) / 18

def __z_score(s, var_s):
    if var_s == 0:
        return 0
    if s > 0:
        return (s - 1) / np.sqrt(var_s)
    return (s + 1) / np.sqrt(var_s) if s < 0 else 0

def __p_value(z, alpha):
    p = 2 * (1 - norm.cdf(abs(z)))
    h = abs(z) > norm.ppf(1 - alpha / 2)
    trend = 'decreasing' if z < 0 and h else 'increasing' if z > 0 and h else 'no trend'
    return p, h, trend

def __mk_probability(p, s):
    """
    Computes the Mann-Kendall probability.
    """
    C = 1 - p / 2
    Cd = C if s <= 0 else p / 2
    return C, Cd

def __sens_estimator_unequal_spacing(x, t):
    """
    Computes Sen's slope for unequally spaced data using a vectorized approach.
    """
    n = len(x)
    if n < 2:
        return np.array([])

    # Create all pairs of indices
    i, j = np.triu_indices(n, k=1)

    # Calculate differences
    x_diff = x[j] - x[i]
    t_diff = t[j] - t[i]

    # Avoid division by zero
    valid_mask = t_diff != 0

    return x_diff[valid_mask] / t_diff[valid_mask]


def _sens_estimator_censored(x, t, cen_type, lt_mult=0.5, gt_mult=1.1, method='lwp'):
    """
    Computes Sen's slope for censored, unequally spaced data.

    This implements the logic from the LWP-TRENDS R script, which is the
    default behavior. An optional, more statistically robust method is also
    provided.

    Args:
        x (np.array): The data values.
        t (np.array): The timestamps.
        cen_type (np.array): The censor types ('lt', 'gt', 'not').
        lt_mult (float): Multiplier for left-censored data.
        gt_mult (float): Multiplier for right-censored data.
        method (str): The method to use for handling ambiguous slopes.
            - 'lwp' (default): Sets ambiguous slopes to 0, mimicking the
              LWP-TRENDS R script.
            - 'nan': Sets ambiguous slopes to np.nan, which is a more
              statistically neutral approach.

    Returns:
        np.array: An array of calculated slopes.
    """
    n = len(x)
    if n < 2:
        return np.array([])

    # Create all pairs of indices
    i_indices, j_indices = np.triu_indices(n, k=1)

    # 1. Calculate raw differences and slopes (for censor checks)
    x_diff_raw = x[j_indices] - x[i_indices]
    t_diff = t[j_indices] - t[i_indices]

    # Avoid division by zero
    valid_t_mask = t_diff != 0
    x_diff_raw = x_diff_raw[valid_t_mask]
    t_diff = t_diff[valid_t_mask]
    i_indices = i_indices[valid_t_mask]
    j_indices = j_indices[valid_t_mask]

    slopes_raw = x_diff_raw / t_diff

    # 2. Modify values for final slope calculation (as per R script)
    x_mod = x.copy().astype(float)
    x_mod[cen_type == 'lt'] *= lt_mult
    x_mod[cen_type == 'gt'] *= gt_mult
    x_diff_mod = x_mod[j_indices] - x_mod[i_indices]
    slopes_mod = x_diff_mod / t_diff

    # 3. Create censor labels for pairs and apply rules
    cen_type_pairs = cen_type[i_indices] + " " + cen_type[j_indices]
    slopes_final = slopes_mod.copy()

    # Determine the value to assign to ambiguous slopes based on the method
    ambiguous_slope_value = 0 if method == 'lwp' else np.nan

    # Rule 1: No slope between two censored values of the same type
    slopes_final[(cen_type_pairs == 'gt gt') | (cen_type_pairs == 'lt lt')] = ambiguous_slope_value

    # Rules 2 & 3: Ambiguous slopes between left-censored and non-censored
    slopes_final[(slopes_raw > 0) & ((cen_type_pairs == 'not lt') | (cen_type_pairs == 'lt not'))] = ambiguous_slope_value

    # Rules 4 & 5: Ambiguous slopes between right-censored and non-censored
    slopes_final[(slopes_raw < 0) & ((cen_type_pairs == 'not gt') | (cen_type_pairs == 'gt not'))] = ambiguous_slope_value

    return slopes_final

def __confidence_intervals(slopes, var_s, alpha):
    """
    Computes the confidence intervals for Sen's slope.
    """
    # Filter out NaN values, which can occur with the 'nan' method for
    # censored slopes.
    valid_slopes = slopes[~np.isnan(slopes)]
    n = len(valid_slopes)

    if n == 0 or var_s == 0:
        return np.nan, np.nan

    # For a two-sided confidence interval
    Z = norm.ppf(1 - alpha / 2)

    # Ranks of the lower and upper confidence limits (1-based)
    C = Z * np.sqrt(var_s)
    M1 = (n - C) / 2
    M2 = (n + C) / 2

    sorted_slopes = np.sort(valid_slopes)

    # Interpolate to find the values at the fractional ranks
    ranks = np.arange(1, n + 1)
    lower_ci = np.interp(M1, ranks, sorted_slopes)
    upper_ci = np.interp(M2, ranks, sorted_slopes)

    return lower_ci, upper_ci


def _aggregate_censored_median(group, is_datetime):
    """
    Computes a robust median for a group of observations which may contain
    censored data, following the LWP-TRENDS R script logic.
    """
    n = len(group)
    if n == 0:
        return pd.DataFrame()  # Changed from None

    # Compute median value
    median_val = group['value'].median()

    # Determine if median is censored (R logic)
    if not group['censored'].any():
        is_censored = False
        cen_type = 'not'
    else:
        # Get maximum censored value
        max_censored = group.loc[group['censored'], 'value'].max()
        is_censored = median_val <= max_censored

        if is_censored:
            # Safely get the most common censor type
            cen_type_mode = group.loc[group['censored'], 'cen_type'].mode()
            if len(cen_type_mode) == 0:
                # All censored values are NaN, default to 'not'
                cen_type = 'not'
                is_censored = False
            else:
                cen_type = cen_type_mode.iloc[0]
        else:
            cen_type = 'not'

    row_data = {
        'value': median_val,
        'censored': is_censored,
        'cen_type': cen_type,
    }

    # Always aggregate time using the median of the original timestamps
    row_data['t_original'] = group['t_original'].median() if is_datetime else np.median(group['t_original'])
    row_data['t'] = np.median(group['t'])

    return pd.DataFrame([row_data])


def _prepare_data(x, t, hicensor):
    """
    Internal helper to prepare and validate data for trend tests.
    """
    if isinstance(x, pd.DataFrame) and all(col in x.columns for col in ['value', 'censored', 'cen_type']):
        data = x.copy()
    elif hasattr(x, '__iter__') and any(isinstance(i, str) for i in x):
        raise TypeError("Input data `x` contains strings. Please pre-process it with `prepare_censored_data` first.")
    else:
        x_proc, _ = __preprocessing(x)
        data = pd.DataFrame({
            'value': x_proc,
            'censored': np.zeros(len(x_proc), dtype=bool),
            'cen_type': np.full(len(x_proc), 'not', dtype=object)
        })

    t_raw = np.asarray(t)
    is_datetime = _is_datetime_like(t_raw)
    t_numeric, _ = __preprocessing(t_raw)
    data['t_original'] = t_raw
    data['t'] = t_numeric

    # Handle missing values
    mask = ~np.isnan(data['value'])
    data_filtered = data[mask].copy()

    # Apply HiCensor rule if requested
    if hicensor and 'lt' in data_filtered['cen_type'].values:
        max_lt_censor = data_filtered.loc[data_filtered['cen_type'] == 'lt', 'value'].max()
        hi_censor_mask = data_filtered['value'] < max_lt_censor
        data_filtered.loc[hi_censor_mask, 'censored'] = True
        data_filtered.loc[hi_censor_mask, 'cen_type'] = 'lt'
        data_filtered.loc[hi_censor_mask, 'value'] = max_lt_censor

    return data_filtered, is_datetime


def _aggregate_by_group(group, agg_method, is_datetime):
    """
    Aggregates a group of data points using the specified method.
    """
    if len(group) <= 1:
        return group

    if agg_method == 'median':
        if group['censored'].any():
            import warnings
            warnings.warn(
                "The 'median' aggregation method uses a simple heuristic for censored data, "
                "which may not be statistically robust. Consider using 'robust_median' for "
                "more accurate censored data aggregation.", UserWarning)
        median_val = group['value'].median()
        is_censored = median_val <= group[group['censored']]['value'].max() if group['censored'].any() else False

        new_row = {
            'value': median_val,
            't_original': group['t_original'].median() if is_datetime else np.median(group['t_original']),
            't': np.median(group['t']),
            'censored': is_censored,
            'cen_type': group.loc[group['censored'], 'cen_type'].mode()[0] if is_censored else 'not'
        }
        return pd.DataFrame([new_row])
    elif agg_method == 'robust_median':
        return _aggregate_censored_median(group, is_datetime)
    elif agg_method == 'middle':
        t_numeric_group = group['t'].to_numpy()
        closest_idx = np.argmin(np.abs(t_numeric_group - np.mean(t_numeric_group)))
        return group.iloc[[closest_idx]]
    return group
