"""
This script provides a modified version of the Mann-Kendall test
and Sen's slope estimator to handle unequally spaced time series data.
"""

from __future__ import division
import datetime
import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata


def prepare_censored_data(x):
    """
    Pre-processes a 1D array-like object containing censored data.

    This function takes an array that may contain a mix of numeric and
    string values (e.g., 10, '<5', '>20') and converts it into a structured
    pandas DataFrame with separate columns for the numeric value, a boolean
    censored flag, and the type of censoring.

    Args:
        x (array-like): A 1D array or list containing the data.

    Returns:
        pandas.DataFrame: A DataFrame with three columns:
            - 'value': The numeric value of the data point.
            - 'censored': A boolean flag, True if the data point was censored.
            - 'cen_type': A string indicating the type of censoring
                          ('lt', 'gt', or 'not').
    """
    values = []
    censored_flags = []
    cen_types = []

    if not hasattr(x, '__iter__') or isinstance(x, str):
        raise TypeError("Input data must be an iterable (e.g., list, numpy array).")

    for item in x:
        if isinstance(item, str):
            item_stripped = item.strip()
            if item_stripped.startswith('<'):
                try:
                    values.append(float(item_stripped[1:]))
                    censored_flags.append(True)
                    cen_types.append('lt')
                except (ValueError, IndexError):
                    raise ValueError(f"Invalid left-censored value format: '{item}'")
            elif item_stripped.startswith('>'):
                try:
                    values.append(float(item_stripped[1:]))
                    censored_flags.append(True)
                    cen_types.append('gt')
                except (ValueError, IndexError):
                    raise ValueError(f"Invalid right-censored value format: '{item}'")
            else:
                try:
                    values.append(float(item_stripped))
                    censored_flags.append(False)
                    cen_types.append('not')
                except ValueError:
                    raise ValueError(f"Could not convert string '{item}' to a float.")
        else:
            try:
                values.append(float(item))
                censored_flags.append(False)
                cen_types.append('not')
            except (ValueError, TypeError):
                 raise ValueError(f"Could not convert non-string value '{item}' to a float.")

    return pd.DataFrame({
        'value': np.array(values, dtype=float),
        'censored': np.array(censored_flags, dtype=bool),
        'cen_type': np.array(cen_types, dtype=object)
    })


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

    if season_type in ['month', 'quarter', 'year', 'day_of_year', 'week_of_year']:
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

def _rle_lengths(a):
    """Calculates the lengths of runs of equal values in an array."""
    if len(a) == 0:
        return np.array([], dtype=int)
    # Pad the array to detect changes at the start and end
    # and find the indices where the array changes
    pad = np.array([False])
    idx = np.flatnonzero(np.concatenate((pad, a[1:] != a[:-1], pad)))
    return np.diff(idx)


def _mk_score_and_var_censored(x, t, censored, cen_type):
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
        # Add a small amount to break ties, treat as uncensored
        max_gt_val = x_mod[gt_mask].max() + 0.1
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

    itot = np.sum(np.triu(signyx * (1 - xplus), k=1))
    kenS = itot
    D = n * (n - 1) / 2.0

    # 5. Variance Calculation (adapted from NADA::cenken)
    varS = n * (n - 1) * (2 * n + 5) / 18.0
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


def _sens_estimator_censored(x, t, cen_type):
    """
    Computes Sen's slope for censored, unequally spaced data.
    This is a Python translation of the GetInterObservationSlopes logic
    from the LWP-TRENDS R script.
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
    x_mod = x.copy()
    x_mod[cen_type == 'lt'] *= 0.5
    x_mod[cen_type == 'gt'] *= 1.1
    x_diff_mod = x_mod[j_indices] - x_mod[i_indices]
    slopes_mod = x_diff_mod / t_diff

    # 3. Create censor labels for pairs and apply rules
    cen_type_pairs = cen_type[i_indices] + " " + cen_type[j_indices]
    slopes_final = slopes_mod.copy()

    # Rule 1: No slope between two censored values of the same type
    slopes_final[(cen_type_pairs == 'gt gt') | (cen_type_pairs == 'lt lt')] = 0

    # Rules 2 & 3: Ambiguous slopes between left-censored and non-censored
    slopes_final[(slopes_raw > 0) & ((cen_type_pairs == 'not lt') | (cen_type_pairs == 'lt not'))] = 0

    # Rules 4 & 5: Ambiguous slopes between right-censored and non-censored
    slopes_final[(slopes_raw < 0) & ((cen_type_pairs == 'not gt') | (cen_type_pairs == 'gt not'))] = 0

    return slopes_final

def __confidence_intervals(slopes, var_s, alpha):
    """
    Computes the confidence intervals for Sen's slope.
    """
    n = len(slopes)
    if n == 0 or var_s == 0:
        return np.nan, np.nan

    # For a two-sided confidence interval
    Z = norm.ppf(1 - alpha / 2)

    # Ranks of the lower and upper confidence limits (1-based)
    C = Z * np.sqrt(var_s)
    M1 = (n - C) / 2
    M2 = (n + C) / 2

    sorted_slopes = np.sort(slopes)

    # Interpolate to find the values at the fractional ranks
    ranks = np.arange(1, n + 1)
    lower_ci = np.interp(M1, ranks, sorted_slopes)
    upper_ci = np.interp(M2, ranks, sorted_slopes)

    return lower_ci, upper_ci
