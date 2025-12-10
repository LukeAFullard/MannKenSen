"""
This script provides plotting utilities for the MannKenSen package.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ._utils import __preprocessing, _get_season_func, _is_datetime_like

def plot_seasonal_distribution(x_old, t_old, period=12, season_type='month', save_path='seasonal_distribution.png'):
    """
    Generates and saves a box plot to visualize the distribution of values
    across different seasons.

    Input:
        x_old: a vector of data
        t_old: a vector of timestamps
        period: seasonal cycle (default 12)
        season_type: For datetime inputs, specifies the type of seasonality.
        save_path: The file path to save the plot.
    Output:
        The file path where the plot was saved.
    """
    x_raw = np.asarray(x_old)
    t_raw = np.asarray(t_old)

    is_datetime = _is_datetime_like(t_raw)

    if is_datetime:
        season_func = _get_season_func(season_type, period)

    mask = ~np.isnan(x_raw)
    x, t = x_raw[mask], t_raw[mask]

    if len(x) < 2:
        print("Not enough data to generate a plot.")
        return None

    if is_datetime:
        seasons = season_func(pd.to_datetime(t))
    else:
        t_numeric = np.asarray(t, dtype=np.float64)
        # Normalize to start from 0 for consistent seasonal calculation
        t_normalized = t_numeric - t_numeric[0]
        seasons = (np.floor(t_normalized) % period).astype(int)

    df = pd.DataFrame({'Value': x, 'Season': seasons})

    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Season', y='Value', data=df)
    plt.title('Distribution of Values Across Seasons')
    plt.xlabel('Season')
    plt.ylabel('Value')

    plt.savefig(save_path)
    plt.close()

    return save_path

def plot_trend(data, results, save_path):
    """
    Generates and saves a plot of the data with the calculated trend line.

    Input:
        data (pd.DataFrame): The DataFrame containing the data, including 'value',
                             'censored', 't', and optionally 't_original'.
        results (namedtuple): The results from original_test or seasonal_test.
        save_path (str): The file path to save the plot.
    """
    if save_path is None:
        return

    plt.figure(figsize=(10, 6))

    # Determine x-axis values (datetime or numeric)
    is_datetime = 't_original' in data.columns and _is_datetime_like(data['t_original'].values)
    x_axis = pd.to_datetime(data['t_original']) if is_datetime else data['t']

    # Scatter plot for censored and non-censored data
    censored_data = data[data['censored']]
    non_censored_data = data[~data['censored']]

    plt.scatter(x_axis[non_censored_data.index], non_censored_data['value'],
                color='blue', label='Non-censored', marker='o')
    plt.scatter(x_axis[censored_data.index], censored_data['value'],
                color='red', label='Censored', marker='x')

    # Trend line and confidence intervals
    if pd.notna(results.slope):
        t_numeric = data['t'].values
        t_min, t_max = t_numeric.min(), t_numeric.max()

        # Trend line
        trend_line = results.slope * np.array([t_min, t_max]) + results.intercept

        # Confidence interval lines
        ymed = data['value'].median()
        tmed = data['t'].median()
        intercept_lower = ymed - results.lower_ci * tmed
        intercept_upper = ymed - results.upper_ci * tmed

        lower_line = results.lower_ci * np.array([t_min, t_max]) + intercept_lower
        upper_line = results.upper_ci * np.array([t_min, t_max]) + intercept_upper

        x_line = pd.to_datetime([t_min, t_max], unit='s') if is_datetime else [t_min, t_max]

        plt.plot(x_line, trend_line, color='black', linestyle='--', label="Sen's Slope")
        plt.fill_between(x_line, lower_line, upper_line, color='gray', alpha=0.3, label='95% CI')

    # Add statistics text box
    stats_text = (f"Trend: {results.trend}\n"
                  f"Slope: {results.slope:.4f}\n"
                  f"P-value: {results.p:.4f}")
    plt.gca().text(0.05, 0.95, stats_text, transform=plt.gca().transAxes,
                   fontsize=10, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))

    plt.title('Trend Analysis')
    plt.xlabel('Time')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)

    plt.savefig(save_path)
    plt.close()
