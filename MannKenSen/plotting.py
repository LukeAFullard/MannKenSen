"""
This script provides plotting utilities for the MannKenSen package.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from ._utils import __preprocessing, _get_season_func, _is_datetime_like, _get_cycle_identifier

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


def plot_inspection_data(data, save_path, value_col, time_col, time_increment, increment_map):
    """
    Creates and saves a 2x2 grid of data inspection plots.
    """
    # 1. Validate data
    if 'censored' not in data.columns:
        raise ValueError(
            "Input data does not appear to be prepared for censored analysis. "
            "Please run `prepare_censored_data` first to add 'censored' and 'cen_type' columns."
        )

    fig, axes = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f"Data Inspection (Time Increment: {time_increment})", fontsize=16)

    # Plot 1: Time series plot (consistent with plot_trend)
    ax1 = axes[0, 0]
    is_datetime = _is_datetime_like(data[time_col].values)
    x_axis = pd.to_datetime(data[time_col]) if is_datetime else data[time_col]

    censored_mask = data['censored']
    ax1.scatter(x_axis[~censored_mask], data.loc[~censored_mask, value_col],
                color='blue', label='Non-censored', marker='o', alpha=0.7)
    if censored_mask.any():
        ax1.scatter(x_axis[censored_mask], data.loc[censored_mask, value_col],
                    color='red', label='Censored', marker='x')
    ax1.set_title('Time Series Plot')
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Value')
    ax1.legend()
    ax1.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.setp(ax1.get_xticklabels(), rotation=30, ha='right')

    # Get season and cycle columns for matrix plots
    season_col = increment_map.get(time_increment)

    if season_col is None or time_increment == 'none' or not is_datetime:
        plot_title = "Cannot generate matrix plots.\n"
        if not is_datetime:
            plot_title += "Matrix plots require datetime-like time column."
        else:
            plot_title += "No suitable time increment found."

        for ax in [axes[0, 1], axes[1, 0], axes[1, 1]]:
            ax.text(0.5, 0.5, plot_title, ha='center', va='center', fontsize=12)
            ax.set_xticks([])
            ax.set_yticks([])
    else:
        plot_df = data.copy()
        # Dynamically determine the cycle (e.g., year for monthly, week for daily)
        plot_df['cycle'] = _get_cycle_identifier(plot_df[time_col], time_increment)

        # Plot 2: Value matrix
        ax2 = axes[0, 1]
        try:
            data_pivot = plot_df.pivot_table(
                values=value_col, index='cycle', columns=season_col, aggfunc='median')
            sns.heatmap(data_pivot, cmap='viridis', ax=ax2, cbar_kws={'label': 'Median Value'})
            ax2.set_title(f'Median Value ({season_col.capitalize()} vs. Cycle)')
            ax2.set_xlabel(season_col.capitalize())
            ax2.set_ylabel('Cycle')
        except Exception as e:
            ax2.text(0.5, 0.5, f"Could not generate value matrix:\n{e}", ha='center', va='center')

        # Plot 3: Censoring matrix
        ax3 = axes[1, 0]
        try:
            cens_pivot = plot_df.pivot_table(
                values='censored', index='cycle', columns=season_col,
                aggfunc='any', fill_value=False)
            sns.heatmap(cens_pivot.astype(int), cmap='coolwarm', ax=ax3,
                       cbar_kws={'label': 'Censored (1=True)'}, vmin=0, vmax=1)
            ax3.set_title('Censoring Status')
            ax3.set_xlabel(season_col.capitalize())
            ax3.set_ylabel('Cycle')
        except Exception as e:
            ax3.text(0.5, 0.5, f"Could not generate censoring matrix:\n{e}", ha='center', va='center')

        # Plot 4: Sample count matrix
        ax4 = axes[1, 1]
        try:
            count_pivot = plot_df.pivot_table(
                values=value_col, index='cycle', columns=season_col,
                aggfunc='count', fill_value=0)
            sns.heatmap(count_pivot, cmap='Blues', ax=ax4,
                       cbar_kws={'label': 'Sample Count'}, annot=True, fmt='g')
            ax4.set_title('Sample Count')
            ax4.set_xlabel(season_col.capitalize())
            ax4.set_ylabel('Cycle')
        except Exception as e:
            ax4.text(0.5, 0.5, f"Could not generate count matrix:\n{e}", ha='center', va='center')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()

    return save_path

def plot_trend(data, results, save_path, alpha):
    """
    Generates and saves a plot of the data with the calculated trend line.

    Input:
        data (pd.DataFrame): The DataFrame containing the data, including 'value',
                             'censored', 't', and optionally 't_original'.
        results (namedtuple): The results from original_test or seasonal_test.
        save_path (str): The file path to save the plot.
        alpha (float): The significance level for the confidence intervals.
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

        # Confidence interval lines, pivoted around the median data point
        ymed = np.nanmedian(data['value'])
        tmed = np.nanmedian(data['t'])

        # Correctly calculate intercepts so lines pass through (tmed, ymed)
        intercept_lower = ymed - results.lower_ci * tmed
        intercept_upper = ymed - results.upper_ci * tmed

        lower_line = results.lower_ci * np.array([t_min, t_max]) + intercept_lower
        upper_line = results.upper_ci * np.array([t_min, t_max]) + intercept_upper

        x_line = pd.to_datetime([t_min, t_max], unit='s') if is_datetime else [t_min, t_max]

        plt.plot(x_line, trend_line, color='black', linestyle='--', label="Sen's Slope")
        ci_label = f'{int((1 - alpha) * 100)}% CI'
        plt.fill_between(x_line, lower_line, upper_line, color='gray', alpha=0.3, label=ci_label)

    # Add statistics text box
    stats_text = (f"Trend: {results.trend}\n"
                  f"Tau: {results.Tau:.4f}\n"
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
