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
        t_numeric, _ = __preprocessing(t)
        seasons = (np.floor(t_numeric - 1) % period).astype(int)

    df = pd.DataFrame({'Value': x, 'Season': seasons})

    plt.figure(figsize=(10, 6))
    sns.boxplot(x='Season', y='Value', data=df)
    plt.title('Distribution of Values Across Seasons')
    plt.xlabel('Season')
    plt.ylabel('Value')

    plt.savefig(save_path)
    plt.close()

    return save_path
