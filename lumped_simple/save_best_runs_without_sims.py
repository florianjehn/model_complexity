# -*- coding: utf-8 -*-
"""
Created on Aug 04 11:08 2017
@author(s): Florian U. Jehn
"""

import pandas as pd
import matplotlib.pyplot as plt


def read_data(filename):
    """
    Reads in the data from a csv file.

    :param filename:
    :return: pd.dataframe
    """
    results = pd.read_csv(filename, delimiter=",")
    objective_functions = 0
    params = 0
    for col_name in results.columns.values.tolist():
        if "sim" in col_name:
            break
        if "like" in col_name:
            objective_functions += 1
        if "par" in col_name:
            params += 1
    results_no_sims = results.ix[:, 0:objective_functions + params]
    return results_no_sims


def save_best_runs(results_no_sims, org_name):
    """
    Saves the best 20 % of the last 10 % of runs Rope has produced.

    :param results_no_sims: results without the simulation data
    :param org_name: name of the original file
    :return: None
    """

    repetitions = results_no_sims.shape[0]
    # Calculate how large the last 10 % of all runs are
    last_ten_percent = int(repetitions * 0.1)
    # calculate how large the best 20 % of the last 10 % are.
    best_20_percent = int(last_ten_percent * 0.2)

    # Get the last 10 % of the dataframe
    last_ten_percent_dataframe = results_no_sims.tail(last_ten_percent)

    # Get the best 20 % of the last 10 %
    sorted_dataframe = last_ten_percent_dataframe.sort_values(by="like2")
    best_runs = sorted_dataframe.tail(best_20_percent)

    best_runs.to_csv(org_name[:-4] + "_best_runs.csv", index=False)

names = ["simple_lumped_2_subsets.csv", "simple_lumped_3_subsets.csv",
         "simple_lumped_4_subsets.csv", "simple_lumped.csv"]


if __name__ == '__main__':
    for name in names:

        results_no_sims = read_data(name)
        save_best_runs(results_no_sims, name)



