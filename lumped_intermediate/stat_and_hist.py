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


def count_NS_over_thresh(results, threshold):
    """
    :return: None
    """
    over_tresh = 0
    for value in results["like1"]:
        if value > threshold:
            over_tresh += 1

    print("The results have {} entries with an NS over {}".format(over_tresh,
                                                                  threshold))


def make_hist(filename, threshold):
    """
    Makes a histogram of the results.

    :return: None
    """
    # First make a histogram of all results
    results = read_data(filename)
    count_NS_over_thresh(results, threshold)
    results["like1"].hist(bins=500)
    plt.xlim([0, 1])
    plt.savefig("histogram_intermediate_lumped_rope_70000_ndir_10000_all_runs"
                ".png")
    plt.close()
    # Then make only a histogram of the best 1000 runs of the last 5000 runs
    last_5000 = results.tail(5000)
    sorted_5000 = last_5000.sort_values(by="like1")
    best_1000 = sorted_5000.tail(1000)
    print("Test")
    best_1000["like1"].hist(bins=20)
    plt.xlim([0, 1])
    plt.savefig(
        "histogram_intermediate_lumped_rope_70000_ndir_10000_best_1000"
        ".png")
    plt.close()


make_hist("intermediate_lumped_rope_70000_first_50000_ndir_10000.csv", 0.65)
