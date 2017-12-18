#/usr/bin/env python
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
    results = read_data(filename)
    count_NS_over_thresh(results, threshold)
    results["like1"].hist(bins=500, normed=1)
    plt.axis([0, 1, 0, 10])
    #plt.show()
    plt.savefig("histogram_lumped_conmplex.png")


make_hist("complex_lumped_1000_rope.csv", 0.6)

