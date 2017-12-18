# -*- coding: utf-8 -*-
"""
Created on Sep 19 09:48 2017
@author(s): Florian U. Jehn
"""
import pandas as pd
import matplotlib.pyplot as plt


names = ["simple_lumped_2_subsets_best_runs.csv",
         "simple_lumped_3_subsets_best_runs.csv",
         "simple_lumped_4_subsets_best_runs.csv",
         "simple_lumped_new_rope_best_runs.csv"]

for name in names:
    dt = pd.read_csv(name)
    max = dt["like2"].max()
    min = dt["like2"].min()
    mean = dt["like2"].mean()
    std = dt["like2"].std()
    amount = dt["like2"].size
    print(max, min, mean, std)

    dt["like2"].hist(bins=20)
    plt.xlabel("Nash-Sutcliffe")
    plt.ylabel("Count")
    plt.xlim(0.67, 0.73)
    plt.ylim(0, 300)
    plt.title("\n\n" + name + "\nmax: {} min: {} mean: {} \nstd: {} amount: {}".format(
        str(round(max, 4)),
        str(round(min, 4)),
        str(round(mean, 4)),
        str(round(std, 4)),
        str(amount)
    ))
    plt.savefig(name + "histogram.png")
    plt.close()







