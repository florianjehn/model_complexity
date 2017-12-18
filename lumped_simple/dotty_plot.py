# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:27:20 2016

@author: Florian Jehn

Produces a simple dotty plot for a .csv output file of the tutorial template
for CMF
Plots each Parameter against the effiency and writes results to .png if wanted
name = name of the .csv

Last edited: 18.05.2016

"""

import pandas as pd
import matplotlib.pyplot as plt


def dotty_plot(name, amount_params):
    """
    reads in a .csv and plots it as a dotty plot
    """
    csv = pd.read_csv(name)
    col_names_param = csv.columns
    # go through all parameters
    for param in col_names_param[:amount_params + 1]:
        # Skip efficiency to not plot it against itself
        if param != col_names_param[0]:

            plt.scatter(csv[param], csv[col_names_param[0]], alpha=0.2)
            plt.xlabel(param)
            plt.ylabel(col_names_param[0])
            plt.savefig(col_names_param[0] + "_" + param+"_dotty.png",
                        dpi=250, bbox_inches="tight")
            plt.close()

    print("done")


# Enter your file name here
dotty_plot("simple_lumped_500_lhs_short.csv", 10)


