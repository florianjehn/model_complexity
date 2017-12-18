# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 11:03:26 2016

@author: Florian
"""
import csv
def spotpy_csv_shrinker(in_name, out_name, thresh_NS):
    """
    Goes through a csv file produced by Spotpy and creates
    a new csv file with all entries above eff_thresh
    """
    # create an empy list to temp save the data
    temp = []
    # read through the file
    with open(in_name, newline="") as csv_in:
        reader = csv.reader(csv_in, delimiter = ",")
        for row in reader:
            # use try to avoid problems with "nan" and "efficieny"
            try:
                # only write those values to the file that are above thresh
                if float(row[0]) > thresh_NS:
                    temp.append(row)
            except ValueError:
                # Add the first row if it pops upp
                if row[0] == "like1":
                    temp.append(row)
                    
    # write the new csv file
    with open(out_name, "w", newline="") as csv_out:
        writer = csv.writer(csv_out, delimiter=",", quoting=csv.QUOTE_MINIMAL)
        for row in temp:
            writer.writerow(row)
                    

spotpy_csv_shrinker("complex_lumped_500_lhs.csv",
                    "complex_lumped_500_lhs_short.csv", 0.0)
