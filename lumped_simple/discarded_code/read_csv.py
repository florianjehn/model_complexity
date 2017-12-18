# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 10:38:38 2017

@author: Florian Jehn
"""
import csv

def read_csv_for_objective_funcs(name):
    """
    reads a from a shortened parameter list csv and transforms it to a 
    list of lists of floats 
    the outer list is just a container for all behavioural runs
    each inner list equals one behavioural run of the model
    the floats being the values for the three objective functions
    the parameterer values after them are ignored    
    """
    model = []
    with open(name, newline="") as csvfile:
        reader = csv.reader(csvfile, delimiter=",")
        for row in reader:
            try:
                logNS = float(row[0])
                pbias = float(row[1])
                rsr = float(row[2])
                model.append([logNS, pbias, rsr])
            except ValueError:
                pass
    return model
            


            
            
def read_csv_for_sim_dis(name, num_params, num_obj_func, day_off = False):
    """
    reads from a shortened simulation list csv and transforms it into a list 
    of lists of floats
    the outer list is just a container for all behavioural runs
    , the inner list are the single behavioural runs and the floats are
    the values of the days
    in: 
    - name = name of the csv
    - num_params = number of parameters of the model
    - num_obj_func = number of objective functions saved in the csv
    - day_off = indicates if the csv file is saved with the day off method
    """
    model = []
    day_off_increment = 0
    if day_off:
        day_off_increment = 1
    
    with open(name, newline = "") as csvfile:
        reader = csv.reader(csvfile, delimiter = ",")
        # skip the first row with the names
        next(reader, None)
        i = 0
        for row in reader:
            values = row[(day_off_increment + num_params + num_obj_func):]
            model.append([])
            for val in values:
                model[i].append(float(val))
            i += 1
                
    return model
            
            
if __name__ == "__main__":
 #   test = read_csv_for_objective_funcs("model1-parameters_short.csv")
    test2 = read_csv_for_sim_dis("model1-simulation_short.csv", 19, 3, 
                                 day_off = True)
    