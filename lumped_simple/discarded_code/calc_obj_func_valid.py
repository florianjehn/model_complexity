# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 15:30:54 2017

@author: Florian Jehn
"""
import matplotlib.pyplot as plt
import numpy as np
import read_csv
import csv
import spotpy


def read_data(csv):
    """
    loads all the simulation files of all models and returns them as a
    dict of lists of float
    
    models = dict with number of the model as key and number of the parameters 
    as value
    """
    # load all the simulation files
    models_sims = {}
        
    models_sims["Lumped_simple"] = read_csv.read_csv_for_sim_dis(csv, 10, 1)
    return models_sims
    

    
    
def write_csv(models_obj_funcs, calib_valid):
    """
    writes the calculated objective functions into csv files (one for each 
    model), in a way that can be understood by the radar plot function.
    
    expects a dict with the number of the model as key and lists of the 
    objective functions as value (for all behavioural runs seperately)   
    """
    for key in models_obj_funcs.keys():
        filename = "model_"+str(key)+"_obj_funcs_"+calib_valid+".csv"
        with open(filename, "w", newline ='') as csvfile:
            writer = csv.writer(csvfile, delimiter = ",", 
                                quoting = csv.QUOTE_MINIMAL)
            writer.writerow(["logNS", "pbias", "rsr"])
            for run in models_obj_funcs[key]:
                writer.writerow(run)
        csvfile.close()
    
    
    
    
def calc_obj_funcs(models_sims, dis_obs, shift_one_day = False):
    """
    calculates the objective functions logNS, pbias and RSR for all behavioural
    runs of all behavioural models
    it expects a dict eith the number of the model as key to lists of the 
    behavioural runs
    returns a dict with the number of the model as key and lists of the 
    objective functions as value (for all behavioural runs seperately)
    """
        
    dis_obs = np.array(dis_obs)
    if shift_one_day:
        dis_obs = dis_obs[1:]

    models_obj_funcs = {}
    for key in models_sims.keys():
        models_obj_funcs[key] = []
        for run in models_sims[key]:
            if shift_one_day:
                run = run[:-1]
            print(len(dis_obs), len(run))
            logNS = spotpy.objectivefunctions.lognashsutcliffe(dis_obs, run)
            pbias = spotpy.objectivefunctions.pbias(dis_obs, run)
            rmse = spotpy.objectivefunctions.rmse(dis_obs, run)
            std = dis_obs.std()
            rsr = rmse/std
            models_obj_funcs[key].append([logNS,pbias,rsr])

            
    return models_obj_funcs

    
    

discharge = np.loadtxt("Q_Kammerzell_1979_1999.txt")
# changes m³/s to mm/day
discharge *= 86400 * 1e3 / (2976.41 * 1e6)

#split into validation and calibration timeperiod
dis_calib = discharge[365:]


# split the simulated timeseries in calibration and validation period
models_sims = read_data("simple_lumped.csv")

models_sims_calib = {}

for key in models_sims.keys():
    models_sims_calib[key] = []
    for run in models_sims[key]:
        models_sims_calib[key].append(run[:2192])

models_obj_funcs_calib = calc_obj_funcs(models_sims_calib, dis_calib)
write_csv(models_obj_funcs_calib, "calibr")
    

