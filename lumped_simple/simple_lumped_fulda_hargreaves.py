#/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Jul 24 10:24 2017
@author(s): Florian U. Jehn
"""
import cmf
import spotpy
from spotpy.parameter import Uniform as param
import os
import datetime
import sys
import numpy as np
from dateutil.relativedelta import relativedelta


class SimpleLumped(object):
    """
    Class which contains the complete model, readeable for Spotpy
    """
    def __init__(self, begin, end):
        """Initializes the model and build the core setup"""
        # tr_soil_GW = Residence time of the water in the soil to the GW
        self.params = [param("tr_soil_out_days", 0., 200.),
                       # tr_GW_out = Residence time in the groundwater to
                       #  the outlet
                       param('V0_soil_mm', 0., 300.),
                       # beta_soil_out = exponent that changes the form of the
                       # flux from the soil to the outlet
                       param("beta_soil_out", 0.3, 8.0),
                       # ETV1 = the Volume that defines that the evaporation
                       # is lowered because of not enough water in the soil
                       param('ETV1_mm', 0., 300.),
                       # fETV0 = factor the ET is multiplied by, when water is
                       #  low
                       param('fETV0', 0., 0.9),
                       # Rate of snow melt
                       param('meltrate_mm_degC_day', 0.01, 15.),
                       # Snow_melt_temp = Temperature at which the snow melts
                       # (needed because of averaged temp
                       param('snow_melt_temp_degC', -3.0, 3.0)
                       ]
        self.begin = begin
        self.end = end

        # load the weather data and discharge data
        prec, temp, temp_min, temp_max, Q = self.loadPETQ()
        self.Q = Q

        # use only one core (faster when model is small
        cmf.set_parallel_threads(1)

        # Generate a cmf project with one cell for a lumped model
        self.project = cmf.project()
        p = self.project

        # Create a cell in the project
        c = p.NewCell(0, 0, 0, 562 * 1e6)

        # Add snow storage
        c.add_storage("Snow", "S")
        cmf.Snowfall(c.snow, c)

        # Add the soil and groundwater layers
        soil = c.add_layer(2.0)

        # Give the storages a initial volume
        c.layers[0].volume = 15

        # Install a calculation of the evaporation
        cmf.HargreaveET(soil, c.transpiration)

        # Create an outlet
        self.outlet = p.NewOutlet("outlet", 10, 0, 0)

        # Create the meteo stations
        self.make_stations(prec, temp, temp_min, temp_max)
        self.project = p


    def set_parameters(self,
                       tr_soil_out_days,

                       V0_soil_mm,

                       beta_soil_out,

                       ETV1_mm,
                       fETV0,

                       meltrate_mm_degC_day,
                       snow_melt_temp_degC,
                       ):
        """
        Creates all connections with the parameter values produced by the
        sampling algorithm.
        """
        # Get all definition from the init method
        p = self.project
        c = p[0]
        outlet = self.outlet
        soil = c.layers[0]

        ETV1_m3 = (ETV1_mm / 1000) * c.area
        V0_soil_m3 = (V0_soil_mm / 1000) * c.area

        # Adjustment of the ET
        c.set_uptakestress(cmf.VolumeStress(ETV1_m3, ETV1_m3 * fETV0))

        # Flux from soil to outlet
        cmf.kinematic_wave(soil, outlet, tr_soil_out_days/V0_soil_m3,
                           V0=V0_soil_m3,
                           exponent=beta_soil_out)

        # # Set parameters of the snow calculations
        cmf.Weather.set_snow_threshold(snow_melt_temp_degC)
        cmf.SimpleTindexSnowMelt(c.snow, soil, c, rate=meltrate_mm_degC_day)

    def loadPETQ(self):
        """
        Loads climata and discharge data from the corresponding files fnQ,
        fnT and fnP
        """
        # Fixed model starting point
        begin = self.begin - relativedelta(years=1)
        step = datetime.timedelta(days=1)
        # empty time series
        prec = cmf.timeseries(begin, step)
        prec.extend(float(Pstr.strip("\n")) for Pstr in open(fnP))

        discharge = cmf.timeseries(begin, step)
        discharge.extend(float(Qstr.strip("\n")) for Qstr in open(fnQ))
        # Convert m3/s to mm/day
        area_catchment = 562.41
        # 86400 = seconds per day
        #discharge *= 86400 * 1e3 / (area_catchment * 1e6)
        temp = cmf.timeseries(begin, step)
        temp_min = cmf.timeseries(begin, step)
        temp_max = cmf.timeseries(begin, step)

        # Go through all lines in the file
        for line in open(fnT):
            columns = line.strip("\n").split('\t')
            if len(columns) == 3:
                temp_max.add(float(columns[0]))
                temp_min.add(float(columns[1]))
                temp.add(float(columns[2]))

        return prec, temp, temp_min, temp_max, discharge

    def make_stations(self, prec, temp, temp_min, temp_max):
        """
        Creates the cmf weather stations
        """
        rainstation = self.project.rainfall_stations.add("Rainfall Station",
                                                         prec, (0, 0, 0))
        self.project.use_nearest_rainfall()

        # Temperature data
        meteo = self.project.meteo_stations.add_station('Grebenau avg', (0,
                                                                         0, 0))
        meteo.T = temp
        meteo.Tmin = temp_min
        meteo.Tmax = temp_max
        self.project.use_nearest_meteo()
        return rainstation

    def run_model(self):
        """
        Starts the model. Used by spotpy
        """
      #  start = datetime.datetime.now()
     #   print("\nStarting new run\n")
        try:
            # Create a solver for differential equations
            solver = cmf.CVodeIntegrator(self.project, 1e-8)

            # New time series for model results
            resQ = cmf.timeseries(self.begin, cmf.day)
            # starts the solver and calculates the daily time steps
            end = self.end
          #  year_model = self.begin
         #   temp = start
            for t in solver.run(self.project.meteo_stations[0].T.begin,
                                end,
                                cmf.day):
                # if t > year_model + relativedelta(years=1):
                #     year_model += relativedelta(years=1)
                #     duration_year = temp - datetime.datetime.now()
                #     temp = datetime.datetime.now()
                #     print("Year took {} seconds to complete".format(-round(
                #         duration_year.total_seconds(), 3)))

                # Fill the results (first year is included but not used to
                # calculate the NS)
                if t >= self.begin:
                    resQ.add(self.outlet.waterbalance(t))

          #  end_run = datetime.datetime.now()
            # print("Current run took {} seconds to finish".format(-round(
            #     (start - end_run).total_seconds(), 2)))
            return resQ

        # Return an nan - array when a runtime error occurs
        except RuntimeError:
            # end_run = datetime.datetime.now()
            # print("Current run took {} seconds to finish".format(-round(
            #     (start - end_run).total_seconds(), 2)))
            return np.array(self.Q[
                            self.begin:self.end + datetime.timedelta(
                                days=1)]) * np.nan

    def simulation(self, vector):
        """
        SpotPy expects a method simulation. This methods calls set_parameters
        and run_models, so SpotPy is satisfied
        """
        paramdict = dict((pp.name, v) for pp, v in zip(self.params, vector))
        self.set_parameters(**paramdict)
        try:
            resQ = self.run_model()
        except KeyboardInterrupt:
            resQ = np.array(self.Q[
                            self.begin:self.end + datetime.timedelta(
                                days=1)])*np.nan
        resQ = resQ / 86400
        return np.array(resQ)

    def evaluation(self):
        """
        For Spotpy
        """
        return np.array(
            self.Q[self.begin:self.end + datetime.timedelta(days=1)])

    def parameters(self):
        """
        For Spotpy
        """
        return spotpy.parameter.generate(self.params)

    def objectivefunction(self, simulation, evaluation):
        """
        For Spotpy
        """
        # Slice the arrays to only use the days of the calibration period
        # for objective function
        # 1827 = 2 * 366 + 3 * 365
        evaluation_calibration = evaluation[:1827]
        evaluation_validation = evaluation[1827:]
        simulation_calibration = simulation[:1827]
        simulation_validation = simulation[1827:]
        ns_calibration = spotpy.objectivefunctions.nashsutcliffe(
                                                        evaluation_calibration,
                                                        simulation_calibration)
        ns_validation = spotpy.objectivefunctions.nashsutcliffe(
                                                        evaluation_validation,
                                                        simulation_validation)
        return [ns_calibration, ns_validation]


if __name__ == '__main__':

    # 1979 is spin up
    # Whole modelling period
    # 1980 till 1984 is used for calibration
    # 1985 till 1989 is used for validation
    begin = 1980
    end = 1989

    # For validation data from Luterbacher is used
    prefix = "simple_lumped"

    # Number of runs
    runs = 15

    # File names of the forcing data
    fnQ = "Q_Kammerzell_1979_1999.txt"
    fnT = "T_kammerzell_1979_1999_max_min_avg.txt"
    fnP = "P_Krigavg_kammerzell_1979_1999.txt"

    # import algorithm
    from spotpy.algorithms import lhs as Sampler
    #Sampler = rope.rope

    # Find out if the model should run parallel (for supercomputer)
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'

    # Create the model
    model = SimpleLumped(datetime.datetime(begin, 1, 1), datetime.datetime(
        end, 12, 31))

    # If there is an command line argument, take its value for the amount of
    #  runs
    if len(sys.argv) > 1:
        runs = int(sys.argv[1])

    # run the model
    if runs:
        sampler = Sampler(model, parallel=parallel, dbname="simple_lumped",
                          dbformat="csv", save_sim=True,save_threshold=0.0)
        sampler.sample(runs)#, subsets=30)


