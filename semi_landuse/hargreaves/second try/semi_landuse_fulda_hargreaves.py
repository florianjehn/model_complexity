# -*- coding: utf-8 -*-
"""
Created on Nov 29 13:46 2017
@author(s): Florian U. Jehn
"""
from cell_template import CellTemplate
import cmf
import datetime
import os
import numpy as np
import spotpy
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt


class SemiDisLanduse:
    def __init__(self, begin: datetime.datetime, end: datetime.datetime,
                 subcatchment_names):
        """

        :param begin:
        :param end:
        """
        project = cmf.project()
        # Add outlet
        self.outlet = project.NewOutlet("Outlet", 50, 0, 0)
        self.project = project
        self.begin = begin
        self.end = end
        self.subcatchment_names = subcatchment_names
        self.dis_eval, self.subcatchments = self.load_data()
        self.params = self.create_params()
        self.cell_list = self.create_cells()
        cmf.set_parallel_threads(1)

    def create_cells(self):
        """
        Creates a cell for every subcatchment and stores them in a list.

        :return: cell_list: List of all cells, so they are more easily
        callable
        """
        cell_list = []
        for sub in self.subcatchments:
            new_cell = CellTemplate(self.project, self.subcatchments[sub],
                                    sub, self.outlet)
            cell_list.append(new_cell)
        return cell_list

    def set_parameters(self, params):
        """
        Initializes all cells with the parameters provided

        :return: None
        """
        for cell in self.cell_list:
            cell.set_parameters(params)

    def load_data(self):
        """

        :return:
        """
        # Size of all subcatchments
        sizes = {"crops": 147.107, "grass": 151.110, "wood": 215.157,
                 "rest": 49.036}

        # Average height of all subcatchments
        heights = {"crops": 363.666, "grass": 439.544, "wood": 478.530,
                   "rest": 318.808}

        # Different input data types (except discharge)
        input_data = ["T_avg", "T_min", "T_max", "prec"]

        # Read in the data for all subcatchments separately
        subcatchments = {}
        for sub in self.subcatchment_names:
            subcatchments[sub] = {"size": sizes[sub]}
            subcatchments[sub]["height"] = heights[sub]
            subcatchments[sub]["data"] = {}

            for data_type in input_data:
                name = data_type + "_kaemmerzell_" + sub + "_79_89.txt"
                timeseries = self.read_timeseries(name)
                subcatchments[sub]["data"][data_type] = timeseries

        dis_eval = self.read_timeseries("dis_eval_kaemmerzell_79_89.txt")

        return dis_eval, subcatchments

    def read_timeseries(self, timeseries_name):
        """
        Loads in a timeseries and returns it

        :param timeseries_name:
        :param convert: Discharge needs to be converted.

        :return: timeseries
        """
        # Fixed model starting point
        # Change this if you want a warm up period other than a year
        begin = self.begin - relativedelta(years=1)
        step = datetime.timedelta(days=1)

        timeseries = cmf.timeseries(begin, step)
        timeseries.extend(float(value.strip("\n")) for value in open(
            timeseries_name))

#        if convert:
#            area_catchment = 562.41
#            timeseries *= 86400 * 1e3 / (area_catchment * 1e6)

        return timeseries

    @staticmethod
    def create_params():
        """
        Creates all the parameters needed.

        :return: List of parameters
        """
        param = spotpy.parameter.Uniform
        params = [param('tr_soil_gw', 1., 400.),
                  # tr_soil_out = residence time from soil to outlet
                  param("tr_soil_out", 1., 200.),
                  # tr_GW_out = Residence time in the groundwater to
                  #  the outlet
                  param('tr_gw_out', 1., 650.),
                  # V0_soil = field capacity for the soil
                  param('V0_soil', 1., 300.),
                  # param("V0_gw", 1, 300.),
                  # beta_soil_GW = Changes the flux curve of the soil
                  # to the groundwater
                  param('beta_soil_gw', 0.5, 6.0),
                  # beta_soil_out = exponent that changes the form of the
                  # flux from the soil to the outlet
                  param("beta_soil_out", 0.5, 7.0),
                  # ETV1 = the Volume that defines that the evaporation
                  # is lowered because of not enough water in the soil
                  param('ETV1', 1., 300.),
                  # fETV0 = factor the ET is multiplied by, when water is
                  #  low
                  param('fETV0', 0.1, 0.9),
                  # Rate of snow melt (for the low region)
                  param('meltrate', 0.01, 12.),
                  # Snow_melt_temp = Temperature at which the snow melts
                  # (needed because of averaged temp (for the low region)
                  param('snow_melt_temp', -3.0, 3.0),
                  # LAI = leaf area index
                  param('LAI', 1., 12),
                  # Canopy Closure
                  param("CanopyClosure", 0.1, 0.9)
                  ]
        return params

    def run_model(self):
        """
        Starts the model. Used by spotpy
        """
#        print("Start new model run at " + str(datetime.datetime.now()))
        try:
            # Create a solver for differential equations
            solver = cmf.CVodeIntegrator(self.project, 1e-8)

            # New time series for model results
            dis_sim = cmf.timeseries(self.begin, cmf.day)
            # starts the solver and calculates the daily time steps
            end = self.end

            for t in solver.run(self.project.meteo_stations[0].T.begin,
                                end, cmf.day):

                # Fill the results (first year is included but not used to
                # calculate the NS)
                # print(self.project.cells[0].layers[0].flux_to(self.outlet,t )
                # )
                if t >= self.begin:
                    dis_sim.add(self.outlet.waterbalance(t))

            return dis_sim
        # Return an nan - array when a runtime error occurs
        except RuntimeError:
            dis_sim = np.array(self.dis_eval[
                            self.begin:self.end + datetime.timedelta(days=1)])\
                      * np.nan
            # ET = np.array(self.dis_eval[
            #                 self.begin:self.end + datetime.timedelta
            # (days=1)])*np.nan
            return dis_sim

    def simulation(self, vector):
        """
        SpotPy expects a method simulation. This methods calls set_parameters
        and run_models, so SpotPy is satisfied
        """
        paramdict = dict((pp.name, v) for pp, v in zip(self.params, vector))
        self.set_parameters(paramdict)
        discharge = self.run_model()
        print("Simulation")
        print(discharge.begin, discharge.end)
        print(type(discharge))
        print(len(discharge))
        discharge = np.array(discharge)
        # CMF outputs discharge in m³/day
        # Measured discharge is in m³/s
        # Divide m³/day by 86400 to get to m³/s
        discharge /= 86400
        # self.discharge = discharge
        return discharge

    def evaluation(self):
        """
        For Spotpy
        """
        # plus one day because as in lists the last entry is not included in
        # datetime objects
        dis_eval = self.dis_eval[self.begin:self.end +
                                 datetime.timedelta(days=1)]
        # print("Evaluation")
        # print(dis_eval.begin, dis_eval.end)
        # print(type(dis_eval))
        # print(len(dis_eval))
        return np.array(dis_eval)

    def parameters(self):
        """
        For Spotpy
        """
        return spotpy.parameter.generate(self.params)

    @staticmethod
    def objectivefunction(simulation, evaluation):
        """
        For Spotpy
        """
        # Slice the arrays to only use the days of the calibration period
        # for objective function
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
        # if ns_calibration > 0:
        #     plt.plot(simulation * 86400, label="simulation_"+str(
        #         ns_calibration),
        #              alpha=0.7)
        #     plt.plot(evaluation * 86400, label="evaluation", alpha=0.7)
        #     plt.plot(model.ET, label="ET", alpha=0.7)
        #     plt.legend()
        #     name = "test_" + str(datetime.datetime.timestamp(
        #         datetime.datetime.now()))+".jpg"
        #     plt.savefig(name, dpi=300)
        #     plt.close()

        return [ns_calibration, ns_validation]


if __name__ == '__main__':
    # 1979 is spin up
    # 1980 till 1984 is used for calibration
    # 1985 till 1989 is used for validation
    begin = 1980
    end = 1989

    prefix = "semi_dis_height"

    runs = 10

    # File names of the forcing data
    subcatchment_names = ["grass", "wood", "rest", "crops"]

    # import algorithm
    from spotpy.algorithms import lhs as sampler

    # Find out if the model should run parallel (for supercomputer)
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'

    # Create the model
    model = SemiDisLanduse(datetime.datetime(begin, 1, 1),
                           datetime.datetime(end, 12, 31),
                           subcatchment_names)
    sampler = sampler(model, parallel=parallel, dbname="semi_dis_landuse",
                      dbformat="csv", save_sim=True, save_threshold=[0, 0])
    sampler.sample(runs)#, subsets=30)
    print(cmf.describe(model.project))
