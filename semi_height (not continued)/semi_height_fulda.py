# -*- coding: utf-8 -*-
"""
Created on Nov 28 14:08 2017
@author(s): Florian U. Jehn

A semi distributed model for the Fulda catchment till Kaemmerzell.
The seperating factor is the height. The structure of lumped_intermediate is
used for both the low and the high region. The snow parameters are
calibrated seperately for both subcatchments.
"""
import cells
import cmf
import datetime
import os
import numpy as np
import spotpy


class SemiDisHeightFulda:
    def __init__(self, begin: datetime.datetime, end: datetime.datetime):
        """

        :param begin:
        :param end:
        """
        project = cmf.project()
        # Add outlet
        self.outlet = project.NewOutlet("Outlet", 50, 0, 0)
        self.project = project
        self.low = cells.LowRegionCell(project)
        self.high = cells.HighRegionCell(project)
        self.begin = begin
        self.end = end
        self.params = self.create_params()

    def set_parameters(self, paramdict: dict):
        """

        :return:
        """
        self.low.set_parameters(paramdict)
        self.high.set_parameters(paramdict)

    @staticmethod
    def create_params():
        """

        :return:
        """
        param = spotpy.parameter.Uniform
        params = [param('tr_soil_gw', 0., 400.),
                  # tr_soil_out = residence time from soil to outlet
                  param("tr_soil_out", 0., 200.),
                  # tr_GW_out = Residence time in the groundwater to
                  #  the outlet
                  param('tr_gw_out', 0., 650.),
                  # V0_soil = field capacity for the soil
                  param('V0_soil', 0., 300.),
                  # beta_soil_GW = Changes the flux curve of the soil
                  # to the groundwater
                  param('beta_soil_gw', 0., 6.0),
                  # beta_soil_out = exponent that changes the form of the
                  # flux from the soil to the outlet
                  param("beta_soil_out", 0., 8.0),
                  # ETV1 = the Volume that defines that the evaporation
                  # is lowered because of not enough water in the soil
                  param('ETV1', 0., 300.),
                  # fETV0 = factor the ET is multiplied by, when water is
                  #  low
                  param('fETV0', 0., 0.9),
                  # Rate of snow melt (for the low region)
                  param('meltrate_low', 0.01, 15.),
                  # Snow_melt_temp = Temperature at which the snow melts
                  # (needed because of averaged temp (for the low region)
                  param('snow_melt_temp_low', -3.0, 3.0),
                  # Rate of snow melt (for the high region)
                  param('meltrate_high', 0.01, 15.),
                  # Snow_melt_temp = Temperature at which the snow melts
                  # (needed because of averaged temp (for the high region)
                  param('snow_melt_temp_high', -3.0, 3.0)
                  ]
        return params

    def run_model(self):
        """
        Starts the model. Used by spotpy
        """

        try:
            # Create a solver for differential equations
            solver = cmf.CVodeIntegrator(self.project, 1e-8)

            # New time series for model results
            dis_sim = cmf.timeseries(self.begin, cmf.day)
            # starts the solver and calculates the daily time steps
            end = self.end
            for t in solver.run(self.project.meteo_stations[0].T.begin, end,
                                cmf.day):
                # Fill the results (first year is included but not used to
                # calculate the NS)
                if t >= self.begin:
                    dis_sim.add(self.outlet.waterbalance(t))
            return dis_sim
        # Return an nan - array when a runtime error occurs
        except RuntimeError:
            return np.array(self.dis_eval[
                            self.begin:self.end])*np.nan

    def simulation(self, vector):
        """
        SpotPy expects a method simulation. This methods calls set_parameters
        and run_models, so SpotPy is satisfied
        """
        paramdict = dict((pp.name, v) for pp, v in zip(self.params, vector))
        self.set_parameters(**paramdict)
        discharge = self.run_model()
        return np.array(discharge)

    def evaluation(self):
        """
        For Spotpy
        """
        return np.array(
            self.Q[self.begin:self.end])

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
        return [ns_calibration, ns_validation]


if __name__ == '__main__':


    # 1979 is spin up
    # 1980 till 1984 is used for calibration
    # 1985 till 1989 is used for validation
    begin = 1980
    end = 1989

    prefix = "semi_dis_height"

    runs = 100

    # File names of the forcing data
    discharge = "Q_Kammerzell_1979_1999.txt"
    temperature = "T_kammerzell_1979_1999_max_min_avg.txt"

    prec_low = "prec_kaemmerzell_low_regions.txt_79_89.txt"
    prec_high = "prec_kaemmerzell_high_regions.txt_79_89.txt"

    # import algorithm
    from spotpy.algorithms import rope as Sampler

    # Find out if the model should run parallel (for supercomputer)
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'

    # Create the model
    model = SemiDisHeightFulda(datetime.datetime(begin, 1, 1),
                               datetime.datetime(end, 12, 31))

    sampler = Sampler(model, parallel=parallel, dbname="intermediate_lumped",
                      dbformat="csv", save_sim=True)
    sampler.sample(runs, subsets=30)
