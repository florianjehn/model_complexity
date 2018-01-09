# -*- coding: utf-8 -*-
"""
Created on Jan 09 10:39 2018
@author(s): Florian U. Jehn
"""

import sys
import datetime
import cmf
import spotpy
from spotpy.parameter import Uniform
import os
import numpy as np
import pandas as pd


class ScalingTester:
    """
    A class to determine how CMF handles different amounts of cells. To test
    this the same model structure is run as a 1, 2, 4 and 8 cell layout. Each
    cell has the same structure and parameters.
    """
    # Catchment area km²
    area = 562.41

    # General storage parameter
    V0 = Uniform(10, 10000, 1000)

    # ET parameter
    fETV1 = Uniform(0.01, 1, 0.2, doc='if V<fETV1*V0, water uptake stress for plants starts')
    fETV0 = Uniform(0, 0.9, 0.2, doc='if V<fETV0*fETV1*V0, plants die of drought')

    # Outflow parameters
    tr = Uniform(0.1, 1000, doc='Residence time of water in storage when V=V0')

    def __init__(self, begin=None, end=None, num_cells=None):
        """
        Initializes the model.

        :param begin: Start year for the calibration
        :param end: Stop year
        :param num_cells: Number of cells used for this layout
        :return: None
        """
        self.dbname = "scaling_tester_" + str(num_cells)

        # load driver data
        self.data = DataProvider("fulda_kaemmerzell_climate_79_89.csv")
        # Create cells and project
        self.project, self.outlet = self.create_project()
        self.num_cells = num_cells
        self.cells = self.create_cells()

        # Add the data and set the parameters with random value, so the
        # complete structure can be described.
        self.data.add_stations(self.project)
        self.begin = begin or self.data.begin
        self.end = end or self.data.end
        self.setparameters()

    def create_project(self):
        """
        Creates and CMF project with an outlet and other basic stuff and
        returns it.
        :return: cmf project and cmf outlet
        """
        # Use only a single thread, that is better for a calibration run and
        # for small models
        cmf.set_parallel_threads(1)

        # make the project
        p = cmf.project()

        # make the outlet
        outlet = p.NewOutlet("outlet", 10, 0, 0)
        return p, outlet

    def create_cells(self):
        """
        Creates a 'num_cells' amount of cells for the project
        :return:
        """
        # Adjust the cellsize to the amount of cells
        area = self.area / self.num_cells
        # Create all the cells!
        cells = []
        for num in range(self.num_cells):
            cells.append(CellTemplate(self.project, self.outlet, area,
                                      num))
        return cells

    def setparameters(self, par=None):
        """
        Sets the parameters for all cells seperately
        :return:
        """
        # Create tje parameters
        par = par or spotpy.parameter.create_set(self)
        # Call all cells
        for cell in self.cells:
            cell.set_parameters(par)

    def runmodel(self):
        """
        Runs the models and saves the results.

        :return: Simulated discharge
        """
        solver = cmf.CVodeIntegrator(self.project, 1e-9)

        # Result timeseries
        res_q = cmf.timeseries(self.begin, cmf.day)

        # Start solver and calculate in daily steps
        for t in solver.run(self.data.begin, self.end, cmf.day):
            res_q.add(self.outlet.waterbalance(t))

        return res_q

    def simulation(self, vector=None):
        """
        Sets the parameters of the model and starts a run
        :return: np.array with runoff in mm/day
        """
        self.setparameters(vector)
        result_q = self.runmodel()
        return np.array(result_q[self.begin:self.end])

    def evaluation(self):
        """Returns the evaluation data"""
        runoff_mm = self.data.runoff_mm(self.area)

        return np.array(
                runoff_mm[self.begin:self.end])

    def objectivefunction(self, simulation, evaluation):
        return spotpy.objectivefunctions.nashsutcliffe(evaluation, simulation)


class CellTemplate:
    """
    Template, which provides
    """
    def __init__(self, project, outlet, area, cell_num):
        self.project = project
        self.outlet = outlet
        self.area = area
        self.cell = self.project.NewCell(cell_num, 0, 0, area * 1e6)
        self.basic_set_up()

    def basic_set_up(self):
        """
        Creates the basic storages, that are to be connected in set_parameters.
        :return:
        """
        # Add layers
        self.cell.add_layer(2.0)
        # Install a connection for the ET
        cmf.HargreaveET(self.cell.layers[0], self.cell.transpiration)

    def set_parameters(self, par):
        """
        Sets the parameters for a cell instance
        :param par: Object with all parameters
        :return: None
        """
        c = self.cell
        out = self.outlet

        # Set uptake stress
        ETV1 = par.fETV1 * par.V0
        ETV0 = par.fETV0 * ETV1
        c.set_uptakestress(cmf.VolumeStress(ETV1, ETV0))

        # Connect layer with outlet
        cmf.LinearStorageConnection(c.layers[0], out, par.tr)


class DataProvider:
    """
    Holds the forcing and calibration data
    """
    def __init__(self, file_name):
        # Load data from file using numpy magic
        data = pd.read_csv(file_name, encoding="ISO-8859-1", sep=";")
        # Delete first row, as it only contains the units
        data = data.iloc[1:]
        data = data.dropna(axis=0)

        def bstr2date(bs):
            """
            Helper function to convert date byte string to datetime object
            """
            return datetime.datetime.strptime(bs, '%d.%m.%Y')

        # Get begin, step and end from the date column
        self.begin = bstr2date(data["date"].iloc[0])
        self.step = bstr2date(data["date"].iloc[1]) - self.begin
        self.end = bstr2date(data["date"].iloc[-1])

        # Read in the data
        self.P = cmf.timeseries.from_sequence(self.begin, self.step,
                                              data["Prec"])
        self.T = cmf.timeseries.from_sequence(self.begin, self.step,
                                              data["tmean"])
        self.Tmin = cmf.timeseries.from_sequence(self.begin, self.step, data["tmin"])
        self.Tmax = cmf.timeseries.from_sequence(self.begin, self.step, data["tmax"])
        self.Q = cmf.timeseries.from_sequence(self.begin, self.step, data["Q"])

    def runoff_mm(self, area):
        sec_per_day = 86400
        mm_per_m = 1000
        return self.Q * sec_per_day / area * mm_per_m

    def add_stations(self, project):
        """
        Creates a rainstation and a meteo station for the cmf project
        :param project: A cmf.project
        :return: rainstation, meteo
        """
        rainstation = project.rainfall_stations.add('Grebenau avg', self.P,
                                                    (0, 0, 0))

        project.use_nearest_rainfall()

        # Temperature data
        meteo = project.meteo_stations.add_station('Grebenau avg',
                                                   (0, 0, 0))
        meteo.T = self.T
        meteo.Tmin = self.Tmin
        meteo.Tmax = self.Tmax

        project.use_nearest_meteo()
        return rainstation, meteo


if __name__ == '__main__':
    # Get sampler
    from spotpy.algorithms import lhs as Sampler

    # Check if we are running on a supercomputer or local
    parallel = 'mpi' if 'OMPI_COMM_WORLD_SIZE' in os.environ else 'seq'

    # Create the model
    model = ScalingTester(num_cells=1)

    # Get number of runs
    if 'SPOTPYRUNS' in os.environ:
        # from environment
        runs = int(os.environ['SPOTPYRUNS'])
    elif len(sys.argv) > 1:
        # from command line
        runs = int(sys.argv[1])
    else:
        # run once
        runs = 50

    # Create the sampler
    sampler = Sampler(model, parallel=parallel, dbname=model.dbname,
                      dbformat='csv', save_sim=True)

    # Print our configuration
    # print(spotpy.describe.describe(sampler))
    # Print the cmf setup
    print(cmf.describe(model.project))

    # Do the sampling
    if runs > 1:
        # Now we can sample with the implemented Monte Carlo algortihm:
        sampler.sample(runs)
    else:
        result = model.simulation()
        for name, value in spotpy.objectivefunctions.calculate_all_functions(
                model.evaluation(), result):
            try:
                print('{:>30.30s} = {:0.6g}'.format(name, value))
            except ValueError:
                print('{:>30.30s} = {}'.format(name, value))
