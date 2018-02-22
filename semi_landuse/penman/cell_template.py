# -*- coding: utf-8 -*-
"""
Created on Nov 28 14:09 2017
@author(s): Florian U. Jehn
"""
import cmf


class CellTemplate:
    """
    Creates all basic parts of a CMF Model
    """
    def __init__(self, project, subcatchment, name, outlet):
        self.name = name
        self.data = subcatchment["data"]
        self.height = subcatchment["height"]
        self.size = subcatchment["size"]
        self.project = project
        self.outlet = outlet
        self.cell = self.project.NewCell(0, 0, self.height, self.size * 1e6)
        self.basic_set_up()
        self.make_meteo_stations()

    def basic_set_up(self):
        """
        Creates all cmf structures which are the same for all cells.

        :param: cell: cmf.project.NewCell()
        :return: None (cell is changed though)
        """
        # Add Snow
        self.cell.add_storage("Snow", "S")
        cmf.Snowfall(self.cell.snow, self.cell)

        # Add layers and give them some water to start with
        self.cell.add_layer(2.0)
        self.cell.layers[0].volume = 10
        self.cell.add_layer(5.0)
        self.cell.layers[1].volume = 40

        # Create a storage for Interception
        I = self.cell.add_storage("Canopy", "C")

        # Install a connection for the ET
        cmf.PenmanMonteithET(self.cell.layers[0], self.cell.transpiration)

    def set_parameters(self, params):
        """
        Sets the parameters for the current cell.

        :param: params: dictionary of parameters.
        :return:
        """
        # Get all definition from the init method
        cell = self.cell
        outlet = self.outlet
        soil = cell.layers[0]
        gw = cell.layers[1]
        
        # EVT1 must be adjusted to cell size
        ETV1 = params["ETV1"]
        ETV1 = (ETV1 / 1000) * cell.area
        
        # V0 must be adjusted to cell size as well
        V0_soil = params["V0_soil"]
        V0_soil = (V0_soil / 1000) * cell.area
        
        # Adjustment of the ET
        cell.set_uptakestress(cmf.VolumeStress(
                                ETV1,
                                ETV1 * params["fETV0"]))

        # Flux from soil to outlet
        cmf.kinematic_wave(soil,
                           outlet,
                           params["tr_soil_out"] / V0_soil,
                           V0=V0_soil,
                           exponent=params["beta_soil_out"])

        # Flux from soil to groundwater
        cmf.kinematic_wave(soil, gw,
                           params["tr_soil_gw"] / V0_soil,
                           V0=V0_soil,
                           exponent=params["beta_soil_gw"])

        # Flux from the  groundwater to the outlet (baseflow)
        cmf.kinematic_wave(gw, outlet, params["tr_gw_out"])

        # Split the rainfall in interception and throughfall
        cmf.Rainfall(cell.canopy, cell, False, True)
        cmf.Rainfall(cell.surfacewater, cell, True, False)

        # Make an overflow for the interception storage
        cmf.RutterInterception(cell.canopy, cell.surfacewater, cell)

        # Transpiration from the plants is added
        cmf.CanopyStorageEvaporation(cell.canopy, cell.evaporation, cell)

        # Sets the paramaters for interception
        cell.vegetation.LAI = params["LAI"]

        # Defines how much throughfall there is (in %)
        cell.vegetation.CanopyClosure = params["CanopyClosure"]

        # # Set parameters of the snow calculations
        cmf.Weather.set_snow_threshold(params["snow_melt_temp"])
        cmf.SimpleTindexSnowMelt(cell.snow, soil, cell,
                                 rate=params["meltrate"])

    def make_meteo_stations(self):
        """
        Initializes the rain station and the meteo station for the current
        cell. The stations are used only for the current cell and have the
        same height as the cell.

        :return:
        """
        rainstation = self.project.rainfall_stations.add(
            "Rain Station " + self.name, self.data["prec"], (0, 0,
            self.height))

        rainstation.use_for_cell(self.cell)

        meteo_station = self.project.meteo_stations.add_station(
            "Meteo Station " + self.name, (0, 0, self.height))

        meteo_station.T = self.data["T_avg"]
        meteo_station.Tmin = self.data["T_min"]
        meteo_station.Tmax = self.data["T_max"]
        meteo_station.Windspeed = self.data["wind"]
        meteo_station.SetSunshineFraction(self.data["sunshine"])
        meteo_station.rHmean = self.data["rel_hum"]

        meteo_station.use_for_cell(self.cell)
