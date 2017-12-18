# -*- coding: utf-8 -*-
"""
Created on Nov 28 14:09 2017
@author(s): Florian U. Jehn
"""
from cell_template import CellTemplate
import cmf


class LowRegionCell(CellTemplate):
    """

    """
    def __init__(self, project: cmf.project):
        """
        Create a new cell. The creation of the basic structures is done by
        the parent class.

        :param project: A CMF project
        """
        self.low_cell = self.project.NewCell(0, 0, 0, 1000)
        super().__init__(self.low_cell)
        self.project = project

    def set_paremeters(self, params: dict):
        """
        Sets all parameters which are unique to the low region.

        :return: None
        """
        super().set_parameters(params)

    def make_rain_station(self, prec: list):
        pass


class HighRegionCell(CellTemplate):
    """

    """
    def __init__(self, project: cmf.project):
        self.high_cell = self.project.NewCell(30, 0, 0, 1000)
        super().__init__(self.high_cell)
        self.project = project

    def set_parameters(self, params: dict):
        """
        Sets all parameters which are unique to the high region

        :return: None
        """
        super().set_parameters(params)

