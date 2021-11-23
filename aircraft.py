#! /usr/bin/env python3
import h5py as h5
import configobj
import numpy as np
from structure import FLEXOPStructure
from aero import FLEXOPAero
import os
import sharpy.sharpy_main


class FLEXOP:

    def __init__(self, case_name, case_route, output_route):
        self.case_name = case_name
        self.case_route = case_route
        self.output_route = output_route

        self.structure = None
        self.aero = None
        self.fuselage = None

        self.settings = None

    def init_structure(self, **kwargs):
        self.structure = FLEXOPStructure(self.case_name, self.case_route, **kwargs)

    def init_aero(self, m, **kwargs):
        self.aero = FLEXOPAero(m, self.structure, self.case_name, self.case_route, **kwargs)
    def set_flight_controls(self, thrust=0., elevator=0., rudder=0.):
        self.structure.set_thrust(thrust)

        if self.aero is not None:
            self.aero.cs_deflection = elevator
            self.aero.rudder_deflection = rudder

    def generate(self):

        if not os.path.isdir(self.case_route):
            os.makedirs(self.case_route)

        self.structure.generate()

        if self.aero is not None:
            self.aero.generate()

    def create_settings(self, settings):
        file_name = self.case_route + '/' + self.case_name + '.sharpy'
        config = configobj.ConfigObj()
        config.filename = file_name
        for k, v in settings.items():
            config[k] = v
        config.write()
        self.settings = settings

    def clean(self):
        """
        clean_test_files

        Removes the previous h5 files

        Args:
            route (string): path of the case
            case_name (string): name of the case

        """
        list_files = ['.fem.h5', '.aero.h5', '.nonlifting_body.h5', '.dyn.h5', '.mb.h5', '.sharpy', '.flightcon.txt']
        for file in list_files:
            path_file= self.route + '/' + self.case_name + file
            if os.path.isfile(path_file):
                os.remove(path_file)

    def run(self):
        sharpy.sharpy_main.main(['', self.case_route + '/' + self.case_name + '.sharpy'])
