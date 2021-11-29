#! /usr/bin/env python3
import h5py as h5
import configobj
import numpy as np
from structure import FLEXOPStructure
from fuselage import FLEXOPFuselage
from aero import FLEXOPAero
from aero import area_ref
import os
import sharpy.sharpy_main


class FLEXOP:

    def __init__(self, case_name, case_route, output_route, datafiles_directory='./aeroelastic_properties'):
        self.case_name = case_name
        self.case_route = case_route
        self.output_route = output_route

        self.structure = None
        self.aero = None
        self.fuselage = None

        self.settings = None

        self.source_directory = os.path.abspath(datafiles_directory)
        print(f'Looking for source files in {self.source_directory}')

    def init_structure(self, **kwargs):
        self.structure = FLEXOPStructure(self.case_name, self.case_route, **kwargs)
        self.structure.source_directory = self.source_directory

    def init_aero(self, m, **kwargs):
        self.aero = FLEXOPAero(m, self.structure, self.case_name, self.case_route, **kwargs)
        self.aero.source_directory = self.source_directory

    def init_fuselage(self, m, **kwargs):
        self.fuselage = FLEXOPFuselage(m, self.structure, self.case_name, self.case_route, **kwargs)

    def set_flight_controls(self, thrust=0., elevator=0., rudder=0.):
        self.structure.set_thrust(thrust)

        if self.aero is not None:
            self.aero.cs_deflection = elevator
            self.aero.rudder_deflection = rudder

    @property
    def reference_area(self):
        return area_ref

    def generate(self):

        if not os.path.isdir(self.case_route):
            os.makedirs(self.case_route)

        self.structure.generate()

        if self.aero is not None:
            self.aero.generate()
        if self.fuselage is not None:
            self.fuselage.generate()

    def create_settings(self, settings):
        file_name = self.case_route + '/' + self.case_name + '.sharpy'
        config = configobj.ConfigObj()
        config.filename = file_name
        for k, v in settings.items():
            config[k] = v
        config.write()
        self.settings = settings

    def clean(self):
        list_files = ['.fem.h5', '.aero.h5', '.nonlifting_body.h5', '.dyn.h5', '.mb.h5', '.sharpy', '.flightcon.txt']
        for file in list_files:
            path_file = self.case_route + '/' + self.case_name + file
            if os.path.isfile(path_file):
                os.remove(path_file)

    def run(self):
        sharpy.sharpy_main.main(['', self.case_route + '/' + self.case_name + '.sharpy'])
