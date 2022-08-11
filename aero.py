#! /usr/bin/env python3
import h5py as h5
import numpy as np
import pandas as pd
from structure import span_main


# TODO:
# - solution for all the global parameters
# - set elastic axis
span_main = 7.07
half_wing_span = span_main*0.5
sweep_LE_main = np.deg2rad(20.)
chord_main_root = 0.471 
chord_main_tip = 0.236

# calculated inputs
x_tip = half_wing_span*np.tan(sweep_LE_main)
sweep_quarter_chord = np.arctan((x_tip+chord_main_tip/4-chord_main_root/4)/(half_wing_span))
sweep_TE_main= np.arctan((x_tip + chord_main_tip - chord_main_root)/(half_wing_span))


# Geometry parameters tail
chord_tail_root = 0.34 #0.21717159844088685
chord_tail_tip = 0.180325

v_tail_angle =  np.deg2rad(35.)
tail_sweep_LE = np.deg2rad(19.51951)
tail_sweep_TE = np.deg2rad(18.0846)
half_tail_span = 1.318355
tail_span = 2*half_tail_span

# Material
material = "reference"
n_stiffness_wing = 17


# Ailerons
numb_ailerons = 4
y_coord_ailerons= np.array([0.862823, 2.820273, 4.301239, 5.653424, 6.928342])/2.

# Elevators
numb_elevators = 2
y_coord_elevators = np.array([0.258501, 0.788428, 1.318355])/2.

# TODO: Adjust
chord_fin = 0.5
ea_main = 0.3
ea_fin = 0.5
ea_tail = 0.5

# reference area
area_ref = 2.54

y_coord_junction = 0.144


class FLEXOPAero:
    def __init__(self, m, structure, case_name, case_route, source_directory, **kwargs):
        """
        
        Key-Word Arguments:
            - cs_deflection (float): Elevator control surface deflection
            - rudder_deflection (float): rudder deflection
            - polars (list(np.array)): 4-column array for AoA (rad), Cl, Cd, Cm of each airfoil polar
        """
        self.m = m
        self.structure = structure

        self.route = case_route
        self.case_name = case_name

        self.cs_deflection = kwargs.get('cs_deflection', 0.)
        self.rudder_deflection = kwargs.get('rudder_deflection', 0.)

        self.chord_main_root = chord_main_root
        self.chord_tail_root = chord_tail_root        
        self.chord_main_tip = chord_main_tip
        self.chord_tail_tip = chord_tail_tip

        self.sweep_LE_main = sweep_LE_main
        self.sweep_TE_main = sweep_TE_main

        self.wing_only = self.structure.wing_only
        self.lifting_only = self.structure.lifting_only

        self.polars = kwargs.get('polars', None)
        self.source_directory = source_directory

    def generate(self):
        tail = not self.wing_only
        n_surfaces = 2
        if not self.wing_only:
            n_surfaces += 2

        structure = self.structure

        self.n_elem = structure.n_elem
        self.n_node_elem = structure.n_node_elem
        self.n_elem_main = structure.n_elem_main
        self.n_node_main = structure.n_node_main
        self.n_elem_fuselage = structure.n_elem_fuselage
        self.n_node_fuselage = structure.n_node_fuselage
        self.n_elem_tail = structure.n_elem_tail
        self.n_node_tail = structure.n_node_tail

        # aero
        self.airfoil_distribution = np.zeros((structure.n_elem, structure.n_node_elem), dtype=int)
        self.surface_distribution = np.zeros((structure.n_elem,), dtype=int) - 1
        self.surface_m = np.zeros((n_surfaces, ), dtype=int)
        m_distribution = 'uniform'
        self.aero_node = np.zeros((structure.n_node,), dtype=bool)
        self.twist = np.zeros((structure.n_elem, structure.n_node_elem))
        self.sweep = np.zeros((structure.n_elem, structure.n_node_elem))
        self.chord = np.zeros((structure.n_elem, structure.n_node_elem,))
        self.elastic_axis = np.zeros((structure.n_elem, structure.n_node_elem,))

        self.junction_boundary_condition_aero = np.zeros((1, n_surfaces), dtype=int) - 1

        ###############
        # Control Surfaces
        ###############
        n_control_surfaces = 4 
        if not self.wing_only:
            n_control_surfaces +=  2# on each side
        if not self.structure.symmetry_condition:
            n_control_surfaces *= 2
        self.control_surface = np.zeros((self.n_elem, self.n_node_elem), dtype=int) - 1
        control_surface_type = np.zeros((n_control_surfaces, ), dtype=int)
        control_surface_deflection = np.zeros((n_control_surfaces, ))
        control_surface_chord = np.zeros((n_control_surfaces, ), dtype=int)
        # aileron 1 right
        control_surface_type[0] = 0
        control_surface_deflection[0] =  0
        control_surface_chord[0] = self.m/4 # 0.25

        # aileron 2 right
        control_surface_type[1] = 0
        control_surface_deflection[1] = 0 
        control_surface_chord[1] = self.m/4 # 0.25

        # aileron 3 right
        control_surface_type[2] = 0
        control_surface_deflection[2] =  0
        control_surface_chord[2] = self.m/4 # 0.25
        
        # aileron 4 right
        control_surface_type[3] = 0
        control_surface_deflection[3] =  np.deg2rad(0)
        control_surface_chord[3] = self.m/4 # 0.25
        # TODO: Setup right elevator chord length
        if self.controllable:
            control_surface_type[:4] = 1
        if not self.wing_only:
            # rudder 1 - used for trim
            control_surface_type[4]  = 0
            control_surface_deflection[4]  = np.deg2rad(self.cs_deflection)
            control_surface_chord[4]  =  self.m/2 # Flexop's elevator cs have a ,chord of 36%. problems with aerogrid discretization
            # # rudder 2
            control_surface_type[5]  = 0
            control_surface_deflection[5]  = np.deg2rad(self.cs_deflection)# np.deg2rad(self.cs_deflection)
            control_surface_chord[5]  = self.m/2  # Flexop@s elevator cs have a ,chord of 36%. problems with aerogrid discretization
       
        if not self.structure.symmetry_condition:
            self.n_cs_right = int(n_control_surfaces/2)
            for i_cs_right in range(self.n_cs_right):
                i_cs_left = self.n_cs_right + i_cs_right
                control_surface_deflection[i_cs_left] = control_surface_deflection[i_cs_right] 
                control_surface_type[i_cs_left] = control_surface_type[i_cs_right] 
                control_surface_chord[i_cs_left] = control_surface_chord[i_cs_right] 
        ###############
        # right wing
        ###############
        we = 0
        wn = 0
        self.set_right_wing(we, wn)

        we += self.n_elem_main
        wn += self.n_node_main
        
        if not self.structure.symmetry_condition:
            self.set_left_wing(we, wn)            
            we += self.n_elem_main
            wn += self.n_node_main - 1
        if not self.wing_only and tail:
            we += self.n_elem_fuselage
            wn += self.n_node_fuselage - 1 - 1
            #
            
            self.set_right_tail(we, wn)
            we += self.n_elem_tail
            wn += self.n_node_tail
            
            if not self.structure.symmetry_condition:
                self.set_left_tail(we, wn)
                we += self.n_elem_tail
                wn += self.n_node_tail

        np.savetxt("control_surfaces.csv", self.control_surface)
        
        with h5.File(self.route + '/' + self.case_name + '.aero.h5', 'a') as h5file:
            airfoils_group = h5file.create_group('airfoils')
            # add one airfoil
            FLEXOP_airfoil = airfoils_group.create_dataset('0', data=np.column_stack(
                self.load_airfoil_data_from_file()))
            naca_airfoil_tail = airfoils_group.create_dataset('1', data=np.column_stack(
                self.generate_naca_camber(P=0, M=0)))
            naca_airfoil_fin = airfoils_group.create_dataset('2', data=np.column_stack(
                self.generate_naca_camber(P=0, M=0)))

            # chord
            chord_input = h5file.create_dataset('chord', data=self.chord)
            chord_input.attrs['units'] = 'm'

            # twist
            twist_input = h5file.create_dataset('twist', data=self.twist)
            twist_input.attrs['units'] = 'rad'

            # sweep
            sweep_input = h5file.create_dataset('sweep', data=self.sweep)
            sweep_input.attrs['units'] = 'rad'

            # airfoil distribution
            h5file.create_dataset('airfoil_distribution', data=self.airfoil_distribution)
            h5file.create_dataset('surface_distribution', data=self.surface_distribution)
            h5file.create_dataset('surface_m', data=self.surface_m)
            h5file.create_dataset('m_distribution', data=m_distribution.encode('ascii', 'ignore'))
            h5file.create_dataset('aero_node', data=self.aero_node)
            h5file.create_dataset('elastic_axis', data=self.elastic_axis)
            h5file.create_dataset('junction_boundary_condition', data=self.junction_boundary_condition_aero)
            h5file.create_dataset('control_surface', data=self.control_surface)
            h5file.create_dataset('control_surface_deflection', data=control_surface_deflection)
            h5file.create_dataset('control_surface_chord', data=control_surface_chord)
            h5file.create_dataset('control_surface_type', data=control_surface_type)

            if self.polars is not None:
                polars_group = h5file.create_group('polars')
                for i_airfoil in range(3):  # there are three airfoils
                    polars_group.create_dataset('{:g}'.format(i_airfoil), data=self.polars[i_airfoil])
    
    def set_right_wing(self, we, wn):
        i_surf = 0
        self.airfoil_distribution[we:we + self.n_elem_main, :] = 0
        self.surface_distribution[we:we + self.n_elem_main] = i_surf
        self.surface_m[i_surf] = self.m

        if self.lifting_only:
            self.aero_node[wn:wn + self.n_node_main] = True
        else:
            self.aero_node[wn:wn + self.n_node_main] = abs(self.structure.y[wn:wn + self.n_node_main]) >= y_coord_junction  
 
        n_node_junctions = int(3 + 2*(self.structure.n_elem_junction_main-1))
        self.junction_boundary_condition_aero[0, i_surf] = 1 # BC at fuselage junction
        temp_chord = np.zeros((self.n_node_main)) + self.chord_main_root
        temp_chord[n_node_junctions:self.n_node_main] = abs(self.structure.y[wn+n_node_junctions:wn +self.n_node_main]*np.tan(self.sweep_LE_main)-(self.chord_main_root + self.structure.y[wn+n_node_junctions:wn + self.n_node_main]*np.tan(self.sweep_TE_main)))
        temp_sweep = np.linspace(0.0, 0*np.pi/180, self.n_node_main)

        node_counter = 0
        global_node_counter = wn

        list_spanwise_shear_center = self.read_spanwise_shear_center()
        for i_elem in range(we, we + self.n_elem_main):
            for i_local_node in range(self.n_node_elem):
                if not i_local_node == 0:
                    node_counter += 1
                inode = node_counter
                if i_local_node == 1:                
                    inode += 1
                elif i_local_node == 2:
                    inode -= 1
                self.chord[i_elem, i_local_node] = temp_chord[inode]
                self.sweep[i_elem, i_local_node] = temp_sweep[node_counter]
                # get jig twist            
                self.twist[i_elem, i_local_node] = -self.get_jigtwist_from_y_coord(self.structure.y[wn + inode])
                self.elastic_axis[i_elem, i_local_node] = list_spanwise_shear_center[self.structure.elem_stiffness[i_elem]]
            global_node_counter += 2

        node_counter = 0
        cs_counter = -1
        cs_surface = False
        for i_elem in range(we, we + self.n_elem_main):
            for i_local_node in [0,2,1]:
                if not i_local_node == 0:
                    node_counter += 1
                if abs(self.structure.y[node_counter]) == y_coord_ailerons[0] and i_local_node == 0:
                    cs_surface = True 
                if cs_surface:
                    if abs(self.structure.y[node_counter]) in y_coord_ailerons:
                        if i_local_node == 0:
                            cs_counter += 1
                    self.control_surface[i_elem, i_local_node] = cs_counter
                if abs(self.structure.y[node_counter]) >= y_coord_ailerons[-1]:
                    cs_surface = False
    def set_left_wing(self, we, wn):
        i_surf = 1
        self.airfoil_distribution[we:we + self.n_elem_main] = 0
        self.surface_distribution[we:we + self.n_elem_main] = i_surf
        self.surface_m[i_surf] = self.m

        if self.lifting_only:
            self.aero_node[wn:wn + self.n_node_main] = True
        else:
            self.aero_node[wn:wn + self.n_node_main] = self.structure.y[wn:wn + self.n_node_main] <= -y_coord_junction

        self.junction_boundary_condition_aero[0, i_surf] = 0 # BC at fuselage junction
        node_counter = 0
        for i_elem in range(we, we + self.n_elem_main):
            for i_local_node in range(self.n_node_elem): 
                self.twist[i_elem, i_local_node] = self.twist[i_elem - we, i_local_node] 
                self.chord[i_elem, i_local_node] = self.chord[i_elem-we, i_local_node]
                self.elastic_axis[i_elem, i_local_node] = self.elastic_axis[i_elem - we, i_local_node]
                self.sweep[i_elem, i_local_node] = self.sweep[i_elem-we, i_local_node] 
                # if self.control_surface[i_elem-we, i_local_node] >= 0:
                #     self.control_surface[i_elem, i_local_node] = self.control_surface[i_elem-we, i_local_node]  + self.n_cs_right + 1


        # For control surfaces setup
        node_counter = 0
        cs_counter = -1
        cs_surface = False
        for i_elem in range(we, we + self.n_elem_main):
            for i_local_node in [0,2,1]:
                if not i_local_node == 0:
                    node_counter += 1
                if abs(self.structure.y[node_counter]) == y_coord_ailerons[0] and i_local_node == 0:
                    cs_surface = True 
                if cs_surface:
                    if abs(self.structure.y[node_counter]) in y_coord_ailerons:
                        if i_local_node == 0:
                            cs_counter += 1
                    self.control_surface[i_elem, i_local_node] = cs_counter + self.n_cs_right
                if abs(self.structure.y[node_counter]) >= y_coord_ailerons[-1]:
                    cs_surface = False

    def set_right_tail(self, we, wn):
        i_surf = 2
        if self.structure.symmetry_condition:
            i_surf -= 1
        
        self.airfoil_distribution[we:we + self.n_elem_tail, :] = 2
        # self.airfoil_distribution[wn:wn + self.n_node_tail] = 0
        self.surface_distribution[we:we + self.n_elem_tail] = i_surf
        self.surface_m[i_surf] = self.m
        # XXX not very elegant
        # aero_node[wn:] = True

        if self.lifting_only:
            self.aero_node[self.structure.index_tail_start] = True
            self.aero_node[wn:wn + self.n_node_tail] = True
        else:
            self.aero_node[wn:wn + self.n_node_tail] = self.structure.y[wn:wn + self.n_node_tail] >= 0.04
        self.junction_boundary_condition_aero[0, i_surf] = 3 # BC at fuselage junction
        
        temp_chord = self.chord_tail_root - abs(self.structure.y[wn:wn + self.n_node_tail]*np.tan(tail_sweep_LE)) + abs(self.structure.y[wn:wn + self.n_node_tail]*np.tan(tail_sweep_TE))

        node_counter = 0
        for i_elem in range(we, we + self.n_elem_tail):
            for i_local_node in range(self.n_node_elem):
                self.twist[i_elem, i_local_node] = -0
        for i_elem in range(we, we + self.n_elem_tail):
            for i_local_node in range(self.n_node_elem):
                if not i_local_node == 0:
                    node_counter += 1
                    if i_local_node == 1:
                        self.chord[i_elem, i_local_node] = temp_chord[node_counter + 1]
                    elif i_local_node == 2:
                        self.chord[i_elem, i_local_node] = temp_chord[node_counter - 1]
                else:
                    self.chord[i_elem, i_local_node] = temp_chord[node_counter]            
                self.elastic_axis[i_elem, i_local_node] = ea_main
        node_counter = wn - 2
        cs_counter = -1
        
        cs_surface = False
        for i_elem in range(we, we + self.n_elem_tail):
            for i_local_node in range(3):
                if not i_local_node == 0:
                    node_counter += 1
                if abs(self.structure.y[node_counter]) == y_coord_elevators[0] and i_local_node == 0:
                    cs_surface = True 
                if cs_surface:
                    if abs(self.structure.y[node_counter]) in y_coord_elevators:
                        if i_local_node == 0:
                            if cs_counter == -1:
                                cs_counter = 4
                            else:
                                cs_counter += 1
                    # TODO: introduce flag!
                    if cs_counter == 5:
                        cs_counter = 4
                    self.control_surface[i_elem, i_local_node] = cs_counter
                if abs(self.structure.y[node_counter]) >= y_coord_elevators[-1]:
                    cs_surface = False
            
    def set_left_tail(self, we, wn):
        i_surf = 3
        self.airfoil_distribution[we:we + self.n_elem_tail, :] = 2
        self.surface_distribution[we:we + self.n_elem_tail] = i_surf
        self.surface_m[i_surf] = self.m

        if self.lifting_only:
            self.aero_node[wn:wn + self.n_node_tail] = True
        else:
            self.aero_node[wn:wn + self.n_node_tail] = self.structure.y[wn:wn + self.n_node_tail] <= -0.04
        self.junction_boundary_condition_aero[0, i_surf] = 2 # BC at fuselage junction
        node_counter = 0
        for i_elem in range(we, we + self.n_elem_tail):
            for i_local_node in range(self.n_node_elem):
                self.twist[i_elem, i_local_node] = -0     
                self.elastic_axis[i_elem, i_local_node] = ea_main
                self.chord[i_elem, i_local_node] = self.chord[i_elem-self.n_elem_tail, i_local_node]

        for i_elem in range(we, we + self.n_elem_tail):
            if self.control_surface[i_elem - self.n_elem_tail, 0] > -1:
                self.control_surface[i_elem, :] = self.control_surface[i_elem - self.n_elem_tail, :] #+ self.n_cs_right

    def get_jigtwist_from_y_coord(self, y_coord):
        y_coord = abs(y_coord)
        df_jig_twist = pd.read_csv(self.source_directory + '/jig_twist.csv',
                                sep=';')
        idx_closest_value = self.find_index_of_closest_entry(df_jig_twist.iloc[:,0], y_coord)
        if self.structure.material == "reference":
            column = 1
        else:
            column = 2
        if idx_closest_value == df_jig_twist.shape[0]:
            idx_adjacent = idx_closest_value - 1 
        elif idx_closest_value == 0 or df_jig_twist.iloc[idx_closest_value,0] < y_coord:
            idx_adjacent = idx_closest_value + 1
        else:
            idx_adjacent = idx_closest_value - 1  
        
        
        jig_twist_interp = df_jig_twist.iloc[idx_closest_value,column] + ((y_coord - df_jig_twist.iloc[idx_closest_value, 0]) 
                                                    / (df_jig_twist.iloc[idx_adjacent, 0] - df_jig_twist.iloc[idx_closest_value,0])
                                                    *(df_jig_twist.iloc[idx_adjacent, column] - df_jig_twist.iloc[idx_closest_value,column]))
        # when the denominator of the interpolation is zero
        if np.isnan(jig_twist_interp):
            jig_twist_interp = df_jig_twist.iloc[idx_closest_value, 1]
        return np.deg2rad(jig_twist_interp)


    def generate_naca_camber(self,M=0, P=0):
        mm = M*1e-2
        p = P*1e-1

        def naca(x, mm, p):
            if x < 1e-6:
                return 0.0
            elif x < p:
                return mm/(p*p)*(2*p*x - x*x)
            elif x > p and x < 1+1e-6:
                return mm/((1-p)*(1-p))*(1 - 2*p + 2*p*x - x*x)

        x_vec = np.linspace(0, 1, 1000)
        y_vec = np.array([naca(x, mm, p) for x in x_vec])
        return x_vec, y_vec

    def load_airfoil_data_from_file(self):
        file = self.source_directory + "/camber_line_airfoils.csv"
        camber_line = pd.read_csv(file, sep = ";")
        return np.array(camber_line.iloc[:,0]), np.array(camber_line.iloc[:,1])

    def find_index_of_closest_entry(self, array_values, target_value):
        return np.argmin(np.abs(array_values - target_value))

    def read_spanwise_shear_center(self):
        reference_shear_center = 0.71 # given by Jurij
        df = pd.read_csv(self.source_directory + '/shear_center.csv',
                                sep=';')
        if self.structure.material == "reference":
            column = 1
        else:
            column = 2
        return (reference_shear_center + df.iloc[:,column]).to_list()