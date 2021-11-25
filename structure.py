#! /usr/bin/env python3
import h5py as h5
import pandas as pd
import numpy as np

# Define all main parameters in the aircraft 
# MODEL GEOMETRY
# beam
span_main = 7.07
sweep_LE_main = np.deg2rad(20.)

chord_root = 0.471  #0.5048 (Graph)
chord_tip = 0.236 

# calculated inputs
x_tip = span_main*0.5*np.tan(sweep_LE_main)
sweep_quarter_chord = np.arctan((x_tip+chord_tip/4-chord_root/4)/(span_main*0.5))

# Fuselage information
length_fuselage = 3.44
offset_fuselage_vertical = 0
offset_wing_nose = 0.8842 + 0.09#0.8822765386648912
offset_tail_nose = 2.86236881559
sigma_fuselage = 10
m_bar_fuselage = 0.3
j_bar_fuselage = 0.1

# Tail
tail_chord_tip = 0.180325
tail_sweep_LE = np.deg2rad(19.51951)
tail_sweep_TE = np.deg2rad(18.0846)
half_tail_span = 1.318355
tail_span = 2*half_tail_span
# TODO: Get correct values

# calculated inputs
tail_x_tip = half_tail_span*np.tan(tail_sweep_LE)
tail_chord_root = tail_x_tip + tail_chord_tip - half_tail_span*np.tan(tail_sweep_TE)
tail_sweep_quarter_chord = np.arctan((tail_x_tip+tail_chord_tip/4-tail_chord_root/4)/(half_tail_span))

v_tail_angle = np.deg2rad(35.)
tail_sweep_quarter_chord
span_tail = 2.5
ea_tail = 0.5
sigma_tail = 10
m_bar_tail = 0.3
j_bar_tail = 0.1

# Material
material = "reference"
n_stiffness_wing = 17


# Ailerons
numb_ailerons = 4
y_coord_ailerons= np.array([0.862823, 2.820273, 4.301239, 5.653424, 6.928342])/2.

# Elevators
numb_elevators = 2
y_coord_elevators = np.array([0.258501, 0.788428, 1.318355])/2.

# Fuselage = 
y_coord_junction = 0.144



class FLEXOPStructure:

    def __init__(self, case_name, case_route, **kwargs):
        self.material = material
        self.sigma = kwargs.get('sigma', 1)
        self.n_elem_multiplier = kwargs.get('n_elem_multiplier', 1.5)

        self.route = case_route
        self.case_name = case_name

        self.thrust = kwargs.get('thrust', 0.)

        self.n_elem = None
        self.n_node = None
        self.n_node_elem = 3

        self.x = None
        self.y = None
        self.z = None

        self.n_elem_main = None
        self.n_elem_root_main = None
        self.n_elem_junction_main = None
        self.n_elem_per_aileron = None
        self.n_elem_tip_main = None

        self.n_elem_fuselage = None

        self.n_elem_tail = None
        self.n_elem_per_elevator = None
        self.n_elem_junction_tail = None

        self.n_node_main = None
        self.n_node_fuselage = None
        self.n_node_tail = None

        self.span_main = span_main
        self.span_tail = span_tail
        self.y_coord_junction = y_coord_junction # Radius fuselage at wing

        self.wing_only = True
        self.lifting_only = True

        self.n_stiffness_wing = n_stiffness_wing
        self.n_ailerons_per_wing = numb_ailerons
        self.n_elev_per_tail_surf = numb_elevators
        self.v_tail_angle = v_tail_angle
        self.tail_sweep_quarter_chord = tail_sweep_quarter_chord
        self.sweep_quarter_chord = sweep_quarter_chord

    def set_thrust(self, value):
        self.thrust = value

    def generate(self):
        # Set Elements
        self.n_elem_junction_main = int(0.5*self.n_elem_multiplier)
        if self.n_elem_junction_main < 1:
            self.n_elem_junction_main = 1
        self.n_elem_root_main = int(1*self.n_elem_multiplier)
        self.n_elem_tip_main = int(0.5*self.n_elem_multiplier)
        if self.n_elem_tip_main < 1:
            self.n_elem_tip_main = 1
        self.n_elem_per_aileron = int(4*self.n_elem_multiplier)
        self.n_elem_per_elevator =  int(3*self.n_elem_multiplier)
        self.n_elem_junction_tail = int(2*self.n_elem_multiplier)
        self.n_elem_main = int(self.n_elem_junction_main + self.n_elem_root_main + self.n_ailerons_per_wing * self.n_elem_per_aileron + self.n_elem_tip_main)
        self.n_elem_tail = int(self.n_elem_junction_tail + self.n_elev_per_tail_surf * self.n_elem_per_elevator)
        self.n_elem_fuselage = 21 #int(8*2*self.n_elem_multiplier)
 
        if self.wing_only:
            self.lifting_only = True

        # lumped masses
        df_lumped_masses = self.read_lumped_masses()
        n_lumped_mass = df_lumped_masses.shape[0]*2
        lumped_mass_nodes = np.zeros((n_lumped_mass, ), dtype=int)
        lumped_mass = np.zeros((n_lumped_mass, ))
        # lumped_mass[0] = 0 # mass_take_off
        lumped_mass_inertia = np.zeros((n_lumped_mass, 3, 3))
        lumped_mass_position = np.zeros((n_lumped_mass, 3))
        x_lumped_mass = 0.606 - offset_wing_nose
                # total number of elements
        self.n_elem = self.n_elem_main + self.n_elem_main
        if not self.wing_only:
            self.n_elem += self.n_elem_fuselage
            self.n_elem += self.n_elem_tail + self.n_elem_tail


        # number of nodes per part
        self.n_node_main = self.n_elem_main*(self.n_node_elem - 1) + 1
        self.n_node_tail = self.n_elem_tail*(self.n_node_elem - 1) + 1
        self.n_node_fuselage = (self.n_elem_fuselage+1)*(self.n_node_elem - 1) -1

        # total number of nodes
        self.n_node = self.n_node_main + self.n_node_main - 1
        if not self.wing_only:
            self.n_node += self.n_node_fuselage - 1
            self.n_node += self.n_node_tail - 1
            self.n_node += self.n_node_tail - 1


        # Aeroelastic properties
        n_stiffness = self.n_stiffness_wing
        if not self.wing_only:
            n_stiffness += 2

        n_mass = n_stiffness
        # TODO: Adjust for paper values
        m_bar_fuselage = 0.3*1.5
        j_bar_fuselage = 0.08

        sigma_tail = 100
        m_bar_tail = 0.3
        j_bar_tail = 0.08

        # beam
        self.x = np.zeros((self.n_node, ))
        self.y = np.zeros((self.n_node, ))
        self.z = np.zeros((self.n_node, ))
        structural_twist = np.zeros((self.n_elem, self.n_node_elem))
        beam_number = np.zeros((self.n_elem, ), dtype=int)
        frame_of_reference_delta = np.zeros((self.n_elem, self.n_node_elem, 3))
        conn = np.zeros((self.n_elem, self.n_node_elem), dtype=int)
        stiffness = np.zeros((n_stiffness, 6, 6))
        elem_stiffness = np.zeros((self.n_elem, ), dtype=int)
        mass = np.zeros((n_mass, 6, 6))
        elem_mass = np.zeros((self.n_elem, ), dtype=int)
        boundary_conditions = np.zeros((self.n_node, ), dtype=int)
        app_forces = np.zeros((self.n_node, 6))


            # Stiffness and mass properties
        list_stiffness_matrix, list_mass_matrix, y_cross_sections = self.load_stiffness_and_mass_matrix_from_matlab_file()
        for i in  range(self.n_stiffness_wing):
            stiffness[i, ...] = list_stiffness_matrix[i] #list_spanwise_stiffness_properties[i]
            mass[i, ...] =  list_mass_matrix[i] #list_spanwise_mass_properties[i]
        if not self.wing_only:
            stiffness[-2, ...] = list_stiffness_matrix[0]*sigma_fuselage
            stiffness[-1, ...] = list_stiffness_matrix[0]*sigma_tail

            mass[-2, ...] = self.generate_mass_matrix(m_bar_fuselage, j_bar_fuselage)
            mass[-1, ...] = self.generate_mass_matrix(m_bar_tail, j_bar_tail)


        ###############
        # right wing
        ###############
        we = 0
        wn = 0
        beam_number[we:we + self.n_elem_main] = 0
        # junction (part without ailerons)
        n_node_junctions = int(3 + 2*(self.n_elem_junction_main-1))
        self.y[wn:wn + n_node_junctions] = np.linspace(0.0, y_coord_junction, n_node_junctions)
        n_node_root = int(3 + 2*(self.n_elem_root_main-1))
        self.y[wn + n_node_junctions:wn + n_node_junctions+n_node_root-1] = np.linspace(y_coord_junction, y_coord_ailerons[0], n_node_root)[1:]
        
        
        # Approach 1: Direct transition from one aileron to another aileron
        n_nodes_per_cs = (self.n_elem_per_aileron)*2+1
        wn_end = 0
        n_node_tip = int(3 + 2*(self.n_elem_tip_main-1))
        for i_control_surface in range(self.n_ailerons_per_wing):
            wn_start = wn +  n_node_junctions -1 + n_node_root- 1 + i_control_surface*(n_nodes_per_cs-1)
            wn_end= wn_start + n_nodes_per_cs
            self.y[wn_start:wn_end] = np.linspace(y_coord_ailerons[i_control_surface], 
                                            y_coord_ailerons[i_control_surface+1], 
                                            n_nodes_per_cs)  
        # Aileron to tip
        self.y[wn_end:wn_end + n_node_tip-1] = np.linspace(y_coord_ailerons[-1], self.span_main*0.5, n_node_tip)[1:]
        self.x[wn+n_node_junctions:wn + self.n_node_main] += (abs(self.y[wn+n_node_junctions:wn + self.n_node_main])-y_coord_junction) * np.tan(self.sweep_quarter_chord)
    
        # Set stiffness, mass. For, and elastic axis (TODO: Change elastic axis to aero?)
        elem_stiffness[we:we + self.n_elem_main] = 0
        elem_mass[we:we + self.n_elem_main] = 0
        node_counter = wn
        # if self.ifting_only and not self.structural_properties_fuselage:
        #     aero_node[wn:wn + n_node_main] = True
        # else:
        #     aero_node[wn:wn + n_node_main] = abs(y[wn:wn + n_node_main]) >= y_coord_junction   
        for ielem in range(self.n_elem_main):
            conn[we + ielem, :] = ((np.ones((3, ))*(we + ielem)*(self.n_node_elem - 1)) +
                                [0, 2, 1])          
            # if all nodes of an element are aero_nodes, the element has aero material properties
            index_position= self.find_index_of_closest_entry(y_cross_sections,self.y[node_counter+ 2])
            
            if y_cross_sections[index_position] < self.y[node_counter+ 2]:
                i_material = index_position + 1
            else:
                i_material = index_position
            
            elem_stiffness[we + ielem] = i_material
            elem_mass[we + ielem] = i_material
            for inode in range(self.n_node_elem):
                frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]  
                # Elastic axis is later multiplied with chord in aero function
                # elastic_axis[ielem, inode] = (0.5 + list_spanwise_shear_center[i_material])
            node_counter += 2

        app_forces[wn] = [0, self.thrust, 0, 0, 0, 0]
        boundary_conditions[0] = 1
        # remember this is in B FoR

        ###############
        # left wing
        ###############
        we += self.n_elem_main
        wn += self.n_node_main

        # outer right wing

        beam_number[we:we + self.n_elem_main - 1] = 1
        # Mirror coordinates from left wing
        self.y[wn:wn + self.n_node_main - 1] = -self.y[1:self.n_node_main]
        self.x[wn:wn + self.n_node_main - 1] = self.x[1:self.n_node_main]
        self.z[wn:wn + self.n_node_main - 1] = self.z[1:self.n_node_main]


        elem_stiffness[we:we + self.n_elem_main] = 0
        elem_mass[we:we + self.n_elem_main] = 0


        for ielem in range(self.n_elem_main):
            conn[we + ielem, :] = ((np.ones((3, ))*(we+ielem)*(self.n_node_elem - 1)) +
                                [0, 2, 1])
            elem_stiffness[we + ielem] =elem_stiffness[ielem]
            elem_mass[we + ielem] = i_material
        
        for inode in range(self.n_node_elem):
            frame_of_reference_delta[we + ielem, inode, :] = [1.0, 0.0, 0.0] 

        conn[we, 0] = 0
        boundary_conditions[wn-1] = -1
        if self.wing_only: #TODO Check why?
            boundary_conditions[-1] = -1
        else:
            boundary_conditions[wn + self.n_node_main-1] = -1
        we += self.n_elem_main
        wn += self.n_node_main - 1


        # Set lumped masses wing
        for imass in range(int(n_lumped_mass/2)):
            lumped_mass[imass] =  df_lumped_masses.iloc[imass, 1]
            lumped_mass_position[imass, 0] = df_lumped_masses.iloc[imass, 2]
            lumped_mass_position[imass, 0] -= (0.5-0.140615385) *chord_root # adjust x=0 at LE
            lumped_mass_position[imass, 1] = df_lumped_masses.iloc[imass, 3]
            lumped_mass_position[imass, 2] = df_lumped_masses.iloc[imass, 4]
            lumped_mass_nodes[imass] = self.find_index_of_closest_entry(self.y[:self.n_node_main], lumped_mass_position[imass,1])
        
            # Copy to symmetric wing
            idx_symmetric = int(n_lumped_mass/2)+imass
            lumped_mass[idx_symmetric] =  lumped_mass[imass]
            lumped_mass_position[idx_symmetric, 0] = lumped_mass_position[imass, 0]
            lumped_mass_position[idx_symmetric, 1] = -lumped_mass_position[imass, 1]
            lumped_mass_position[idx_symmetric, 2] =  lumped_mass_position[imass, 2] 
            lumped_mass_nodes[idx_symmetric] = self.n_node_main -1 + lumped_mass_nodes[imass]

        if not self.wing_only:
            ###############
            # fuselage
            ###############
            beam_number[we:we + self.n_elem_fuselage] = 2
            x_fuselage = np.linspace(0.0, length_fuselage, self.n_node_fuselage) - offset_wing_nose
            z_fuselage = np.linspace(0.0, offset_fuselage_vertical, self.n_node_fuselage)
            idx_junction = self.find_index_of_closest_entry(x_fuselage, self.x[0])
            x_fuselage = np.delete(x_fuselage, idx_junction)
            z_fuselage = np.delete(z_fuselage, idx_junction)
            self.x[wn:wn + self.n_node_fuselage-1] = x_fuselage 
            self.z[wn:wn + self.n_node_fuselage-1] = z_fuselage
            adjust = False

            node_fuselage_conn = False
            for ielem in range(self.n_elem_fuselage):
                conn[we + ielem, :] = ((np.ones((3,))*(we + ielem)*(self.n_node_elem - 1)) +
                                    2 + [0, 2, 1]) - 1
                if adjust:
                    conn[we + ielem, :] -= 1
                else:
                    if node_fuselage_conn:
                        conn[we + ielem, 0] = 0
                    elif (conn[we + ielem, :] ==  wn+idx_junction).any():
                        adjust_elem = False
                        for idx_node in [0, 2, 1]:               
                            if adjust_elem:
                                conn[we + ielem, idx_node] -= 1  

                            elif conn[we + ielem, idx_node] ==  wn+idx_junction:
                                adjust = True
                                adjust_elem = True
                                conn[we + ielem, idx_node] = 0
                                if idx_node == 1:
                                    node_fuselage_conn = True
                for inode in range(self.n_node_elem):
                    frame_of_reference_delta[we + ielem, inode, :] = [0.0, 1.0, 0.0]
            # setup lumped mass position
            wn_lumped_mass = wn + self.find_index_of_closest_entry(self.x[wn:wn + self.n_node_fuselage-1], x_lumped_mass)
            lumped_mass_nodes[0] = wn_lumped_mass
            lumped_mass_position[0, 0] = self.x[wn_lumped_mass]
            lumped_mass_position[0, 1] = self.y[wn_lumped_mass]
            lumped_mass_position[0, 2] = self.z[wn_lumped_mass]
            boundary_conditions[wn] = - 1


            elem_stiffness[we:we + self.n_elem_fuselage] = n_stiffness - 2
            elem_mass[we:we + self.n_elem_fuselage] = n_stiffness - 2
            index_tail_start = wn + self.find_index_of_closest_entry(self.x[wn:wn + self.n_node_fuselage-1], offset_tail_nose-offset_wing_nose)
            we += self.n_elem_fuselage
            wn += self.n_node_fuselage - 1
            boundary_conditions[wn - 1] = -1
            boundary_conditions[index_tail_start] = 1

            ###############
            # right tail
            ###############
            beam_number[we:we + self.n_elem_tail] = 3
            self.x[wn:wn + self.n_node_tail - 1] = self.x[index_tail_start]
            wn_right_tail_start = wn
            n_node_junctions = int(3 + 2*(self.n_elem_junction_tail-1))
            self.y[wn:wn + n_node_junctions - 1] = np.linspace(0.0, y_coord_elevators[0], n_node_junctions)[:-1]
            # Approach 1: Direct transition from one aileron to another aileron
            n_nodes_per_cs = (self.n_elem_per_elevator)*2+1

            for i_control_surface in range(self.n_elev_per_tail_surf):
                wn_start = wn +  n_node_junctions - 1 + i_control_surface*(n_nodes_per_cs-1)
                wn_end= wn_start + n_nodes_per_cs
                self.y[wn_start:wn_end] = np.linspace(y_coord_elevators[i_control_surface], 
                                                    y_coord_elevators[i_control_surface+1], 
                                                    n_nodes_per_cs)
            self.x[wn:wn + self.n_node_tail - 1]  += abs(self.y[wn:wn + self.n_node_tail - 1])* np.tan(self.tail_sweep_quarter_chord)
            self.z[wn:wn + self.n_node_tail - 1] = self.z[index_tail_start]
            self.z[wn:wn + self.n_node_tail - 1] += self.y[wn:wn + self.n_node_tail - 1] * np.tan(self.v_tail_angle)

            elem_stiffness[we:we + self.n_elem_tail] = n_stiffness - 1
            elem_mass[we:we + self.n_elem_tail] = n_stiffness - 1
            for ielem in range(self.n_elem_tail):
                conn[we + ielem, :] = ((np.ones((3, ))*(we + ielem)*(self.n_node_elem - 1)) +
                                    [0, 2, 1])
                for inode in range(self.n_node_elem):
                    frame_of_reference_delta[we + ielem, inode, :] = [-1.0, 0.0, 0.0]     
            conn[we, 0] =  index_tail_start 
            boundary_conditions[wn + self.n_node_tail - 2] = -1
            we += self.n_elem_tail
            wn += self.n_node_tail - 1
            ###############
            # left tail
            ###############
            beam_number[we:we + self.n_elem_tail] = 4
            self.x[wn:wn + self.n_node_tail - 1] = self.x[index_tail_start]
            self.y[wn:wn + self.n_node_tail - 1] = -self.y[wn_right_tail_start:wn_right_tail_start + self.n_node_tail - 1]

            self.x[wn:wn + self.n_node_tail - 1]  += abs(self.y[wn:wn + self.n_node_tail - 1])* np.tan(self.tail_sweep_quarter_chord)
            self.z[wn:wn + self.n_node_tail - 1] = self.z[index_tail_start]
            self.z[wn:wn + self.n_node_tail - 1] += abs(self.y[wn:wn + self.n_node_tail - 1]) * np.tan(self.v_tail_angle)


            elem_stiffness[we:we + self.n_elem_tail] = n_stiffness - 1
            elem_mass[we:we + self.n_elem_tail] = n_stiffness - 1
            for ielem in range(self.n_elem_tail):
                conn[we + ielem, :] = ((np.ones((3, ))*(we + ielem)*(self.n_node_elem - 1)) +
                                    [0, 2, 1])
                for inode in range(self.n_node_elem):
                    frame_of_reference_delta[we + ielem, inode, :] = [1.0, 0.0, 0.0]

                node_counter += 2
            conn[we, 0] =  index_tail_start 
            boundary_conditions[wn + self.n_node_tail - 2] = -1
            we += self.n_elem_tail
            wn += self.n_node_tail - 1

        with h5.File(self.route + '/' + self.case_name + '.fem.h5', 'a') as h5file:
            h5file.create_dataset('coordinates', data=np.column_stack((self.x, self.y, self.z)))
            h5file.create_dataset('connectivities', data=conn)
            h5file.create_dataset('num_node_elem', data=self.n_node_elem)
            h5file.create_dataset('num_node', data=self.n_node)
            h5file.create_dataset('num_elem', data=self.n_elem)
            h5file.create_dataset('stiffness_db', data=stiffness)
            h5file.create_dataset('elem_stiffness', data=elem_stiffness)
            h5file.create_dataset('mass_db', data=mass)
            h5file.create_dataset('elem_mass', data=elem_mass)
            h5file.create_dataset('frame_of_reference_delta', data=frame_of_reference_delta)
            h5file.create_dataset('structural_twist', data=structural_twist)
            h5file.create_dataset('boundary_conditions', data=boundary_conditions)
            h5file.create_dataset('beam_number', data=beam_number)
            h5file.create_dataset('app_forces', data=app_forces)
            h5file.create_dataset('lumped_mass_nodes', data=lumped_mass_nodes)
            h5file.create_dataset('lumped_mass', data=lumped_mass)
            h5file.create_dataset('lumped_mass_inertia', data=lumped_mass_inertia)
            h5file.create_dataset('lumped_mass_position', data=lumped_mass_position)


    def load_stiffness_and_mass_matrix_from_matlab_file(self):
        import matlab.engine
        eng = matlab.engine.start_matlab()
        # Load data from file
        if self.material == "reference":
            file = '../01_case_files/FlexOp_Data_Jurij/dynamics_reference.mat'
        else:
            file = '../01_case_files/FlexOp_Data_Jurij/dynamics_tailored.mat'

        D = eng.load(file)
        matrices_cross_stiffness = np.array(D['dynamics'][0]['str']['elm']['C'])
        matrices_cross_mass = np.array(D['dynamics'][0]['str']['elm']['A'])
        matrices_cross_moment_of_inertia = np.array(D['dynamics'][0]['str']['elm']['I'])

        nodal_coordinates =np.array(D['dynamics'][0]['str']['xyz'])
        N_nodes = int(np.array(D['dynamics'][0]['str']['Nnode']))

        # Transform data
        coords = np.zeros((N_nodes, 3))
        counter = 0
        for irow in range(N_nodes):
            # skip first row
            coords[irow, :] = np.transpose(nodal_coordinates[counter:counter+3])
            counter += 3

        list_stiffness_matrix = []
        list_mass_matrix = []

        counter = 0
        inertia_counter = 0
        row_counter = 0
        while counter < matrices_cross_stiffness.shape[0]:
            list_stiffness_matrix.append(np.array(matrices_cross_stiffness[counter:counter+6, :]))
            mass_matrix = np.zeros((6,6))
            # mass distribution
            mass = float(matrices_cross_mass[row_counter])
            for i in range(3):
                mass_matrix[i,i] = mass
            mass_matrix[3:,3:] = matrices_cross_moment_of_inertia[inertia_counter:inertia_counter+3,:3]
            list_mass_matrix.append(mass_matrix)
            # TODO More elegant solution
            counter += 6
            inertia_counter += 3
            row_counter += 1
        return list_stiffness_matrix, list_mass_matrix, coords[1:,1]   

    def generate_mass_matrix(self, m_bar, j_bar):
        np.diag([m_bar, m_bar, m_bar, 
                j_bar, 0.5*j_bar, 0.5*j_bar])

    def find_index_of_closest_entry(self, array_values, target_value):
        return (np.abs(array_values - target_value)).argmin()

    def read_lumped_masses(self):
        file = '../01_case_files/flexOp_data/lumped_masses.csv'
        df = pd.read_csv(file, sep=';')
        print(df.head())
        return df