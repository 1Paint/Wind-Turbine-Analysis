# Test for wind history velocities of 0.
# Test for 1-line history.

"""This module holds a class that loads constants from files.
All values are in SI units unless stated otherwise.
"""

import numpy as np

class Constants():

    def open_design(self, file):
        design = np.loadtxt(file)
        
        # Shaft Geometry and Design Constraints
        self.L = design[0,:]  # section lengths of the shafts
        self.L = np.insert(self.L,0,0)  # correspond length 1 with L[1], etc.
        self.D = design[1,:]  # section diameters of the shafts
        self.D = np.insert(self.D,0,0)
        self.sf_o = design[2,0]  # safety factor against yielding
        self.sf_s = design[2,1]  # stress safety factor against fatigue
        self.life = design[2,2]  # desired service life in days
        self.twist_max = design[2,3]  # maximum allowable angle of twist
        self.R = design[2,4:7]  # radii of the top, bottom, crankshaft gear
        self.R = np.insert(self.R,0,0)
        self.r_prop = 1.5  # radius of the propeller
        
        L = self.L
        self.distance = {}
        # Calculate the distances from the left edge of the main shaft.
        # 'distance_1' is the distance of the end of L1 from the left edge, etc.
        self.distance[0] = 0
        self.distance[1] = L[1] 
        self.distance[2] = L[2]+L[1] 
        self.distance[3] = L[3]+L[2]+L[1] 
        self.distance[4] = L[4]+L[3]+L[2]+L[1] 
            
        # Calculate distances from the left edge of the crankshaft.
        self.distance[5] = L[5]
        self.distance[6] = L[6]+L[5]
        self.distance[7] = L[7]+L[6]+L[5]

        self.length = {}
        # Create arrays of points across the shaft sections. Each point
        # represents the distance from the left edge and is used to 
        # calculate stresses, etc, at that point.
        self.length[0] = np.linspace(0, 0, 100)
        self.length[1] = np.linspace(0, self.distance[1], 100)
        self.length[2] = np.linspace(self.distance[1], self.distance[2], 100)
        self.length[3] = np.linspace(self.distance[2], self.distance[3], 100)
        self.length[23] = np.linspace(self.distance[1], self.distance[3], 100)
        self.length[4] = np.linspace(self.distance[3], self.distance[4], 100)
        self.length[5] = np.linspace(0, self.distance[5], 100)
        self.length[6] = np.linspace(self.distance[5], self.distance[6], 100)
        self.length[7] = np.linspace(self.distance[6], self.distance[7], 100)

    def open_prop(self, file):
        mat_prop = np.loadtxt(file)
        
        # Material Properties
        self.mass_density = mat_prop[0]  # material density
        self.young_E = 10**9*mat_prop[1]  # Young's modulus
        self.poisson_v = mat_prop[2]  # Poisson's ratio
        self.sigma_o = 10**6*mat_prop[3]  # yield strength
        self.sigma_f = 10**6*mat_prop[4]  # Basquin coefficient
        self.exp_b = mat_prop[5]  # Basquin exponent

    def open_load(self, file):
        load = np.loadtxt(file)
        
        # Load Factors
        self.lamb = load[0]  # ratio of propeller tip to wind speed (linear)
        self.rho_air = load[1]  # air density
        self.cp = load[2]  # coefficient of power
        self.theta = load[3]  # angle of force acting on the crankshaft

    def open_wind(self, file):
        wind_history = np.loadtxt(file)
        
        # Wind History
        self.v_wind = wind_history[:,0]  # velocities of the wind
        self.v_wind = np.insert(self.v_wind,0,0)
        self.duration = 3600*wind_history[:,1]  # durations of the wind
        self.duration = np.insert(self.duration,0,0)