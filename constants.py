"""This module holds a class that loads constants from files.
All values are in SI units unless stated otherwise.
"""

import numpy as np

class Constants():

    def open_design(self, file):
        """Load the wind turbine design specifications.
        
        L: lengths of the shafts' sections.
        D: diameters of the shafts' sections.
        R: radii of the top, bottom, and crank gears.
        sf_o: safety factor against yielding.
        sf_s: safety factor against stress-based fatigue.
        life: desired service life in days.
        twist_max: maximum allowable angle of twist.
        r_prop: radius of the propeller.
        """
        design = np.loadtxt(file)
        
        # Shaft Geometry and Design Constraints
        self.L = np.array([None])  # correspond length 1 with L[1], etc.
        self.L = np.append(self.L, design[0,:])
        self.D = np.array([None])
        self.D = np.append(self.D, design[1,:])
        self.R = np.array([None])
        self.R = np.append(self.R, design[2,4:7])
        self.sf_o = design[2,0]
        self.sf_s = design[2,1]
        self.life = design[2,2]
        self.twist_max = design[2,3]
        self.r_prop = 1.5
        
        L = self.L
        self.distance = {}
        # Calculate the distances from the left edge of the main shaft.
        # 'distance_1' is the distance from the left edge to the end of L1.
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
        """Load the material properties.
        
        mass_denesity: material density.
        young_E: Young's modulus.
        poisson_v: Poisson's ratio.
        sigma_o: yield strength.
        sigma_f: Basquin coefficient.
        exp_b: Basquin exponent.
        """
        mat_prop = np.loadtxt(file)
        
        # Material Properties
        self.mass_density = mat_prop[0]
        self.young_E = 10**9*mat_prop[1]
        self.poisson_v = mat_prop[2]
        self.sigma_o = 10**6*mat_prop[3]
        self.sigma_f = 10**6*mat_prop[4]
        self.exp_b = mat_prop[5]

    def open_load(self, file):
        """Load the load factors.
        
        lamb: ratio of linear speeds of propeller tip to wind.
        rho_air: density of the incoming air.
        cp: coefficient of power.
        theta: force angle (relative to the horizontal) on the crankshaft.
        """
        load = np.loadtxt(file)
        
        # Load Factors
        self.lamb = load[0]
        self.rho_air = load[1]
        self.cp = load[2]
        self.theta = load[3]

    def open_wind(self, file):
        """Load the wind history.
        
        v_wind: velocities of the wind.
        duration: durations of corresponding wind velocities.
        """
        wind_history = np.loadtxt(file)
        
        # Wind History
        self.v_wind = np.array([None])
        self.duration = np.array([None])

        try:
            self.v_wind = np.append(self.v_wind, wind_history[:,0])
            self.duration = np.append(self.duration, 3600*wind_history[:,1])
        # If file contains 1 line, do not call for >1 rows.
        except IndexError:
            self.v_wind = np.append(self.v_wind, wind_history[0])
            self.duration = np.append(self.duration, wind_history[1])
        