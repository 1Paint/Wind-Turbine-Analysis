"""This project builds off a project given by Cornell's Spring 2013 MAE
2120 course. Given external settings, material properties, and the
design for a wind turbine, students were asked to solve for and
summarize mechanical factors faced by sections of the turbine's shafts.

Originally coded as a single, standalone file in MATLAB, this program
utilizes object-oriented programming via Python and integrates multiple
libraries in order to solve for, plot, display, and iterate via a GUI
turbine shaft geometries and safety factors.
"""
from __future__ import division

import sys
import numpy as np

from math import pi, sqrt, cos, sin
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from PyQt4 import QtGui, QtCore

from formulas import *
from constants import Constants

class Solutions(object):
    """This class houses values of safety factors and shaft specifications."""
    def __init__(self, constants):
        self.c = constants
        count = len(self.c.L)
        
        self.mass = [None]*count
        self.sf_yield = [None]*count
        self.sf_fatigue = [None]*count
        self.max_twists = [None]*count
        self.life = [None]*count
        self.values_set = [None]*len(self.c.v_wind)
     
    def solve(self):
        """Solve for safety factors against fatigue and yielding, the
        maximum angle of twist, and days until failure for the shafts.
        """
        c = self.c
        # Check for invalid geometry. Run if all dimension are valid.
        if 0 in c.L[1:]:
            print "You cannot have lengths of 0!"
        elif 0 in c.D[1:]:
            print "You cannot have diameters of 0!"
        else:
            # Store relevant data for sections of the main shaft and crankshaft
            # under wind forces of varietal magnitudes and durations.
            for i in range(1, len(c.v_wind)):
                self.values_set[i] = self.get_values(c.v_wind[i], 
                                                        c.duration[i], i, c)
                
            # Find the set of values corresponding to the strongest wind force
            # and store corresponding safety factors against yielding and the
            # maximum angles of twist.
            max_wind_index = np.where(c.v_wind==max(c.v_wind))[0][0]
            self.sf_yield[1:] = self.values_set[max_wind_index][7:14]
            self.max_twists[1:] = self.values_set[max_wind_index][14:21]

                
            #######################
            # Days Until Failures #
            #######################

            # Calculate the days until sections of the shaft fail under given
            # wind conditions.
            sum_n = [0]*len(c.L)
            for i in range(1, len(c.L)):  # for each section
                for k in range(1, len(c.v_wind)):  # for each wind
                    sum_n[i] += self.values_set[k][i-1]
                if sum_n[i] == 0:
                    self.life[i] = 'i' # inf
                else:
                    self.life[i] = calc_life(sum_n[i])

            #################################
            # Safety Factor Against Fatigue #
            #################################
                    
            # Calculate the safety factor against fatigue for sections of the
            # shafts.
            for i in range(1, len(c.L)):
                if self.life[i] == 'i':
                    self.sf_fatigue[i] = 'i' # inf
                else:
                    self.sf_fatigue[i] = calc_sf_fatigue(self.life[i],
                                                         c.life, c.exp_b)

            ########
            # Mass #
            ########

            # Create an array containing the masses of sections L1-L7 of the
            # shafts.
            for i in range(1, len(c.L)):
                self.mass[i] = (calc_mass_cyl(c.L[i], c.D[i], c.mass_density))

    def get_values(self, velocity, time, i, c):
        """Return a list containing safety factors against yielding,
        maximum angles of twist, and values required to solve for 
        safety factors against fatigue and days until failure.
        
        This method analyzes the shafts for a given wind velocity and 
        duration and obtains the maximum angles of twist undergone by 
        sections of the shaft, as well as safety factors against 
        yielding and the ratios of shaft revolutions to revolutions
        until failure.
        
        Args:
            velocity: The wind velocity.
            time: The corresponding wind duration.
            i: An index marking the wind history.
            c: Constants loaded from files.
            
        Returns:
            A list containing ratios of cycles to cycles until failure,
            safety factors against yielding, and the maximum angle of
            twist faced by the 7 shaft sections.
        """
        # power generated by the wind
        p_wind = c.cp*0.5*c.rho_air*(pi*c.r_prop**2)*velocity**3  

        # tangential velocity of the propeller tip
        u = c.lamb*velocity  
        
        if velocity != 0:
            # force_+x of the wind on the main shaft
            f_wind = p_wind/velocity
            
            # angular velocity of the propeller
            omega_prop = u/c.r_prop
        
            # torque_+x acting on the main shaft 
            t_in = p_wind/omega_prop
            
        # 0 velocity results in 0 force, 0 rotation, and 0 torque.
        else:
            f_wind = 0
            omega_prop = 0
            t_in = 0

        # torque_+x acting on the main shaft
        t_gear_main = -t_in

        # force_+y acting on the right end of the main shaft
        f_gear_main = t_gear_main/c.R[1]

        # reaction forces at point 2
        rxn_2_x = -f_wind
        rxn_2_y = -((c.L[2]+c.L[3]+c.L[4])*f_gear_main)/(c.L[2]+c.L[3])
        
        # reaction forces at point 1
        rxn_1_y = -rxn_2_y-f_gear_main

        # torque acting on the gear of crankshaft
        t_gear_crank = -t_gear_main*(c.R[2]/c.R[1])

        # force_+y acting on the gear of the crankshaft
        f_gear_crank = -f_gear_main

        # torque acting on the crankshaft by the pump rod
        t_crank = -t_gear_crank
        
        # Calculate the force generated by the pump rod on the crankshaft
        # and decompose it into orthogonal forces.
        f_pump = abs(t_crank/(c.R[3]*cos(c.theta*pi/180)))
        f_pump_z = -f_pump*cos(c.theta*pi/180)
        f_pump_y = f_pump*sin(c.theta*pi/180)
        
        # reaction forces at point 3
        rxn_3_y = ((c.L[7]*f_gear_crank)-(c.L[6]*f_pump_y))/(c.L[6]+c.L[5])
        rxn_3_z = (-f_pump_z*c.L[6])/(c.L[5]+c.L[6])

        # reaction forces at point 4
        rxn_4_y = -f_gear_crank-rxn_3_y-f_pump_y
        rxn_4_z = -f_pump_z-rxn_3_z
         
        ###################
        # Bending Moments #
        ###################
        
        bending_z = {}
        # Calculate z+ moments across the main shaft (lengths 1, 2&3, 4).
        bending_z[1] = 0*c.length[1]
        bending_z[2] = bending_z[1][-1]+(rxn_1_y*(c.length[2]-c.distance[1]))
        bending_z[3] = bending_z[2][-1]+(rxn_1_y*(c.length[3]-c.distance[2]))
        bending_z[23] = bending_z[1][-1]+(rxn_1_y*(c.length[23]-c.distance[1]))
        bending_z[4] = bending_z[23][-1]+((rxn_2_y+rxn_1_y)
                                            *(c.length[4]-c.distance[3]))
        # Calculate z+ bending moments across the crankshaft.
        bending_z[5] = rxn_3_y*c.length[5]
        bending_z[6] = bending_z[5][-1]+((rxn_3_y+f_pump_y)
                                            *(c.length[6]-c.distance[5]))
        bending_z[7] = bending_z[6][-1]+((rxn_3_y+f_pump_y+rxn_4_y)
                                            *(c.length[7]-c.distance[6]))

        bending_y = {}
        # Calculate y+ bending moments across the crankshaft.
        bending_y[5] = rxn_3_z*c.length[5]
        bending_y[6] = bending_y[5][-1]+((rxn_3_z+f_pump_z)
                                            *(c.length[6]-c.distance[5]))
        bending_y[7] = bending_y[6][-1]+((rxn_3_z+f_pump_z+rxn_4_z)
                                            *(c.length[7]-c.distance[6]))

        bending_mag = {}
        # Calculate the magnitudes of bending moments across the crankshaft.
        for i in range(5, 8):
            bending_mag[i] = calc_magnitude(bending_z[i], bending_y[i])

        #########################
        # Normal Stress (Axial) #
        #########################
        
        sigma_axial = {}
        # Calculate the normal stresses from axial loads across the main shaft.
        for i in range (1, 5):
            sigma_axial[i] = calc_axial_stress(-f_wind, c.D[i])

        ###########################
        # Normal Stress (Bending) #
        ###########################

        sigma_bending_x = {}
        # Calculate the normal stresses from bending across the main shaft.
        for i in range(1, 5):
            sigma_bending_x[i] = abs(calc_bend_stress(bending_z[i], c.D[i]))
        # Calculate the normal stresses from bending across the crankshaft.
        for i in range(5, 8):
            sigma_bending_x[i] = abs(calc_bend_stress(bending_mag[i], c.D[i]))

        ####################
        # Effective Stress #
        ####################

        sigma_eff = {}
        # Calculate the effective stresses across the main shaft.
        for i in range(1, 5):
            sigma_eff[i] = calc_octa_stress(
                                sigma_bending_x[i]+abs(sigma_axial[i]),
                                t_in, c.D[i])
        # Calculate the effective stresses across the crankshaft.
        sigma_eff[5] = calc_octa_stress(sigma_bending_x[5], 0, c.D[5])
        for i in range(6, 8):
            sigma_eff[i] = calc_octa_stress(
                                sigma_bending_x[i],
                                t_gear_crank, c.D[i])

        ##################################
        # Safety Factor against Yielding #
        ##################################

        sf_yield = {}
        # Calculate safety factors across both shafts.
        for i in range(1, 8):
            sf_yield[i] = calc_sf_yield(c.sigma_o, max(sigma_eff[i]))

        ######################
        # Mean Cyclic Stress #
        ######################

        sigma_avg = {}
        # Store the mean cyclic stress across the main shaft.
        for i in range(1, 5):
            sigma_avg[i] = sigma_axial[i]
        # Store the mean cyclic stress across the crankshaft
        for i in range(5, 8):
            sigma_avg[i] = 0  # due to no axial loads--cycles around 0 stress

        ####################
        # Stress Amplitude #
        ####################

        sigma_amp = {}
        # Obtain the amplitudes of cyclic stresses.
        for i in range(1, 8):
            sigma_amp[i] = max(sigma_bending_x[i])

        ###################################################
        # Equivalent Completely Reversed Stress Amplitude #
        ###################################################

        sigma_rev = {}
        # Calculate the equivalent completely reversed stress amplitudes across
        # the shafts.
        for i in range(1, 8):
            sigma_rev[i] = calc_ECRSA(sigma_amp[i], sigma_avg[i], c.sigma_f)

        ########################
        # Cycles until Failure #
        ########################

        n_failure = {}
        # Calculate the cycles to failure for the shaft sections.
        for i in range(1, 8):
            n_failure[i] = calc_cycles_fail(sigma_rev[i], c.sigma_f, c.exp_b)

        ###########################################
        # Ratio of Cycles to Cycles until Failure #
        ###########################################

        n_ratio = {}
        # Calculate the ratio for the main shaft.
        n_given_main = (omega_prop/(2*pi))*time  # number of cycles
        for i in range(1, 5):
            n_ratio[i] = calc_cycle_ratio(n_given_main, n_failure[i])
        # Calculate the ratio for the crankshaft.
        n_given_crank = n_given_main/(c.R[2]/c.R[1])  # number of cycles
        for i in range(5, 8):
            n_ratio[i] = calc_cycle_ratio(n_given_crank, n_failure[i])

        ##########################
        # Maximum Angle of Twist #
        ##########################

        twist_max = {}
        # Calculate the maximum angle of twist for the main shaft.
        for i in range(1, 5):
            twist_max[i] = max(abs(calc_max_twist(
                                    t_in, c.length[i]-c.distance[i-1],
                                    c.D[i], c.young_E, c.poisson_v)))
        # Calculate the maximum angle of twist for the crankshaft.
        twist_max[5] = max(abs(calc_max_twist(
                                    0, c.length[5],
                                    c.D[5], c.young_E, c.poisson_v)))
        for i in range(6, 8):
            twist_max[i] = max(abs(calc_max_twist(
                                    t_crank, c.length[i]-c.distance[i-1],
                                    c.D[i], c.young_E, c.poisson_v)))

        ##########
        # Values #
        ##########

        # Store the necessary information required to calculate the safety
        # factor against stress-based fatigue.
        answer = [n_ratio[1], n_ratio[2], n_ratio[3], n_ratio[4], 
                  n_ratio[5], n_ratio[6], n_ratio[7], 
                  sf_yield[1], sf_yield[2], sf_yield[3], sf_yield[4], 
                  sf_yield[5], sf_yield[6], sf_yield[7], 
                  twist_max[1], twist_max[2], twist_max[3], twist_max[4], 
                  twist_max[5], twist_max[6], twist_max[7]
                  ]
                              
        return answer

class InfoPlot(FigureCanvas):
    """This class takes inputted statics and mechanics solutions and
    plots them. Shown are the main shaft and crankshaft geometries and
    safety factors as well as design specifications.
    """
    def __init__(self, solution, constants, normalized, parent=None):
        self.solution = solution
            
        fig = Figure(tight_layout=True)
        FigureCanvas.__init__(self, fig)
        
        self.plot_info(fig, constants, normalized)

    def plot_info(self, fig, constants, normalized):
        """Plot the shafts and corresponding information.
        
        This method takes in a set of constants used to calculate where
        shafts and corresponding infomation are plotted. The value of
        'normalized' determines whether or not shaft information is
        plotted directly under the plots of the shafts.
        """
        c = constants
        
        # initialize shaft coordinates
        x = [0]*len(c.L)
        y = [0]*len(c.L)
        main_distance = 0
        crank_distance = 0
        
        # reference values used to scale plotted shafts and text
        main_D_max = max(c.D[1:5])
        crank_D_max = max(c.D[5:8])
        main_length = sum(c.L[1:5])
        crank_length = sum(c.L[5:8])
        
        # reference values for plotting info depending on viewing property
        if normalized == False:
            main_info_location = 0
            crank_info_location = 0
        else:
            main_info_location = -main_length/8
            crank_info_location = -crank_length/6
        
        # Draw each shaft section and display relevant info.
        for i in range(1,8):

            failure = False  # whether or not safety factors are exceeded
            
            # Plot-scaling values depend on what shaft is being plotted.
            if i <= 4:
                a = main_D_max
            else:
                a = crank_D_max
                
            # Create a subplot for the main shaft.
            if i == 1:
                self.axes = fig.add_subplot(211)
                self.axes.hold(True)
                
                shaft_distance = main_distance
                info_location = main_info_location

            # Create a subplot for the crankshaft.
            if i == 5:
                self.axes = fig.add_subplot(212)
                self.axes.hold(True)
                
                shaft_distance = crank_distance
                info_location = crank_info_location
                
            # Location of shaft info depends on viewing property.
            if normalized == False:
                info_location = shaft_distance+c.L[i]/2
            else:
                if i <= 4:
                    info_location += main_length/4
                else:
                    info_location += crank_length/3
            
            # Print shaft geometry onto the plot.
            self.axes.text(info_location, -(a/1.5)-(a/3.2),
                           'Diameter: %.2g m' % (c.D[i]),
                           ha='center', fontsize=11)
            self.axes.text(info_location, -(a/1.5)-2*(a/3.2),
                           'Length: %.2g m' % (c.L[i]),
                           ha='center', fontsize=11)
            
            # Print shaft safety factor against yielding.
            if self.solution.sf_yield[i] == 'i':
                sf_yield = 'inf'
                color = 'k'
            else:
                sf_yield = str('%.2f' % (self.solution.sf_yield[i]))
                if self.solution.sf_yield[i] < c.sf_o:
                    color = 'r'
                    failure = True
                else:
                    color = 'k'
            self.axes.text(info_location, -(a/1.5)-3*(a/3.2),
                           'Xo: %s' % (sf_yield), color=color,
                           ha='center', fontsize=11)
            
            # Print shaft safety factor against fatigue.
            if self.solution.sf_fatigue[i] == 'i':
                sf_fatigue = 'inf'
                color = 'k'
            else:
                sf_fatigue = str('%.2f' % (self.solution.sf_fatigue[i]))
                if self.solution.sf_fatigue[i] < c.sf_s:
                    color = 'r'
                    failure = True
                else:
                    color = 'k'
            self.axes.text(info_location, -(a/1.5)-4*(a/3.2),
                           'Xs: %s' % (sf_fatigue), color=color,
                           ha='center', fontsize=11)
            
            # Print shaft maximum angle of twist.
            if self.solution.max_twists[i] == 'i':
                twist = 'inf'
                color = 'r'
            else:
                twist = str('%.2g' % (self.solution.max_twists[i]))
                if self.solution.max_twists[i] > c.twist_max:
                    color = 'r'
                    failure = True
                else:
                    color = 'k'
            self.axes.text(info_location, -(a/1.5)-5*(a/3.2),
                           'Twist: %s rad' % (twist), color=color,
                           ha='center', fontsize=11)
            
            # Print shaft service life.
            if self.solution.life[i] == 'i':
                service = 'inf'
                color = 'k'
            else:
                service = str('%.3g' % (self.solution.life[i]))
                if self.solution.life[i] < c.life:
                    color = 'r'
                    failure = True
                else:
                    color = 'k'
            self.axes.text(info_location, -(a/1.5)-6*(a/3.2),
                           'Life: %s days' % (service), color=color,
                           ha='center', fontsize=11)
            
            # Draw the shafts. 
            if failure == True:
                color = 'r'
            else:
                color = 'k'
            x[i] = np.linspace(shaft_distance, shaft_distance+c.L[i], 2)
            y[i] = np.linspace(c.D[i]/2, c.D[i]/2, 2)
              
            # Plot horizontal lines
            self.axes.plot(x[i], y[i], '-', linewidth=2.0, color=color)
            self.axes.plot(x[i], -y[i], '-', linewidth=2.0, color=color)
            
            # Plot vertical lines
            self.axes.plot((shaft_distance, shaft_distance),
                            (-c.D[i]/2, c.D[i]/2),
                            '-', linewidth=2.0, color=color)
                            
            shaft_distance += c.L[i]
            
            self.axes.plot((shaft_distance, shaft_distance),
                            (-c.D[i]/2, c.D[i]/2),
                            '-', linewidth=2.0, color=color)
                            
            if i == 4:
                # Configure the plot axes and titles of the main shaft.
                self.configure_plot(self.axes, shaft_distance, a, c, i)
                
            if i == 7:
                # Configure the plot axes and titles of the crankshaft.
                self.configure_plot(self.axes, shaft_distance, a, c, i)
                            
    def configure_plot(self, axes, shaft_distance, a, c, i):
        """Configure the plot axes and titles and display specifications.
        
        This method takes in:
            axes: The subplot being made.
            shaft_distance: The length of the shaft.
            a: A value to scale text to the plot.
            c: A set of constants. Specifically mass is needed.
            i: The iteration number to determine what shaft is being plotted. 
        """
        xmin = 0
        xmax = shaft_distance
        ymin = 3*-a
        ymax = 3*a
        axes.axis([xmin, xmax, ymin, ymax])
        if i == 4:
            axes.set_title('Main Shaft')  
            mass = sum(self.solution.mass[1:5])
        else:
            axes.set_title('Crankshaft')
            mass = sum(self.solution.mass[5:8])
        axes.set_xlabel('Length from Left Edge [m]')
        axes.set_ylabel('Radius [m]')
        
        # Plot design specifications.
        shaft_info = ("Total Mass: %.2f kg\nMinimum Required Safety Factor"
                     " Against Yielding, Xo: %.2g\nMinimum Required Safety"
                     " Factor Against Fatigue, Xs: %.2g\nMaximum Allowable"
                     " Angle of Twist: %.2g rad\nDesired Service Life: %.2g"
                     " days" % (mass, c.sf_o, c.sf_s, c.twist_max, c.life))                    
        axes.text(xmax/2, a, shaft_info, ha='center', fontsize=11)

class EditWindow(QtGui.QDialog):
    """This is a window which allows the user to change the current 
    lengths and diameters of sections of the shafts. When dimensions
    are updated, the main window refreshes and shows info corresponding
    to the updated design.
    """
    def __init__(self, constants, parent=None):
        QtGui.QDialog.__init__(self)
        self.setModal(True)
        self.setFixedWidth(200)
        self.parent = parent

        self.c = constants
        self.L = [0]*len(self.c.L)
        self.D = [0]*len(self.c.D)
        
        vlayout = QtGui.QVBoxLayout() 
        line = [None]*8

        # Add fields to edit main shaft sections.
        main_header = QtGui.QLabel('Main Shaft')
        vlayout.addWidget(main_header)
        for i in range(1,5):
            line[i] = QtGui.QHBoxLayout()
            line[i].addWidget(QtGui.QLabel("L%i" % i))
            self.L[i] = QtGui.QLineEdit(str(self.c.L[i]))
            line[i].addWidget(self.L[i])
            line[i].addWidget(QtGui.QLabel("D%i" % i))
            self.D[i] = QtGui.QLineEdit(str(self.c.D[i]))
            line[i].addWidget(self.D[i])
            vlayout.addLayout(line[i])  
        vlayout.addWidget(QtGui.QLabel(""))

        # Add fields to edit crankshaft sections.
        crank_header = QtGui.QLabel('Crankshaft')
        vlayout.addWidget(crank_header)
        for i in range(5,8):
            line[i] = QtGui.QHBoxLayout()
            line[i].addWidget(QtGui.QLabel('L%i' % i))
            self.L[i] = QtGui.QLineEdit(str(self.c.L[i]))
            line[i].addWidget(self.L[i])
            line[i].addWidget(QtGui.QLabel('D%i' % i))
            self.D[i] = QtGui.QLineEdit(str(self.c.D[i]))
            line[i].addWidget(self.D[i])
            vlayout.addLayout(line[i]) 
        vlayout.addWidget(QtGui.QLabel(""))
        
        # Add update and cancel buttons.
        update_button = QtGui.QPushButton('Update')
        update_button.clicked.connect(self.update)
        cancel_button = QtGui.QPushButton('Cancel')
        cancel_button.clicked.connect(self.cancel)
        
        buttonBox = QtGui.QDialogButtonBox()
        buttonBox.addButton(update_button, QtGui.QDialogButtonBox.ActionRole)
        buttonBox.addButton(cancel_button, QtGui.QDialogButtonBox.ActionRole)
        
        button_line = QtGui.QHBoxLayout()
        button_line.addStretch(0)
        button_line.addWidget(buttonBox)
        button_line.addStretch(0)
        vlayout.addLayout(button_line)
        
        self.setLayout(vlayout)
    
    def update(self):
        """Update the plot in the main window."""
        # Change lengths and dimensions to those inputted.
        for i in range(1, len(self.c.L)):
            self.c.L[i] = float(self.L[i].text())
            self.c.D[i] = float(self.D[i].text())
        self.parent.reload = True
        self.parent.refresh()
        self.close()

    # Close the window.
    def cancel(self):
        self.close()
        
class AppWindow(QtGui.QMainWindow):
    """This is the main application window displaying the shaft
    geometries and relevant safety factors. A toolbar allows the 
    user to load & save files, edit the configuration, and alter views.
    """
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)
        self.setWindowTitle("Wind Turbine Shaft Analysis")
        self.setFixedSize(1280,720)
  
        self.create_actions()
        self.create_menu()
        
        # Setup the main application.
        self.main_widget = QtGui.QWidget(self)
        self.layout = QtGui.QVBoxLayout(self.main_widget)
        self.main_widget.setFocus()
        self.setCentralWidget(self.main_widget)
        
        # Initialize constants and statuses of the application.
        self.init_values()
        
        # Update the main window with solutions.
        self.update_solution()
        
        # Tell the user what files are still required.
        self.update_status()
                
    def init_values(self):
        """Initialize statuses of whether or not files have been
        loaded, whether shaft information displays are normalized,
        and whether or not the plot needs to be refreshed. Also
        initialize a class holding constants.
        """
        # Holds the states of whether the necessary files have been loaded. 
        # The 4 required files are: 
        # - designs specifications
        # - material properties
        # - load factors
        # - wind history.
        self.loaded = [False, False, False, False]
        
        # This value determines how details of the shaft sections are displayed.
        # A value of 'True' scales text uniformly; A value of 'False' places
        # text underneath plotted shaft sections.
        self.normalized = False
        
        # This value determines whether solutions have to be (re)solved.
        # Prevents recalculating solutions if user only wants to change views.
        self.reload = True
        
        # Initialize constants.
        self.c = Constants()
        
    def create_actions(self):
        """Create menu actions."""
        # File Menu Actions
        self.act_obt_design = QtGui.QAction('Obtain Design', self,
                shortcut='Ctrl+1', triggered=self.open_design)     
        self.act_obt_prop = QtGui.QAction('Obtain Material Properties', self, 
                shortcut='Ctrl+2', triggered=self.open_prop)      
        self.act_obt_load = QtGui.QAction('Obtain Load Factors', self,
                shortcut='Ctrl+3', triggered=self.open_load)       
        self.act_obt_wind = QtGui.QAction('Obtain Wind History', self,
                shortcut='Ctrl+4', triggered=self.open_wind)
        self.act_save_design = QtGui.QAction('Save Design', self,
                shortcut='Ctrl+S', triggered=self.save_design)
        self.act_quit = QtGui.QAction('Quit', self,
                shortcut='Ctrl+Q', triggered=self.quit_file)        

        # Edit Menu Actions
        self.act_edit_design = QtGui.QAction('Edit Design', self,
                shortcut='Ctrl+E', triggered=self.edit_design)     
        
        # View Menu Actions
        self.act_normalize_view = QtGui.QAction('Normalize View', self,
                                            shortcut='Ctrl+V', checkable=True,
                                            triggered=self.normalize_view)
        self.act_normalize_view.setStatusTip('Uniformly reposition shaft '
                                'information. The order of information still '
                                'corresponds to the order of shaft sections '
                                'i.e. the first set of information corresponds'
                                ' to section 1, the second to section 2, etc.')
        
        # Help Menu Actions
        self.act_about = QtGui.QAction('About', self,
                shortcut='Ctrl+A', triggered=self.about)
                
    def create_menu(self):
        """Create the application menu with File, Edit, View, Help
        and About tabs.
        """
        file_menu = QtGui.QMenu('File', self)
        file_menu.addAction(self.act_obt_design)
        file_menu.addAction(self.act_obt_prop)
        file_menu.addAction(self.act_obt_load)
        file_menu.addAction(self.act_obt_wind)
        file_menu.addSeparator()
        file_menu.addAction(self.act_save_design)
        file_menu.addSeparator()
        file_menu.addAction(self.act_quit)
        
        edit_menu = QtGui.QMenu('Edit', self)
        edit_menu.addAction(self.act_edit_design)
        
        view_menu = QtGui.QMenu('View', self)
        view_menu.addAction(self.act_normalize_view)
        
        help_menu = QtGui.QMenu('Help', self)
        help_menu.addAction(self.act_about)
        
        menubar = self.menuBar()
        menubar.addMenu(file_menu)
        menubar.addMenu(edit_menu)
        menubar.addMenu(view_menu)
        menubar.addMenu(help_menu)
        
    def update_solution(self):
        """Plot shaft sections and related info to the application
        window.
        """
        if self.loaded == [True, True, True, True]:
            if self.reload == True:
                self.sol = Solutions(self.c)
                self.sol.solve()
                infoplot = InfoPlot(self.sol, self.c, 
                                    self.normalized, self.main_widget)
                self.layout.addWidget(infoplot)
                self.reload = False
            else:
                infoplot = InfoPlot(self.sol, self.c, 
                                    self.normalized, self.main_widget)
                self.layout.addWidget(infoplot)
        
    def delete_solution(self):
        """Delete the current plot of the shaft and infomation."""
        for i in range(self.layout.count()): 
            self.layout.itemAt(i).widget().setParent(None)
        
    def normalize_view(self):
        """Scale the info dislayed beneath each shaft section so that
        text does not interfer with each other or the plot.
        """
        self.normalized = not self.normalized
        self.refresh()
        
    def refresh(self):
        """Refresh the plot by deleting current solutions and updating
        with another set of solutions.
        """
        self.delete_solution()
        self.update_solution()
        self.update_status()
            
    def update_status(self):
        """Update the status bar and display what files are still
        required.
        """
        self.statusBar().showMessage(self.check_loaded())
        
    def status_format_error(self):
        """Display an error message on the status bar if the format of
        a file is not usable."""
        self.statusBar().showMessage("ERROR. Please check formatting of file.")
        
    def open_design(self):
        """Load shaft geometries from a file."""
        try:
            file = QtGui.QFileDialog.getOpenFileName(self, 'Open design')
            self.c.open_design(str(file))
            self.loaded[0] = True
            self.reload = True
            self.refresh()
        except IndexError or ValueError:
            self.status_format_error()
        
    def open_prop(self):
        """Load material properties from a file."""
        try:
            file = QtGui.QFileDialog.getOpenFileName(self, 'Open material properties')
            self.c.open_prop(str(file))
            self.loaded[1] = True
            self.reload = True
            self.refresh()
        except IndexError or ValueError:
            self.status_format_error()
        
    def open_load(self):
        """Load external/environmental factors from a file."""
        try:
            file = QtGui.QFileDialog.getOpenFileName(self, 'Open load factors')
            self.c.open_load(str(file))
            self.loaded[2] = True
            self.reload = True
            self.refresh()
        except IndexError or ValueError:
            self.status_format_error()
        
    def open_wind(self):
        """Load wind history from a file."""
        #try:
        file = QtGui.QFileDialog.getOpenFileName(self, 'Open wind history')
        self.c.open_wind(str(file))
        self.loaded[3] = True
        self.reload = True
        self.refresh()
        #except IndexError or ValueError:
            #self.status_format_error()
    
    def check_loaded(self):
        """Check to see what files have been loaded and return a
        message stating which files are still needed.
        """
        if False in self.loaded:
            msg0 = "Required: "
        else:
            msg0 = ""
            
        if self.loaded[0] == False:
            msg1 = "Design specifications. "
        else:
            msg1 = ""
            
        if self.loaded[1] == False:
            msg2 = "Material properties. "
        else:
            msg2 = ""
            
        if self.loaded[2] == False:
            msg3 = "Load factors. "
        else:
            msg3 = ""
            
        if self.loaded[3] == False:
            msg4 = "Wind history."
        else:
            msg4 = ""
            
        message = msg0 + msg1 + msg2 + msg3 + msg4
        return message
        
    def save_design(self):
        """Save the current design by storing shaft lengths and
        diameters into a file. Original design specifications are
        stored in order to preserve the design file format.
        """
        if self.loaded[0] == False:
            self.statusBar().showMessage("Please load design specifications first.")
        else:
            c = self.c
            design = (('%g %g %g %g %g %g %g\n'
                      '%g %g %g %g %g %g %g\n'
                      '%g %g %g %g %g %g %g') %
            (c.L[1], c.L[2], c.L[3], c.L[4], c.L[5], c.L[6], c.L[7],
             c.D[1], c.D[2], c.D[3], c.D[4], c.D[5], c.D[6], c.D[7],
             c.sf_o, c.sf_s, c.life, c.twist_max, c.R[1], c.R[2], c.R[3]))
            
            file = QtGui.QFileDialog.getSaveFileName(self, 'Save design',   
                                                        selectedFilter='*.txt')
            fname = open(file, 'w')
            fname.write(design)
            fname.close()
            
    def edit_design(self):
        """Open a new window to allow users to edit the lengths and
        diameters of the shaft sections.
        """
        if self.loaded[0] == False:
            self.statusBar().showMessage("Please load design specifications first.")
        else:
            self.window = EditWindow(self.c, parent=self)
            self.window.show()
            self.window.setWindowTitle('Edit Design [m]')

    def about(self):
        QtGui.QMessageBox.about(self, 'About', "This project builds off a"
                                " project given by Cornell's Spring 2013 MAE"
                                " 2120 course. Given external settings,"
                                " material properties, and the design for a"
                                " wind turbine, students were asked to solve"
                                " for and summarize mechanical factors faced"
                                " by sections of the turbine's shafts."
                                " Originally coded as a single, standalone"
                                " file in MATLAB, this program utilizes"
                                " object-oriented programming via Python and"
                                " integrates multiple libraries in order to"
                                " solve for, plot, display, and iterate via a"
                                " GUI turbine shaft geometries and safety"
                                " factors.")

    def quit_file(self):
        self.close()

    def closeEvent(self, ce):
        self.quit_file()
      
if __name__ == "__main__":

    App = QtGui.QApplication(sys.argv)
    
    aw = AppWindow()
    aw.setWindowTitle('Wind Turbine Shaft Analysis')
    aw.show()
    sys.exit(App.exec_())


    

    
    
