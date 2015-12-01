"""This module includes statics and mechanics formulas. 
All functions use values in SI units.
"""

from __future__ import division
from math import pi,sqrt

def calc_mass_cyl(length, diameter, density):
    """Return the mass of a cylinder.
    
    Args:
        length: The length of the cylinder.
        diameter: The diameter of the cylinder.
        density: The density of the cylinder's material.
        
    Returns:
        The mass of the cylinder.
    """
    mass = density*pi*length*(diameter/2)**2
    return mass
    
def calc_magnitude(a, b):
    """Return the magnitude of two orthogonal vectors
    
    Args:
        a = vector 1
        b = vector 2
        
    Returns:
        The magnitude of a + b, analagous to c in Pythagoras' theorem.
    """
    c = (a**2 + b**2)**0.5
    return c
      
def calc_axial_stress(force, diameter):
    """Return the compressive normal stress faced by a cylinder due to
    an axial load.
    
    Args:
        force: The axial force acting on the cylinder.
        diameter: The diameter of the cylinder.
        
    Returns:
        The normal stress faced by the cylinder.
    """
    sigma = force/(pi*(diameter/2)**2)
    return sigma

def calc_bend_stress(moment, diameter):
    """Return the normal stress faced by a cylinder due to bending.
    
    Args:
        moment: The moment acting on the cylinder.
        diameter: The diameter of the cylinder.
        
    Returns:
        The normal stress faced by the cylinder.
    """
    I = (pi/4)*(diameter/2)**4  # first moment of inertia
    sigma = (moment*(diameter/2))/I
    return sigma

def calc_octa_stress(sigma, torque, diameter):
    """Return the effective/octahedral stress faced by a cylinder 
    according to the von Mises yielding criterion.
    
    Args:
        sigma: The normal stress faced by the cylinder.
        torque: The torque acting on the cylinder.
        diameter: The diameter of the cylinder.
        
    Returns:
        The effective stress acting on a cylinder.
    """
    J = (pi/2)*(diameter/2)**4  # second moment of inertia
    tau = (torque*(diameter/2))/J  # shear stress due to torsion
    octa = (1/sqrt(2))*(2*sigma**2+6*tau**2)**0.5
    return octa

def calc_ECRSA(stress_amplitude, mean_stress, basquin_coef):
    """Return the equivalent completely reversed stress amplitude.
    
    The value returned by this function is required to calculate how
    quickly a material fails from cyclic loading.
    
    Args:
        stress_amplitude: The maximum stress faced by the material.
        mean_stress: The average stress faced.
        basquin_coef: The Basquin coefficient, a material property.
        
    Returns:
        The equivalent completely reversed stress amplitude.
    """
    ecrsa = stress_amplitude/(1-(mean_stress/basquin_coef))
    return ecrsa
	
def calc_cycles_fail(ecrsa, basquin_coef, exp_b):
    """Return the number of cycles until a material fails.
    
    Args:
        ecrsa: The equivalent completely reversed stress amplitude.
        basquin_coef: A The Basquin coefficient, a material property.
        exp_b: The Basquin exponent, a material property.
        
    Returns:
        The number of cycles until a material fails.
    """
    if ecrsa == 0:
        n_fail = 'i'  # infinite due to no cyclic loading.
    else:
        n_fail = (1/2)*(ecrsa/basquin_coef)**(1/exp_b)
    return n_fail
    
def calc_cycle_ratio(given_cycles, failure_cycles):
    """Return the ratio of total cycles to cycles until a material
    fails.
    
    Args:
        given_cycles: The number of cycles faced.
        failure_cycles: The number of cycles until failure.
        
    Returns:
        The ratio of cycles to cycles until failure.
    """
    if failure_cycles == 'i':
        cycle_ratio = 0
    else:
        cycle_ratio = given_cycles/failure_cycles
    return cycle_ratio
	
def calc_sf_yield(yield_strength, stress):
    """Return the safety factor against yielding.
    
    Args:
        yield_strength: The stress required to cause yielding.
        stress: The stress acting on the material.
        
    Returns:
        The  material's safety factor against yielding.
    """
    sf_y = yield_strength/stress
    return sf_y
    
def calc_sf_fatigue(days_to_failure, service_life_days, exp_b):
    """Return the safety factor against fatigue.
    
    Args:
        days_to_failure: The number of days until material failure.
        service_life_days: The desired service life.
        exp_b: The Basquin exponent, a material property.
        
    Returns:
        The safety factor against fatigue from to cyclic loading.
    """
    sf_fatigue = (days_to_failure/service_life_days)**(-exp_b)
    return sf_fatigue
    
def calc_max_twist(torque, length, diameter, young, poisson):
    """Return the maximum angle of twist of a cylinder.
    
    Args:
        torque: The torque acting on the cylinder.
        length: The length of the cylinder.
        diameter: The diameter of the cylinder.
        young: Young's modulus.
        poisson: Poisson's ratio.
        
    Returns:
        The angle of twist faced by the cylinder.
    """
    J = (pi/2)*(diameter/2)**4  # second moment of inertia
    G = young/(2*(1+poisson))  # modulus of rigidity
    twist = (torque*length)/(J*G)
    return twist

def calc_life(sum_n):
    """Return the number of days a material will last before failing.
    
    Args:
        sum_n: The sum of ratios of cycles faced by a material vs
        cycles until failure.
        
    Returns:
        The number of days a material will last before failing.
    """
    life = 1/sum_n
    return life