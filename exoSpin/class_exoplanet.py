'''
ExoSpin run script - Exoplanet Class


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
import astropy.constants as c
import pickle

from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.stats import gaussian_kde
from scipy.stats import uniform
from scipy.integrate import quad

# --------------------------------------------------------------------------------------------------------------------------------------------------

class Exoplanet():

    def __init__(self, io, radius, vsini, omega_o,period,mass):
        self.io = io
        self.radius = radius
        self.vsini = vsini
        self.omega_o = omega_o
        self.mass = mass
        self.period = period

    #def get_orbital_inclination(self):

    #def get_radius(self):

    #def get_rot_velocity(self):

    #def get_start_inclination(self):

    #def set_orbital_inclination(self, io):

    #def set_radius(self, radius):

    #def set_rot_velocity(self, vsini):

    #def set_start_inclination(self, omega_o):

    #def easy_obliquity(self):

    #def complex_obliquity(self):

    #def get_true_obliquity(self):

    #def get_projected_obliquity(self):

    #def plot_all(self):


        

    