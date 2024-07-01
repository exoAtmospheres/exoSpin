'''
ExoSpin run script - Distrib

This run script is make sure that everything is set well (distribution size, limits, etc...)


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

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Order

#def order_data (io,vsini,radius,omega_o, v_limit, P):
