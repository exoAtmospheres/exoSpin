'''
ExoSpin run script

This run script is a way to get the important informations for computing obliquity
from user's data.


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports
import sys, os, glob
import numpy as np
from scipy.interpolate import interp1d
from class_exoplanet import *

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Initialize ExoSpin
print()
print('Initializing ExoSpin ...' )
print()
print('-> ExoSpin Configuration')
print()
io_input = input('Where are your data for orbital inclination (° unit)? ')
io_file = open(io_input, "r")
io_samp = np.loadtxt(io_file, skiprows=1)
radius_input= input('Where are your data for radius (R_Jup unit)? ')
radius_file = open(radius_input, "r")
radius_samp = np.loadtxt(radius_file, skiprows=1)
vsini_input= input('Where are your data for rotational velocity (km/s unit)? ')
vsini_file = open(vsini_input, "r")
vsini_samp = np.loadtxt(vsini_file, skiprows=1)
omega_o_input = input('Where are your data for star inclination (° unit)? ')
omega_o_file = open(omega_o_input, "r")
omega_o_samp = np.loadtxt(omega_o_file, skiprows=1)
P = input('What is the period of the exoplanet? (h unit)? ')
M = input('What is the mass of the exoplanet? (M_Jup unit)? ')

print()
print('ExoSpin Computing ...')
print()

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Rearrange Data


# ------------------------------------------------------------------------------------------------------------------------------------------------
## Run ExoSpin
print()
print('Which method of computing do you want? (easy/complex)')

if input()=='easy':
    print()
    print('... computing ...')
    #exoplanet_input.easy_obliquity()
else :
    print()
    print('... computing ...')
    #exoplanet_input.complex_obliquity

print()
print('Computing done ^^ !')
print()

