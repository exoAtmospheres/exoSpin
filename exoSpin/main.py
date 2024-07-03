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

exoplanet_name = input('What\'s the name of your beautiful exoplanet? ')
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

exoplanet = Exoplanet(exoplanet_name, io_samp, radius_samp, vsini_samp, omega_o_samp, P, M) 

print()
print('ExoSpin Computing ...')
print()

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Computing ExoSpin
print()

a = input('Which method of computing do you want? (easy/complex) ')

while a!='easy' and a!='complex':
    print()
    print('You need to choose a method of computing!')
    print()
    a = input('Which method of computing do you want? (easy/complex) ')

if a=='easy':
    print()
    print('... easy method computing ...')
    print()
    print('... spin axis computing ...')
    print()
    print('... obliquity computing ...')
elif a=='complex' :
    print()
    print('... complex method computing ...')
    print()
    print('... spin axis computing ...')
    print()
    print('... obliquity computing ...')

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Plot ExoSpin

print()
print('There is the plot for the obliquity of ' + exoplanet.planet_name)
print()

# exoplanet/data/rp_vrot.txt