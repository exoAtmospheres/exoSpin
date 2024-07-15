'''
ExoSpin - Obliquity Run Script

This run script is a way to get the exoplanet obliquity from several parameters of the exoplanet.
The function obliquity() does the run script.


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports
import sys, os, glob
import numpy as np
from scipy.interpolate import interp1d
from exoplanet_class import *


def obliquity():
    """
    Compute the obliquity from exoplanet parameters and plot the obliquity.

    Returns:
        (Exoplanet): an Exoplanet class object.

    """
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
    radius_samp = np.loadtxt(radius_file, skiprows=1,usecols=(1,))
    vsini_input= input('Where are your data for rotational velocity (km/s unit)? ')
    vsini_file = open(vsini_input, "r")
    vsini_samp = np.loadtxt(vsini_file, skiprows=1,usecols=(2,))
    omega_o_input = input('Where are your data for sky projected inclination (° unit)? ')
    omega_o_file = open(omega_o_input, "r")
    omega_o_samp = np.loadtxt(omega_o_file, skiprows=1)
    
    P = input('What is the period of the exoplanet? (h unit)? ')
    P = float(P)
    M = input('What is the mass of the exoplanet? (M_Jup unit)? ')
    M = float(M)

    exoplanet = Exoplanet(exoplanet_name, io_samp, radius_samp, vsini_samp, omega_o_samp, P, M) 

    print('-> ExoSpin Computing')
    print()


    # ------------------------------------------------------------------------------------------------------------------------------------------------
    ## Computing ExoSpin

    a = input('Which method of computing do you want? (easy/complex) ')

    while a!='easy' and a!='complex':
        print()
        print('You need to choose a method of computing!')
        print()
        a = input('Which method of computing do you want? (easy/complex) ')

    if a == 'easy':
        print('Easy method computing ...')
    
    else :
        print('Complex method computing ...') 

    exoplanet.spin_axis_data()
    exoplanet.proj_obli_data()
    exoplanet.true_obli_data()


    # ------------------------------------------------------------------------------------------------------------------------------------------------
    ## Plot ExoSpin

    print()
    print('-> ExoSpin Plot')
    print()

    if a == "easy":
        obli_exoplanet = exoplanet.plot_obli('easy')

    else:
        obli_exoplanet = exoplanet.plot_obli('complex')

    return exoplanet

