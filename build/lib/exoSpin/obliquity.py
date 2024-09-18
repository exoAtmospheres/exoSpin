'''
ExoSpin - Obliquity Run Script

This run script is a way to get the exoplanet obliquity from several parameters of the exoplanet.
The function obliquity() does the run script.


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports

import numpy as np
from .exoplanet_class import *


def obliquity(exoplanet_name, io, radius, vsini, omega_o, P, M):
    """
    From exoplanet data, the function computes the obliquity of the planet and returns an Exoplanet obkect.

    Args:
        exoplanet_name (str): Planet's name.
        io (str or list): Path for the orbital inclination data file.
                                    It can be a path or a list where the first element is the mean value, the second the standard deviation, and the third the length of it. From it, a random distribution will be generated.
        radius (str or list): Path for the radius data file.
                                    It can be a path or a list where the first element is the mean value, the second the standard deviation, and the third the length of it. From it, a random distribution will be generated.
        vsini (str or list): Path for the rotational velocity data file.
                                    It can be a path or a list where the first element is the mean value, the second the standard deviation, and the third the length of it. From it, a random distribution will be generated.
        omega_o (str or list): Path for the sky projected inclination data file.
                                    It can be a path or a list where the first element is the mean value, the second the standard deviation, and the third the length of it. From it, a random distribution will be generated.
        P (float): Rotational period of the planet.
        M (float): Mass of the planet.

    Returns:
        (Exoplanet): An Exoplanet object.

    """
    # ------------------------------------------------------------------------------------------------------------------------------------------------
    ## Initialize ExoSpin
    print()
    print('Initializing ExoSpin ...' )
    print()
    print('-> ExoSpin Configuration')
    print()

    if isinstance(io,list):
        io_samp = np.random.normal(loc=io[0], scale=io[1], size=io[2])

    else:
        io_file = open(io, "r")
        io_samp = np.loadtxt(io_file, skiprows=1)

    if isinstance(radius,list):
        radius_samp = np.random.normal(loc=radius[0], scale=radius[1], size=radius[2])
    else:
        radius_file = open(radius, "r")
        radius_samp = np.loadtxt(radius_file, skiprows=1,usecols=(1,))

    if isinstance(vsini,list):
        vsini_samp = np.random.normal(loc=vsini[0], scale=vsini[1], size=vsini[2])

    else:
        vsini_file = open(vsini, "r")
        vsini_samp = np.loadtxt(vsini_file, skiprows=1,usecols=(2,))

    if isinstance(omega_o,list):
        omega_o_samp = np.random.normal(loc=omega_o[0], scale=omega_o[1], size=omega_o[2])

    else:
        omega_o_file = open(omega_o, "r")
        omega_o_samp = np.loadtxt(omega_o_file, skiprows=1)

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

