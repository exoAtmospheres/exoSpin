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
io_input = 'exoplanet/data/io.dat' #input('Where are your data for orbital inclination (° unit)? ')
io_file = open(io_input, "r")
io_samp = np.loadtxt(io_file, skiprows=1)
radius_input= 'exoplanet/data/rp_vrot.txt' #input('Where are your data for radius (R_Jup unit)? ')
radius_file = open(radius_input, "r")
radius_samp = np.loadtxt(radius_file, skiprows=1,usecols=(1,))
vsini_input= 'exoplanet/data/rp_vrot.txt' #input('Where are your data for rotational velocity (km/s unit)? ')
vsini_file = open(vsini_input, "r")
vsini_samp = np.loadtxt(vsini_file, skiprows=1,usecols=(2,))
omega_o_input = 'exoplanet/data/io.dat' #input('Where are your data for star inclination (° unit)? ')
omega_o_file = open(omega_o_input, "r")
omega_o_samp = np.loadtxt(omega_o_file, skiprows=1)
P = 2.1 #input('What is the period of the exoplanet? (h unit)? ')
# P = float(P)
M = 10 # input('What is the mass of the exoplanet? (M_Jup unit)? ')
# M = float(M)

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

    exoplanet.spin_axis_data()

    print('... obliquity computing ...')
    
    exoplanet.proj_obli_data()
    exoplanet.true_obli_data()



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

io_hist = exoplanet.plot_hist('Orbital inclination','green')
io_pdf = exoplanet.plot_pdf('Orbital inclination','green')

vsini_hist = exoplanet.plot_hist('Rotational velocity','blue')
vsini_pdf = exoplanet.plot_pdf('Rotational velocity','blue')

radius_hist = exoplanet.plot_hist('Radius','red')
radius_pdf = exoplanet.plot_pdf('Radius','red')

ip_hist = exoplanet.plot_hist('Spin axis', '#a7e6d7')
ip_pdf = exoplanet.plot_pdf('Spin axis - easy','#0d98ba')

pro_hist = exoplanet.plot_hist('Projected obliquity','purple')
pro_pdf = exoplanet.plot_pdf('Projected obliquity - easy','purple')

true_hist = exoplanet.plot_hist('True obliquity','orange')
true_pdf = exoplanet.plot_pdf('True obliquity - easy','orange')


# exoplanet/data/rp_vrot.txt