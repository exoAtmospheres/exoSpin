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
from pdf_functions import *

# --------------------------------------------------------------------------------------------------------------------------------------------------

class Exoplanet():

    def __init__(self,planet_name, io, radius, vsini, omega_o, period, mass):
        
        self.planet_name = planet_name

        # Checking if io is in rad or deg :
        if io[-1]<=np.pi:
            io=np.rad2deg(io)
        self.io = io

        # Checking if omega_o is in rad or deg :
        if omega_o[-1]<=np.pi:
            omega_o=np.rad2deg(omega_o)
        self.omega_o = omega_o

        # Setting units
        radius = radius * u.Rjup
        self.radius = radius

        vsini = vsini * u.km/u.s
        self.vsini = vsini

        mass = mass * u.Mjup
        self.mass = mass

        period = period * u.hr
        self.period = period

        velocity = 2*np.pi*radius.value/period 
        velocity = velocity.to(u.km/u.s)
        self.velocity = velocity

        # Setting limits
        self.v_lim = np.sqrt(c.G * mass*u.Mjup/(radius.value.max()*u.Rjup))
        self.P_lim = 2*np.pi*(radius.value**(3/2))/(np.sqrt(c.G*mass*u.Mjup))

        # Setting unknown parameters
        self.ip = None

        self.proj_obli = None

        self.true_obli = None


    def plot_hist(self,arg,color_graph):

        bins = 200

        # Which data to want to plot ?

        if arg=='Orbital inclination':
            fig = plt.figure()
            y, x, _ = plt.hist(self.io, bins=bins, density=True, color=color_graph,label='Distribution of $i_o$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Orbital inclination  \n $i_p$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Radius':
            fig = plt.figure()
            y, x, _ = plt.hist(self.radius, bins=bins, density=True, color=color_graph,label='Distribution of $R$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Radius  \n $R$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '$R_{Jup}$')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.show()
            return fig

        elif arg=='Rotational velocity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.vsini, bins=bins, density=True, color=color_graph,label='Distribution of $\\nu sin (i_p)$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Rotational velocity  \n $\\nu sin (i_p)$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '$km.s^{-1}$')
            plt.xlabel('Velocity ($km.s^{-1}$)})')
            plt.show()
            return fig

        elif arg=='Star inclination':
            fig = plt.figure()
            y, x, _ = plt.hist(self.omega_o, bins=bins, density=True, color=color_graph,label='Distribution of $\Omega_o$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Orbital inclination  \n $\Omega_o$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Spin axis':
            fig = plt.figure()
            y, x, _ = plt.hist(self.ip, bins=bins, density=True, color=color_graph,label='Distribution of $i_p$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution -  Spin axis  \n $i_p$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Projected obliquity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.proj_obli, bins=bins, density=True, color=color_graph,label='Distribution of $|i_p-i_o|$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Projected obliquity  \n $|i_p-i_o|$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='True obliquity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.proj_obli, bins=bins, density=True, color=color_graph,label='Distribution of $|i_p-i_o|$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Projected obliquity  \n $|i_p-i_o|$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        else:
            return None


    def plot_pdf(self,arg,color_graph):

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,1000)
            io_pdf = pdf(kde(self.io),angles)
            fig, ax = plt.subplot()
            ax.plot(angles,io_pdf,color=color_graph,label='PDF of $R$')
            ax.title('PDF - Orbital inclination  \n $i_p$ = '+ (str(angles[np.argmax(io_pdf)],2))+ '°')
            ax.xlabel('Degree (°)')
            ax.legend()
            plt.show()
            return fig, ax

        elif arg=='Radius':
            meters = np.linspace(0,4,1000)
            radius_pdf = pdf(kde(self.radius),meters)
            fig, ax = plt.subplot()
            ax.plot(meters,radius_pdf,color=color_graph,label='PDF of $R$')
            ax.xlabel('Meters ($R_{Jup}$)')
            ax.ylabel('PDF - Radius \n $R$ = '+ (str(round(angles[np.argmax(radius_pdf)],2))) + '°')
            ax.legend()
            plt.show()
            return fig, ax

        elif arg=='Rotational velocity':
            velocities = np.linspace(0,self.v_lim,1000)
            vsini_pdf = pdf(kde(self.vsini),velocities)
            fig, ax = plt.subplot()
            ax.plot(meters,radius_pdf,color=color_graph,label='PDF of $\\nu sin (i_p)$')
            plt.title('PDF - Rotational velocity  \n $\\nu sin (i_p)$ = '+ (str(round(velocities[np.argmax(vsini_pdf)],2)))+ '$km.s^{-1}$')
            plt.xlabel('Velocity ($km.s^{-1}$)})')
            ax.legend()
            plt.show()
            return fig, ax

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,1000)
            omega_o_pdf = pdf(kde(self.omega_o),angles)
            fig, ax = plt.subplot()
            ax.plot(angles,omega_o_pdf_pdf,color=color_graph,label='PDF of $\Omega_o$')
            ax.title('PDF - Start inclination  \n $\Omega_o$ = '+ (str(angles[np.argmax(omega_o_pdf)],2))+ '°')
            ax.xlabel('Degree (°)')
            ax.legend()
            plt.show()
            return fig, ax

        if arg=='Orbital inclination':
            user_input = input("Which computing method? (easy/complex)")
            while user_input!='easy' and user_input!='complex':
                print()
                print('You have to choose a method!')
                print()
                user_input = input("Which computing method? (easy/complex)")
            if user_input == 'easy':
                angles = np.linspace(0,180,1000)
                ip_pdf = pdf(kde(self.oip),angles)
                fig, ax = plt.subplot()
                ax.plot(angles,ip_pdf,color=color_graph,label='PDF of $i_p$')
                ax.title('PDF - Spin axis  \n $i_p$ = '+ (str(angles[np.argmax(ip_pdf)],2))+ '°')
                ax.xlabel('Degree (°)')
                ax.legend()
                plt.show()
                return fig, ax
            else:
                Lv = kde(self.velocity)
                Lu = kde(self.vsini)
                v_range = np.linspace(0,self.v_lim,self.ip.size)                                     
                n_points = 100                                                                          
                ip_pdf = ip_complex_pdf(Lv,Lu,v_range,n_points)
                ## Plot
                fig, ax = plt.subplot()
                ax.plot(angles,ip_pdf,color=color_graph,label='PDF of $i_p$')
                ax.title('PDF - Spin axis  \n $i_p$ = '+ (str(angles[np.argmax(ip_pdf)],2))+ '°')
                ax.xlabel('Degree (°)')
                ax.legend()
                plt.show()
                return fig, ax

        if arg=='Projected inclination':
            user_input = input("Which computing method? (easy/complex)")
            while user_input!='easy' and user_input!='complex':
                print()
                print('You have to choose a method!')
                print()
                user_input = input("Which computing method? (easy/complex)")
            if user_input == 'easy':
                angles = np.linspace(0,180,1000)
                pro_obli_pdf = pdf(kde(self.proj_obli),angles)
                fig, ax = plt.subplot()
                ax.plot(angles,pro_obli_pdf,color=color_graph,label='PDF of $|i_p-i_o|$')
                ax.title('PDF - Spin axis  \n $i_p$ = '+ (str(angles[np.argmax(ip_pdf)],2))+ '°')
                ax.xlabel('Degree (°)')
                ax.legend()
                plt.show()
                return fig, ax
            else:
                return None



    #def plot_all(self,arg):

    #def spin_axis_data(self):  

    #def proj_obli_data(self):

    #def true_obli_data(self):





        

    