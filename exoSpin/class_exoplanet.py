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
        self.radius = radius * u.Rjup

        self.vsini = vsini * u.km/u.s

        self.mass = mass * u.Mjup

        self.period = period * u.hr

        velocity = 2*np.pi*self.radius/self.period 
        velocity=velocity.to(u.km/u.s)
        self.velocity = velocity

        # Setting limits
        v_lim = np.sqrt(c.G *self.mass/(self.radius.max()))
        v_lim = v_lim.to(u.km/u.s)
        self.v_lim = v_lim
        
        P_lim = 2*np.pi*(self.radius**(3/2))/(np.sqrt(c.G*self.mass))
        P_lim = P_lim.to(u.hr)
        self.P_lim = P_lim

        # Setting unknown parameters
        self.ip = None

        self.proj_obli = None

        self.true_obli = None

    def spin_axis_data(self):
        P_sample = (self.period > self.P_lim)
        # Set vel and visini with P > P_limit condition and v < v_limit
        self.vsini=self.vsini[P_sample]
        self.velocity=self.velocity[P_sample]
        # Velocity limitation  due to centrifugal force and gravitationnal force
        v_sample= (self.vsini < self.velocity) 

        # Generate ip histogram
        sin_ip = self.vsini/self.velocity
        ip = np.arcsin(sin_ip[v_sample])
        ip = np.concatenate((ip.value, np.pi-ip.value))
        self.ip=np.rad2deg(ip)

    def proj_obli_data(self):

        if len(self.io) > self.ip.size:
            self.io =np.random.choice(self.io,self.ip.size)
        else:
            self.ip = np.random.choice(self.ip, self.io.size)

        self.proj_obli = np.abs(self.io-self.ip)

    def true_obli_data(self):
        omega_p = np.random.uniform(0, 180, self.ip.size)
        psi_op = np.arccos(np.cos(np.deg2rad(self.io))*np.cos(np.deg2rad(self.ip))+np.sin(np.deg2rad(self.io))*np.sin(np.deg2rad(self.ip))*np.cos(np.deg2rad(omega_p)))
        self.true_obli = np.rad2deg(psi_op)

    def plot_hist(self,arg,color_graph):

        bins = 200

        # Which data to want to plot ?

        if arg=='Orbital inclination':
            fig = plt.figure()
            y, x, _ = plt.hist(self.io, bins=bins, density=True, color=color_graph,label='Distribution of $i_o$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Orbital inclination  \n $i_o$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Radius':
            fig = plt.figure()
            y, x, _ = plt.hist(self.radius, bins=bins, density=True, color=color_graph,label='Distribution of $R$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Radius  \n $R$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' $R_{Jup}$')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.show()
            return fig

        elif arg=='Rotational velocity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.vsini, bins=bins, density=True, color=color_graph,label='Distribution of $\\nu sin (i_p)$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Rotational velocity  \n $\\nu sin (i_p)$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' $km.s^{-1}$')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.show()
            return fig

        elif arg=='Star inclination':
            fig = plt.figure()
            y, x, _ = plt.hist(self.omega_o, bins=bins, density=True, color=color_graph,label='Distribution of $\Omega_o$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Orbital inclination  \n $\Omega_o$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Spin axis':
            fig = plt.figure()
            y, x, _ = plt.hist(self.ip, bins=bins, density=True, color=color_graph,label='Distribution of $i_p$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution -  Spin axis  \n $i_p$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='Projected obliquity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.proj_obli, bins=bins, density=True, color=color_graph,label='Distribution of $|i_p-i_o|$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - Projected obliquity  \n $|i_p-i_o|$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        elif arg=='True obliquity':
            fig = plt.figure()
            y, x, _ = plt.hist(self.true_obli, bins=bins, density=True, color=color_graph,label='Distribution of $\Psi_{op}$')
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            plt.title('Distribution - True obliquity  \n $\Psi_{op}$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °')
            plt.xlabel('Degree (°)')
            plt.show()
            return fig

        else:
            return None

    def plot_pdf(self,arg,color_graph):

        n_easy = 1000
        n_complex = 100

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,n_easy)
            io_pdf = pdf(kde(self.io),angles)
            fig = plt.figure()
            plt.plot(angles,io_pdf,color=color_graph,label='PDF of $i_o$')
            plt.title('PDF - Orbital inclination  \n $i_p$ = '+ str(round(angles[np.argmax(io_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            return fig

        elif arg=='Radius':
            meters = np.linspace(0,4,n_easy)
            radius_pdf = pdf(kde(self.radius.value),meters)
            fig = plt.figure()
            plt.plot(meters,radius_pdf,color=color_graph,label='PDF of $R$')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.title('PDF - Radius \n $R$ = '+ str(round(meters[np.argmax(radius_pdf)],2)) + ' $R_{Jup}$')
            plt.legend()
            return fig

        elif arg=='Rotational velocity':
            velocities = np.linspace(0,self.v_lim.value,n_easy)
            vsini_pdf = pdf(kde(self.vsini.value),velocities)
            fig = plt.figure()
            plt.plot(velocities,vsini_pdf,color=color_graph,label='PDF of $\\nu sin (i_p)$')
            plt.title('PDF - Rotational velocity  \n $\\nu sin (i_p)$ = '+ str(round(velocities[np.argmax(vsini_pdf)],2))+ ' $km.s^{-1}$')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.legend()
            plt.show()
            return fig

        elif arg=='Star inclination':
            angles = np.linspace(0,180,n_easy)
            omega_o_pdf = pdf(kde(self.omega_o),angles)
            fig = plt.figure()
            plt.plot(angles,omega_o_pdf_pdf,color=color_graph,label='PDF of $\Omega_o$')
            plt.title('PDF - Start inclination  \n $\Omega_o$ = '+ (str(angles[np.argmax(omega_o_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig

        elif arg=='Spin axis - easy':
            n = self.ip.size
            ip_half = self.ip[:n//2]
            angles = np.linspace(0,90,n_easy//2)
            ip_pdf = pdf(kde(ip_half),angles)
            ip_pdf = np.concatenate((ip_pdf,ip_pdf[::-1]))
            new_angles = np.linspace(0,180,n_easy)
            fig = plt.figure()
            plt.plot(new_angles,ip_pdf,color=color_graph,label='PDF of $i_p$')
            plt.title('PDF - Spin axis  \n $i_p$ = '+ str(round(angles[np.argmax(ip_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig
        elif arg=='Spin axis - complex':
            v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)
            ## Plot
            fig = plt.figure()
            plt.plot(angles,ip_pdf,color=color_graph,label='PDF of $i_p$')
            plt.title('PDF - Spin axis  \n $i_p$ = '+ str(round(angles[np.argmax(ip_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig

        elif arg=='Projected obliquity - easy':
            angles = np.linspace(0,180,n_easy)
            pro_obli_pdf = pdf(kde(self.proj_obli),angles)
            fig = plt.figure()
            plt.plot(angles,pro_obli_pdf,color=color_graph,label='PDF of $|i_p-i_o|$')
            plt.title('PDF - Projected obliquity  \n $|i_p-i_o|$ = '+ str(round(angles[np.argmax(pro_obli_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig
        elif arg=='Projected obliquity - complex': 
            v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)    
            pro_obli_pdf = proj_obli_complex_pdf(kde(np.deg2rad(self.io)),ip_pdf,v_range,n_complex)
            ## Plot
            fig = plt.figure()
            plt.plot(angles,ip_pdf,color=color_graph,label='PDF of $|i_p-i_o|$')
            plt.title('PDF - Projected obliquity  \n $|i_p-i_o|$ = '+ str(round(angles[np.argmax(pro_obli_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig

        elif arg=='True obliquity - easy':
            angles = np.linspace(0,180,n_easy)
            true_obli_pdf = pdf(kde(self.true_obli),angles)
            fig = plt.figure()
            plt.plot(angles,true_obli_pdf,color=color_graph,label='PDF of $\Psi_{op}$')
            plt.title('PDF - True obliquity  \n $\Psi_{op}$ = '+ str(round(angles[np.argmax(true_obli_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig

        elif arg=='True obliquity - complex':       
            v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)    
            pro_obli_pdf = proj_obli_complex_pdf(kde(np.deg2rad(self.io)),ip_pdf,v_range,n_complex)
            ## Plot
            fig = plt.figure()
            plt.plot(angles,ip_pdf,color=color_graph,label='PDF of $\Psi_{op}$')
            plt.title('PDF - True obliquity  \n $\Psi_{op}$ = '+ str(round(angles[np.argmax(true_obli_pdf)],2))+ ' °')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            return fig

    def get_pdf(self, arg):
        
        n_easy = 1000
        n_complex = 100

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,n_easy)
            io_pdf = pdf(kde(self.io),angles)
            return io_pdf

        elif arg=='Radius':
            meters = np.linspace(0,4,n_easy)
            radius_pdf = pdf(kde(self.radius),meters)
            return radius_pdf

        elif arg=='Rotational velocity':
            velocities = np.linspace(0,self.v_lim,n_easy)
            vsini_pdf = pdf(kde(self.vsini),velocities)
            return vsini_pdf

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,n_easy)
            omega_o_pdf = pdf(kde(self.omega_o),angles)
            return omega_o_pdf

        if arg=='Spin axis':
            user_input = input('Which method of computing do you want? (easy/complex) ')
            while user_input!='easy' and user_input!='complex':
                print()
                print('You have to choose a method!')
                print()
                user_input = input('Which method of computing do you want? (easy/complex) ')
            if user_input == 'easy':
                angles = np.linspace(0,180,n_easy)
                ip_pdf = pdf(kde(self.ip),angles)
                return ip_pdf
            else:
                v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
                ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)
                return ip_pdf

        if arg=='Projected obliquity':
            user_input = input("Which computing method? (easy/complex)")
            while user_input!='easy' and user_input!='complex':
                print()
                print('You have to choose a method!')
                print()
                user_input = input("Which computing method? (easy/complex)")
            if user_input == 'easy':
                angles = np.linspace(0,180,n_easy)
                pro_obli_pdf = pdf(kde(self.proj_obli),angles)
                return fig, ax
            else:                  
                v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
                ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)    
                pro_obli_pdf = proj_obli_complex_pdf(kde(np.deg2rad(self.io)),ip,n_complex)
                return pro_obli_pdf


        if arg=='True obliquity':
            user_input = input("Which computing method? (easy/complex)")
            while user_input!='easy' and user_input!='complex':
                print()
                print('You have to choose a method!')
                print()
                user_input = input("Which computing method? (easy/complex)")
            if user_input == 'easy':
                angles = np.linspace(0,180,n_easy)
                true_obli_pdf = pdf(kde(self.true_obli),angles)
                return fig, ax
            else:                  
                print('Working on')

    #def plot_all(self,arg):










        

    