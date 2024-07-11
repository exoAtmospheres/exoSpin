'''
ExoSpin run script - Exoplanet Class


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports

import sys, os, glob
import matplotlib.patches as patches
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
from plot_class import Plot

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

    ## Computing methods

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
    
        self.omega_o =np.random.choice(self.omega_o,self.ip.size)
        omega_p = np.random.uniform(0, 180, self.ip.size)
        self.lambda_ = self.omega_o-omega_p
    
        true_obli = np.arccos(np.cos(np.deg2rad(self.ip))*np.cos(np.deg2rad(self.io))+np.sin(np.deg2rad(self.ip))*np.sin(np.deg2rad(self.io))*np.cos(np.deg2rad(self.lambda_)))
        true_obli = np.rad2deg(true_obli)
        self.true_obli = true_obli

    ## Plot methods

    def hist(self,arg,color_graph):

        bins = 200
        if color_graph == None:
            color_graph = '#74D0F1'

        if arg=='Orbital inclination':
            y, x, _ = plt.hist(self.io, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Orbital inclination  \n $i_o$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
            xlabel='Degree (°)'
            plot = Plot('Histogram' , self.io , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Radius':
            y, x, _ = plt.hist(self.io, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Radius  \n $R$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' $R_{Jup}$'
            xlabel='Length ($R_{Jup}$)'
            plot = Plot('Histogram' , self.radius , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Rotational velocity':
            y, x, _ = plt.hist(self.io, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Rotation velocity  \n $\\nu sin (i_p)$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' $km.s^{-1}$'
            xlabel='Velocity ($km.s^{-1}$)'
            plot = Plot('Histogram' , self.vsini , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Sky projected inclination':
            y, x, _ = plt.hist(self.omega_o, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Sky projected inclination  \n $\Omega_o$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
            xlabel='Degree (°)'
            plot = Plot('Histogram' , self.omega_o , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Spin axis':
            y, x, _ = plt.hist(self.ip, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Spin axis  \n $\i_p$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
            xlabel='Degree (°)'
            plot = Plot('Histogram' , self.ip , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Projected obliquity':
            y, x, _ = plt.hist(self.proj_obli, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Projected obliquity  \n $|i_p-i_o|$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
            xlabel='Degree (°)'
            plot = Plot('Histogram' , self.proj_obli , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='True obliquity':
            y, x, _ = plt.hist(self.true_obli, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - True obliquity  \n $\Psi_{op}$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
            xlabel='Degree (°)'
            plot = Plot('Histogram' , self.omega_o , None , xlabel , None , color_graph , title)
            return plot

        else:
            print('Error - You need to set argument')
            return None

    def pdf(self,arg,color_graph):

        n_easy = 1000
        n_complex = 100
        if color_graph == None:
            color_graph = '#74D0F1'

        if arg=='Orbital inclination':
            angles = np.linspace(0,180,n_easy)
            io_pdf = pdf(kde(self.io),angles)
            title = 'PDF - Orbital inclination  \n $i_p$ = '+ str(round(angles[np.argmax(io_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , io_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Radius':
            meters = np.linspace(0,4,n_easy)
            radius_pdf = pdf(kde(self.radius.value),meters)
            title = 'PDF - Radius  \n $R$ = '+ str(round(meters[np.argmax(radius_pdf)],2))+ ' °'
            xlabel = 'Length ($R_{Jup}$)'
            ylabel = 'PDF'
            plot = Plot('PDF' , meters , radius_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Rotational velocity':
            velocities = np.linspace(0,self.v_lim.value,n_easy)
            vsini_pdf = pdf(kde(self.vsini.value),velocities)
            title = 'PDF - Rotational velocity  \n $\\nu sin (i_p)$ = '+ str(round(velocities[np.argmax(vsini_pdf)],2))+ ' °'
            xlabel = 'Velocity ($km.s^{-1}$)'
            ylabel = 'PDF'
            plot = Plot('PDF' , velocities , vsini_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Sky projected inclination':
            angles = np.linspace(0,180,n_easy)
            omega_o_pdf = pdf(kde(self.omega_o),angles)
            title = 'PDF - Sky projected inclination  \n $\Omega_o$ = '+ str(round(angles[np.argmax(omega_o_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , omega_o_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Spin axis - easy':
            n = self.ip.size
            ip_half = self.ip[:n//2]
            angles = np.linspace(0,90,n_easy//2)
            ip_pdf = pdf(kde(ip_half),angles)
            ip_pdf = np.concatenate((ip_pdf,ip_pdf[::-1]))
            new_angles = np.linspace(0,180,n_easy)
            title = 'PDF - Spin axis  \n $i_p$ = '+ str(round(angles[np.argmax(ip_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , new_angles , ip_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Spin axis - complex':
            angles = np.linspace(0,180,n_complex)
            v_range = np.linspace(0,self.v_lim.value,self.ip.size)                                                                                                             
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)
            ## Plot
            title = 'PDF - Spin axis  \n $i_p$ = '+ str(round(angles[np.argmax(ip_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , ip_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Projected obliquity - easy':
            angles = np.linspace(0,180,n_easy)
            pro_obli_pdf = pdf(kde(self.proj_obli),angles)
            title = 'PDF - Projected obliquity  \n $|i_p-i_o|$ = '+ str(round(angles[np.argmax(pro_obli_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , pro_obli_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='Projected obliquity - complex': 
            angles = np.linspace(0,180,n_complex)
            v_range = np.linspace(0,self.v_lim.value,self.ip.size)                                                                                                             
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)    
            pro_obli_pdf = proj_obli_complex_pdf(kde(np.deg2rad(self.io)),ip_pdf,n_complex)
            ## Plot
            title = 'PDF - Projected obliquity  \n $|i_p-i_o|$ = '+ str(round(angles[np.argmax(pro_obli_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , pro_obli_pdf , xlabel , ylabel , color_graph , title)
            return plot


        elif arg=='True obliquity - easy':
            angles = np.linspace(0,180,n_easy)
            true_obli_pdf = pdf(kde(self.true_obli),angles)
            title = 'PDF - True obliquity  \n $\Psi_{op}$ = '+ str(round(angles[np.argmax(true_obli_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , true_obli_pdf , xlabel , ylabel , color_graph , title)
            return plot

        elif arg=='True obliquity - complex': 

            io_kde = kde(np.deg2rad(self.io))
            lambda_kde = kde(np.deg2rad(self.lambda_)) 
            angles = np.linspace(0,180,n_complex)     
            v_range = np.linspace(0,self.v_lim.value,self.ip.size)                                                                                                             
            
            ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)  
            io_pdf = pdf(io_kda,angles)
            lambda_pdf = pdf(lambda_kde,angles)
            true_obli_pdf = true_obli_complex_pdf(io_pdf,ip_pdf,lambda_pdf,n_complex)
            ## Plot
            title = 'PDF - True obliquity  \n $\Psi_{op}$ = '+ str(round(angles[np.argmax(true_obli_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , true_obli_pdf , xlabel , ylabel , color_graph , title)
            return plot

    def plot_obli(self,arg):

        angles = np.linspace(0,180,1000)

        if arg == 'easy':

            io_plot = self.pdf('Orbital inclination','#000099')
            ip_plot = self.pdf('Spin axis - easy','#E73E01')
            proj_plot = self.pdf('Projected obliquity - easy','#A100A1')
            true_plot = self.pdf('True obliquity - easy','#87E990')

            fig_final = plt.figure()
            fig, axs = plt.subplots(2,2, figsize=(10,8))
            axs[0,0].plot(io_plot.x,io_plot.y,color=io_plot.color,label='$i_o$')

            axs[0,0].tick_params(axis='y', colors=io_plot.color)
            axs[0,0].set_title('Orbital inclination and companion spin axis')
            axs[0,0].set(ylabel='PDF')

            axs2 = axs[0,0].twinx()
            axs2.plot(ip_plot.x,ip_plot.y,color=ip_plot.color,label='$i_p$')
            axs2.spines['left'].set_color(io_plot.color)
            axs2.spines['right'].set_color(ip_plot.color)
            axs2.tick_params(axis='y', colors=ip_plot.color)
            axs2.set(ylabel='PDF')

            lines_1, labels_1 = axs[0, 0].get_legend_handles_labels()
            lines_2, labels_2 = axs2.get_legend_handles_labels()
            axs2.legend(lines_1 + lines_2, labels_1 + labels_2)


            #Diagram
            star = plt.Circle((0.35, 0.5), 0.07, color='#ffd319', ec='#bb6f1e', lw=2, label=self.planet_name + '\'s star')
            orbit = patches.Ellipse((0.5,0.5), 0.70,0.45, ec='black', lw=2, linestyle = 'dotted', fill=False )
            exoplanet = plt.Circle((0.85,0.5), 0.03, color = '#133984', ec='#c3dbff', lw=1, label = self.planet_name)

            ip_max = np.deg2rad(angles[np.argmax(io_plot.x)])

            # Define the spin axis line through the planet
            spin_length = 0.05
            x0, y0 = 0.85, 0.5
            x1 = x0 + spin_length * np.cos(ip_max)
            y1 = y0 + spin_length * np.sin(ip_max)
            x2 = x0 - spin_length * np.cos(ip_max)
            y2 = y0 - spin_length * np.sin(ip_max)


            axs[0,1].set_title('Diagram of ' + self.planet_name + '\'s system')
            axs[0,1].add_patch(star)
            axs[0,1].add_patch(orbit)
            axs[0,1].add_patch(exoplanet)
            axs[0,1].plot([x2, x1], [y2, y1], color=ip_plot.color, label='Spin Axis')
            axs[0,1].set_xlim(0, 1)
            axs[0,1].set_ylim(0, 1)
            axs[0,1].set_aspect('equal')
            axs[0,1].set_axis_off()
            axs[0,1].legend()


            axs[1,0].plot(proj_plot.x,proj_plot.y,color=proj_plot.color,label='$|i_p-i_o|$')
            axs[1,0].set_title('Projected obliquity')
            axs[1,0].set(ylabel='PDF')
            axs[1,0].set(xlabel='Degree (°)')
            axs[1,0].legend()

            axs[1,1].plot(true_plot.x,true_plot.y,color=true_plot.color,label='$\Psi_{op}$')
            axs[1,1].set_title('True obliquity')
            axs[1,1].set(xlabel='Degree (°)')
            axs[1,1].legend()
            plt.show()

        elif arg == 'complex':

            io_plot = self.pdf('Orbital inclination','#dad45e')
            ip_plot = self.pdf('Spin axis - complex', 'green')
            proj_plot = self.pdf('Projected obliquity - complex','#dad45e')
            true_plot = self.pdf('True obliquity - easy','#d04648')



            fig_final = plt.figure()
            fig, axs = plt.subplots(2,2,figsize=(15,5))
            axs[0,0].plot(io_plot.x,io_plot.y,color=io_plot.color,label='PDF of $i_o$')
            axs[0,0].plot(ip_plot.x,ip_plot.y,color=ip_plot.color,label='PDF of $i_p$')
            axs[0,0].set_title('Orbital inclination and spin axis')
            axs[0,0].set(ylabel='PDF')
            plt.legend()

            axs[1,0].plot(proj_plot.x,proj_plot.y,color=proj_plot.color,label='PDF of $|i_p-i_o|$')
            axs[1,0].set_title('Projected obliquity')
            axs[1,0].set(ylabel='PDF')
            axs[1,0].set(xlabel='Degree (°)')
            plt.legend()

            axs[1,1].plot(true_plot.x,true_plot.y,color=true_plot.color,label='PDF of $\Psi_{op}$')
            axs[1,1].set_title('True obliquity')
            axs[1,1].set(xlabel='Degree (°)')
            plt.legend()

            plt.figure(fig_final.number)
            plt.show()

        else : 
            return None

    ## Set methods

    def set_data(self,arg,input_):
        if arg=='Orbital inclination':
            self.io = input_
        elif arg=='Spin axis':
            self.ip = input_
        elif arg=='Rotational velocity':
            self.io = input_
        elif arg=='Radius':
            self.radius = input_
        elif arg=='Sky projected inclination':
            self.omega_o = input_
        elif arg=='Mass':
            self.mass = input_
        elif arg=='Period':
            self.period = input_
        elif arg=='Planet name':
            self.planet_name = input_
        else:
            return None

    # def get_data(self, arg):

    #     if arg=='Orbital inclination':
    #         return self.io

    #     elif arg=='Radius':
    #         return self.radius

    #     elif arg=='Rotational velocity':
    #         return self.vsini

    #     elif arg=='Star inclination':
    #         return self.omega_o
        
    #     elif arg=='Spin axis':
    #         return self.ip
        
    #     elif arg=='Projected obliquity':
    #         return self.proj_obli
        
    #     elif arg=='True obliquity':
    #         return self.true_obli

    # def get_pdf(self, arg):
        
    #     n_easy = 1000
    #     n_complex = 100

    #     if arg=='Orbital inclination':
    #         angles = np.linspace(0,180,n_easy)
    #         io_pdf = pdf(kde(self.io),angles)
    #         return io_pdf

    #     elif arg=='Radius':
    #         meters = np.linspace(0,4,n_easy)
    #         radius_pdf = pdf(kde(self.radius),meters)
    #         return radius_pdf

    #     elif arg=='Rotational velocity':
    #         velocities = np.linspace(0,self.v_lim,n_easy)
    #         vsini_pdf = pdf(kde(self.vsini),velocities)
    #         return vsini_pdf

    #     if arg=='Orbital inclination':
    #         angles = np.linspace(0,180,n_easy)
    #         omega_o_pdf = pdf(kde(self.omega_o),angles)
    #         return omega_o_pdf

    #     if arg=='Spin axis':
    #         user_input = input('Which method of computing do you want? (easy/complex) ')
    #         while user_input!='easy' and user_input!='complex':
    #             print()
    #             print('You have to choose a method!')
    #             print()
    #             user_input = input('Which method of computing do you want? (easy/complex) ')
    #         if user_input == 'easy':
    #             angles = np.linspace(0,180,n_easy)
    #             ip_pdf = pdf(kde(self.ip),angles)
    #             return ip_pdf
    #         else:
    #             v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
    #             ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)
    #             return ip_pdf

    #     if arg=='Projected obliquity':
    #         user_input = input("Which computing method? (easy/complex)")
    #         while user_input!='easy' and user_input!='complex':
    #             print()
    #             print('You have to choose a method!')
    #             print()
    #             user_input = input("Which computing method? (easy/complex)")
    #         if user_input == 'easy':
    #             angles = np.linspace(0,180,n_easy)
    #             pro_obli_pdf = pdf(kde(self.proj_obli),angles)
    #             return fig, ax
    #         else:                  
    #             v_range = np.linspace(0,self.v_lim,self.ip.size)                                                                                                             
    #             ip_pdf = ip_complex_pdf(kde(self.velocity),kde(self.vsini),v_range,n_complex)    
    #             pro_obli_pdf = proj_obli_complex_pdf(kde(np.deg2rad(self.io)),ip,n_complex)
    #             return pro_obli_pdf


    #     if arg=='True obliquity':
    #         user_input = input("Which computing method? (easy/complex)")
    #         while user_input!='easy' and user_input!='complex':
    #             print()
    #             print('You have to choose a method!')
    #             print()
    #             user_input = input("Which computing method? (easy/complex)")
    #         if user_input == 'easy':
    #             angles = np.linspace(0,180,n_easy)
    #             true_obli_pdf = pdf(kde(self.true_obli),angles)
    #             return fig, ax
    #         else:                  
    #             print('Working on')











        

    