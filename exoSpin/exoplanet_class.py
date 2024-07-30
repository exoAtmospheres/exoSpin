'''
ExoSpin - Exoplanet Class


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

from scipy.signal import savgol_filter
from scipy.stats import gaussian_kde
from scipy.stats import uniform

from exoSpin.pdf_functions import *
from exoSpin.plot_class import Plot

# --------------------------------------------------------------------------------------------------------------------------------------------------

class Exoplanet():

    def __init__(self,planet_name, io, radius, vsini, omega_o, period, mass):
        """
        Initialize every parameters for the Exoplanet class object.

        Args:
            planet_name (String) : Name of the exoplanet.
            io (numpy.ndarray) : 1D array representing data of the orbital inclination system.
            radius (astropy.units.quantity.Quantity) : 1D array representing data of the exoplanet's radius.
            vsini (astropy.units.quantity.Quantity) : 1D array representing data of the exoplanet's rotational velocity.
            omega_o (numpy.ndarray) : 1D array representing data of the exoplanet's sky projected inclination.
            period (astropy.units.quantity.Quantity) : a quantity representing the exoplanet's period.
            mass (astropy.units.quantity.Quantity) : a quantity representing the exoplanet's mass.

        Attributes:
            planet_name (String) : Name of the exoplanet.
            io (numpy.ndarray) : 1D array representing data of the orbital inclination system.
            radius (astropy.units.quantity.Quantity) : 1D array representing data of the exoplanet's radius.
            vsini (astropy.units.quantity.Quantity) : 1D array representing data of the exoplanet's rotational velocity.
            omega_o (numpy.ndarray) : 1D array representing data of the exoplanet's sky projected inclination.
            period (astropy.units.quantity.Quantity) : a quantity representing the exoplanet's period.
            mass (astropy.units.quantity.Quantity) : a quantity representing the exoplanet's mass.
            velocity (numpy.ndarray) : 1D array representing data of the exoplanet's velocity computed with radius and period data.
            v_lim (astropy.units.quantity.Quantity) : a quantity representing the velocity break limit of the exoplanet.
            P_lim (astropy.units.quantity.Quantity) : a quantity representing the period break limit of the exoplanet.
            lambda_ (numpy.ndarray) : 1D array representing data of the sky projected obliquity.
            ip (numpy.ndarray) : 1D array representing data of the exoplanet's spin axis.
            proj_obli (numpy.ndarray) : 1D array representing data of the exoplanet's projected obliquity.
            true_obli (numpy.ndarray) : 1D array representing data of the exoplanet's true obliquity.
            ip_pdf_saved (Plot) : a Plot that contains PDF data for spin axis

        """

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
        self.lambda_ = None
        
        self.ip = None

        self.proj_obli = None

        self.true_obli = None

        self.ip_plot_saved = None

    ## Computing methods

    def spin_axis_data(self):
        """
        Compute the spin axis data of the exoplanet
        """

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
        """
        Compute the projected obliquity data of the exoplanet
        """

        if len(self.io) > self.ip.size:
            self.io =np.random.choice(self.io,self.ip.size)
        else:
            self.ip = np.random.choice(self.ip, self.io.size)

        self.proj_obli = np.abs(self.io-self.ip)

    def true_obli_data(self):
        """
        Compute the true obliquity data of the exoplanet
        """
    
        self.omega_o =np.random.choice(self.omega_o,self.ip.size)
        omega_p = np.random.uniform(0, 180, self.ip.size)
        self.lambda_ = omega_p 
    
        true_obli = np.arccos(np.cos(np.deg2rad(self.ip))*np.cos(np.deg2rad(self.io))+np.sin(np.deg2rad(self.ip))*np.sin(np.deg2rad(self.io))*np.cos(np.deg2rad(self.lambda_)))
        true_obli = np.rad2deg(true_obli)
        self.true_obli = true_obli

    ## Plot methods

    def hist(self,arg,color_graph = 'blue'):
        """
        Compute and save histogram plot of a exoplanet parameter.

        Args:
            arg (String): a String that represents which histogram to plot.
                        {
                        'Orbital inclination',
                        'Radius',
                        'Rotational velocity',
                        'Sky projected inclination',
                        'Sky projected spin-orbit',
                        'Spin axis',
                        'Project obliquity',
                        'True obliquity'
                        }
            color_graph (String): a string to set the graph color.

        Returns:
            Plot : a Plot object that has every important parameters of the histogram 

        Raises:
            If the arg is not in the excepted values.    
        """

        if arg not in {
            'Orbital inclination',
            'Radius',
            'Rotational velocity',
            'Sky projected inclination',
            'Sky projected spin-orbit',
            'Spin axis',
            'Project obliquity',
            'True obliquity'
        }:
            raise ValueError("The arg value is not in the expected values")
                        
        bins = 200

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
            y, x, _ = plt.hist(self.radius, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Radius  \n $R$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' $R_{Jup}$'
            xlabel='Length ($R_{Jup}$)'
            plot = Plot('Histogram' , self.radius , None , xlabel , None , color_graph , title)
            return plot

        elif arg=='Rotational velocity':
            y, x, _ = plt.hist(self.vsini, bins=bins, density=True)
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

        elif arg=='Sky projected obliquity':
            y, x, _ = plt.hist(self.lambda_, bins=bins, density=True)
            plt.close()
            x_max =  x[np.where(y == y.max())][0]
            x_err    =  np.std(x_max)
            title='Distribution - Sky projected obliquity  \n $\lambda$ = '+ (str(round(x_max,2)))+ '$\pm$'+ str(x_err) + ' °'
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

    def pdf(self,arg,color_graph = 'blue'):
        """
        Compute and save the PDF plot of a exoplanet parameter.

        Args:
            arg (String): a String that represents which histogram to plot.
                        {
                        'Orbital inclination',
                        'Radius',
                        'Rotational velocity',
                        'Sky projected inclination',
                        'Sky projected spin-orbit',
                        'Spin axis - easy',
                        'Spin axis - complex',
                        'Projected obliquity - easy',
                        'Projected obliquity - complex',
                        'True obliquity - easy',
                        'True obliquity - complex'
                        }
            color_graph (String): a string to set the graph color.

        Returns:
            Plot : a Plot object that has every important parameters of the histogram 

        Raises:
            If the arg is not in the excepted values.    
        """

        if arg not in {
            'Orbital inclination',
            'Radius',
            'Rotational velocity',
            'Sky projected inclination',
            'Sky projected spin-orbit',
            'Spin axis - easy',
            'Spin axis - complex',
            'Projected obliquity - easy',
            'Projected obliquity - complex',
            'True obliquity - easy',
            'True obliquity - complex',
        }:
            raise ValueError("The arg value is not in the expected values")
               
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
    
        elif arg=='Sky projected obliquity':
            angles = np.linspace(0,180,n_easy)
            lambda_pdf = pdf(kde(self.lambda_),angles)
            title = 'PDF - Sky projected obliquity  \n $\lambda$ = '+ str(round(angles[np.argmax(lambda_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , lambda_pdf , xlabel , ylabel , color_graph , title)
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
            self.ip_pdf_saved = plot
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
            angles = self.ip_pdf_saved.x                                                                                                             
            ip_pdf = self.ip_pdf_saved.y    
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

            io_kde = kde(self.io)
            lambda_kde = kde(self.lambda_) 
            angles = np.linspace(0,180,2*n_easy)     
            v_range = np.linspace(0,self.v_lim.value,self.ip.size)                                                                                                             
            
            ip_pdf = self.ip_pdf_saved.y   
            io_pdf = pdf(io_kde,angles)
            lambda_pdf = pdf(lambda_kde,angles)

            angles_new = np.linspace(0,180,n_complex)

            true_obli_pdf = true_obli_complex_pdf(io_pdf,ip_pdf,lambda_pdf,2*n_easy)
            ## Plot
            title = 'PDF - True obliquity  \n $\Psi_{op}$ = '+ str(round(angles[np.argmax(true_obli_pdf)],2))+ ' °'
            xlabel = 'Degree (°)'
            ylabel = 'PDF'
            plot = Plot('PDF' , angles , true_obli_pdf , xlabel , ylabel , color_graph , title)
            return plot

    def plot_obli(self,arg):
        """
        Plot the obliquity of the exoplanet.

        Args:
            arg (String): a string that describe which histogram to plot.
                        {
                        'easy',
                        'complex'
                        }
            color_graph (String): a string to set the graph color.

        Raises:
            If  arg is not in the excepted values
        """

        if arg not in {
            'easy',
            'complex'
        }:
            raise ValueError("The arg value is not in the expected values")

        if arg == 'easy':

            io_plot = self.pdf('Orbital inclination','#000099')
            ip_plot = self.pdf('Spin axis - easy','#E73E01')
            proj_plot = self.pdf('Projected obliquity - easy','#A100A1')
            true_plot = self.pdf('True obliquity - easy','#87E990')

            print('Plot - Obliquity of ' + self.planet_name)

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
            exoplanet = plt.Circle((0.85,0.5), 0.03, color = '#c3dbff', ec='#133984', lw=1, label = self.planet_name)

            ip_max = np.deg2rad(ip_plot.y[np.argmax(ip_plot.x)])

            # Define the spin axis line

            spin_line = 0.05
            x0, y0 = 0.85, 0.5
            x1 = x0 + spin_line * np.cos(ip_max)
            y1 = y0 + spin_line * np.sin(ip_max)
            x2 = x0 - spin_line * np.cos(ip_max)
            y2 = y0 - spin_line * np.sin(ip_max)

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

            io_plot = self.pdf('Orbital inclination','#000099')
            ip_plot = self.pdf('Spin axis - complex','#E73E01')
            proj_plot = self.pdf('Projected obliquity - complex','#A100A1')
            true_plot = self.pdf('True obliquity - complex','#87E990')

            print('Plot - Obliquity of ' + self.planet_name)

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
            exoplanet = plt.Circle((0.85,0.5), 0.03, color = '#c3dbff', ec='#133984', lw=1, label = self.planet_name)

            ip_max = np.deg2rad(ip_plot.y[np.argmax(ip_plot.x)])

            # Define the spin axis line

            spin_line = 0.05
            x0, y0 = 0.85, 0.5
            x1 = x0 + spin_line * np.cos(ip_max)
            y1 = y0 + spin_line * np.sin(ip_max)
            x2 = x0 - spin_line * np.cos(ip_max)
            y2 = y0 - spin_line * np.sin(ip_max)


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

    ## Set method

    def set_data(self,arg,data):
        """
        Set new data for an exoplanet parameter. 

        Args:
            arg (String): a string that describe which histogram to plot.
                        {
                        'Orbital inclination',
                        'Radius',
                        'Rotational velocity',
                        'Sky projected inclination',
                        'Sky projected spin-orbit',
                        'Spin axis',
                        'Project obliquity',
                        'True obliquity'
                        }
            data: a input that contains the data to be set

        Raises:
            If the arg is not in the excepted values and if data and self attributes don't have the same type
        """

        if arg not in {
            'Orbital inclination',
            'Radius',
            'Rotational velocity',
            'Sky projected inclination',
            'Sky projected obliquity',
            'Mass',
            'Period'
            'Planet name'

        }:
            raise ValueError("The arg value is not in the expected values")

        if arg=='Orbital inclination':
            if type(self.io)==type(data):
                self.io = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Rotational velocity':
            if type(self.vsini)==type(data):
                self.vsini = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Radius':
            if type(self.radius)==type(data):
                self.radius = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Sky projected inclination':
            if type(self.omega_o)==type(data):
                self.omega_o = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Mass':
            if type(self.mass)==type(data):
                self.mass = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Period':
            if type(self.period)==type(data):
                self.period = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 
        elif arg=='Planet name':
            if type(self.plna)==type(data):
                self.planet_name = data
            else : 
               raise ValueError("Data input and self attribute don't have the same type") 