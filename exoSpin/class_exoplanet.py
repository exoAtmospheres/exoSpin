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

# --------------------------------------------------------------------------------------------------------------------------------------------------

class Exoplanet():

    def __init__(self,planet_name, io, radius, vsini, omega_o, period, mass):
        
        self.name = planet_name

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

        self.velocity = 2*np.pi*radius.value/period

        # Setting limits
        self.v_lim = np.sqrt(c.G * mass*u.Mjup/(radius.value.max()*u.Rjup))
        self.P_lim = 2*np.pi*(radius.value**(3/2))/(np.sqrt(c.G*mass*u.Mjup))

        # Setting unknows parameters
        self.ip = None
        self.ip_pdf = None

        self.proj_obli = None
        self.proj_obli_pdf = None

        self.true_obli = None
        self.true_obli_pdf = None

    # def order_data(self):

    def get_pdf(self,arg):
        if arg == 'orbital_inclination': 
            return gaussian_kde(self.io,bw_method='scott')

        elif arg == 'rotation_velocity': 
            return gaussian_kde(self.vsini,bw_method='scott')
        
        elif arg == 'radius': 
            return gaussian_kde(self.radius,bw_method='scott')

        elif arg == 'star_inclination': 
            return gaussian_kde(self.omega_o,bw_method='scott')

        elif arg == 'spin_axis_easy': 
            return gaussian_kde(self.ip,bw_method='scott')

        elif arg == 'spin_axis_complex':
            print('Working on..')

        elif arg == 'projected_obliquity_easy': 
            return gaussian_kde(self.ip,bw_method='scott')
        
        elif arg == 'projected_obliquity_complex':
            print('Working on..')
        
        elif arg == 'true_obliquity_easy':
            return gaussian_kde(self.true_obli,bw_method='scott')
        
        elif arg == 'true_obliquity_complex':
            print('Working on..')

    def get_orbital_inclination(self,arg):
        bins = 200
        # Compute io pdf
        angles = np.linspace(0,180,1000)
        io_kde = self.get_pdf('orbital_inclination')
        io_pdf = io_kde(angles)

        if arg=='pdf':
            plt.figure()
            plt.plot(angles,io_pdf,color='green',label='PDF of $i_o$')
            plt.xlabel('Degree (°)')
            plt.ylabel('PDF - Orbital inclination \n $i_p$ = '+ (str(round(angles[np.argmax(io_pdf)],2))) + '°')
            plt.legend()
            plt.show()
        
        elif arg=='distribution':
            plt.figure()
            y, x, _ = plt.hist(self.io, bins=bins, density=True, color='#50C878',label='Distribution of $i_o$')
            ip_max =  x[np.where(y == y.max())][0]
            ip_err    =  np.std(ip_max)
            plt.title('Distribution - Orbital inclination  \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
            plt.xlabel('Degree (°)')
            plt.show()
        else : 
            plt.figure()
            y, x, _ = plt.hist(self.io, bins=bins, density=True, color='#50C878',label='Distribution of $i_po')
            ip_max =  x[np.where(y == y.max())][0]
            ip_err    =  np.std(ip_max)
            plt.title('Orbital inclination  \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
            plt.xlabel('Degree (°)')
            plt.plot(angles,io_pdf,color='green',label='PDF of $i_o$')
            plt.xlabel('Degree (°)')
            plt.legend()
            plt.show()
            
    def get_radius(self,arg):
        
        # Compute radius pdf
        meters = np.linspace(0,4,self.radius.size)
        radius_kde = self.get_pdf('radius')
        radius_pdf = radius_kde(meters)

        if arg=='pdf':
            plt.figure()
            plt.plot(meters,radius_pdf,color='#666666',label='PDF of $R$')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.ylabel('PDF - Radius \n $R$ = '+ (str(round(angles[np.argmax(radius_pdf)],2))) + '°')
            plt.legend()
            plt.show()
        
        elif arg=='distribution':
            plt.figure()
            y, x, _ = plt.hist(self.radius.value, bins=bins, density=True, color='#50C878',label='Distribution of $R$')
            radius_max =  x[np.where(y == y.max())][0]
            radius_err    =  np.std(radius_max)
            plt.title('Distribution - PDF - Radius \n $R$ = '+ (str(round(radius_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.show()
        else : 
            plt.figure()
            y, x, _ = plt.hist(self.radius.value, bins=bins, density=True, color='#50C878',label='Distribution of $R$')
            radius_max =  x[np.where(y == y.max())][0]
            radius_err    =  np.std(radius_max)
            plt.title('Distribution - PDF - Radius \n $R$ = '+ (str(round(radius_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.plot(meters,radius_pdf,color='#666666',label='PDF of $R$')
            plt.xlabel('Meters ($R_{Jup}$)')
            plt.legend()
            plt.show()
 
    def get_rot_velocity(self,arg):

        # Period limitation due to centrifugal force and gravitationnal force
        P_sample = (self.P > self.P_lim)

        # Set vel and visini with P > P_limit condition and v < v_limit
        self.vsini.value=self.vsini.value[P_sample]

        # Compute vsini pdf
        velocities = np.linspace(0,self.v_lim,1000)
        vsini_kde = self.get_pdf('rotational_velocity')
        vsini_pdf_pdf = io_kde(velocities)

        if arg=='pdf':
            plt.figure()
            plt.plot(velocities,vsini_pdf,color='#666666',label='PDF of $R$')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.ylabel('PDF - Rotational velocity \n $\\nu \\sin i_p$ = '+ (str(round(angles[np.argmax(vsini_pdf)],2))) + '°')
            plt.legend()
            plt.show()
        
        elif arg=='distribution':
            plt.figure()
            y, x, _ = plt.hist(self.vsini.value, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
            vsini_max =  x[np.where(y == y.max())][0]
            vsini_err    =  np.std(vsini_max)
            plt.title('Distribution - PDF - Rotational velocity \n $\\nu \\sin i_p$ = '+ (str(round(vsini_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.show()
        else : 
            plt.figure()
            y, x, _ = plt.hist(self.vsini.value, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
            vsini_max =  x[np.where(y == y.max())][0]
            vsini_err    =  np.std(vsini_max)
            plt.title('Distribution - PDF - Rotational velocity \n $\\nu \\sin i_p$ = '+ (str(round(vsini_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.plot(velocities,vsini_pdf,color='#666666',label='PDF of $R$')
            plt.legend()   
            plt.show()

    def get_start_inclination(self,arg):
        
        # Compute omega_o pdf
        angles = np.linspace(0,180,1000)
        omega_o_kde = self.get_pdf(self,'star_inclination')
        omega_o_pdf = omega_o_kde(angles)

        if arg=='pdf':
            plt.figure()
            plt.plot(angles,omega_o_pdf,color='#666666',label='PDF of $R$')
            plt.xlabel('Degree (°)')
            plt.ylabel('PDF - Star inclination \n $\Omega_o$ = '+ (str(round(angles[np.argmax(omega_o_pdf)],2))) + '°')
            plt.legend()
        
        elif arg=='distribution':
            plt.figure()
            y, x, _ = plt.hist(self.vsini.value, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
            vsini_max =  x[np.where(y == y.max())][0]
            vsini_err    =  np.std(vsini_max)
            plt.title('Distribution - PDF - Rotational velocity \n $\\nu \\sin i_p$ = '+ (str(round(vsini_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.xlabel('Velocity ($km.s^{-1}$)')
        else : 
            plt.figure()
            y, x, _ = plt.hist(self.vsini.value, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
            vsini_max =  x[np.where(y == y.max())][0]
            vsini_err    =  np.std(vsini_max)
            plt.title('Distribution - PDF - Rotational velocity \n $\\nu \\sin i_p$ = '+ (str(round(vsini_max,2)))+ '$\pm$'+ str(radius_err) + '°')
            plt.xlabel('Velocity ($km.s^{-1}$)')
            plt.plot(velocities,vsini_pdf,color='#666666',label='PDF of $R$')
            plt.legend()   

    def set_orbital_inclination(self, io):
        if io[-1]>3.15:
            io=np.rad2deg(io)
        self.io = io

    def set_radius(self, radius):
        radius = radius * u.Rjup
        self.radius = radius

    def set_rot_velocity(self, vsini):
        vsini = vsini * u.km/s
        self.vsini = vsini
    
    def set_start_inclination(self, omega_o):
        if omega_o[-1]<=np.pi:
            omega_o=np.rad2deg(omega_o)
        self.omega_o = omega_o

    def compute_spin_axis(self):
        # Period limitation due to centrifugal force and gravitationnal force
        P_sample = (self.P > self.P_lim)

        # Set vel and visini with P > P_limit condition and v < v_limit
        self.vsini.value=self.vsini.value[P_sample]
        self.velocity=self.velocity[P_sample]

        v_sample= (self.vsini < self.velocity) 
        sin_ip = self.vsini/self.velocity

        self.ip=np.rad2deg(np.arcsin(sin_ip[v_sample]))
        n_sample=2*np.sum(v_sample)

    def get_spin_axis(self,arg,method):

        # Compute ip pdf

        if method =='easy':
            angles = np.linspace(0,180,1000)
            ip_kde = self.get_pdf(self.ip,'orbital_inclination')
            ip_pdf = ip_kde(angles)

            if arg=='pdf':
                plt.figure()
                plt.plot(angles,ip_pdf,color='green',label='PDF of $i_p$')
                plt.xlabel('Degree (°)')
                plt.ylabel('PDF - Spin axis \n $i_p$ = '+ (str(round(angles[np.argmax(ip_pdf)],2))) + '°')
                plt.legend()
        
            elif arg=='distribution':
                plt.figure()
                y, x, _ = plt.hist(self.ip, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
                ip_max =  x[np.where(y == y.max())][0]
                ip_err    =  np.std(ip_max)
                plt.title('Distribution - Spin axis \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
                plt.xlabel('Degree (°)')
            else : 
                plt.figure()
                y, x, _ = plt.hist(self.io, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
                ip_max =  x[np.where(y == y.max())][0]
                ip_err    =  np.std(ip_max)
                plt.title('Spin axis  \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
                plt.xlabel('Degree (°)')
                plt.plot(angles,ip_pdf,color='green',label='PDF of $i_p$')
                plt.xlabel('Degree (°)')
                plt.legend()

        elif method == 'complex':
            self.velocity=self.velocity.to(u.km/u.s)

            ### Initialization of important parameters
            Lv = gaussian_kde(self.velocity.value,bw_method='scott')
            Lu = gaussian_kde(self.vsini.value,bw_method='scott')
            v_range = np.linspace(0,self.v_limit.value,sizeof(self.velocity))         # Velocity interval [0;v_limit]
            n_points = 100                                          # The more the point we take the more the integral value is accurate, however computing time is longer. 
            angles_rad = np.linspace(0,np.pi,n_points)              # Angle interval [0;pi]
            cos_ip_complex_pdf = np.zeros_like(angles_rad)

            ### Integral calculation
            for k, cos_k in enumerate (np.cos(angles_rad)):
                int_dv = Lv(v_range)*Lu(v_range*np.sqrt(1-cos_k*cos_k))
                cos_ip_complex_pdf[k] = np.trapz(int_dv,v_range)

            ### Normalization of cos_ip PDF
            cos_ip_complex_pdf /= np.trapz(cos_ip_complex_pdf,angles_rad)

            ### PDF of ip
            ip_complex_rad_pdf = cos_ip_complex_pdf*np.abs(np.sin(angles_rad)) 

            ### Normalization of ip
            angles_comp = angles_rad*180/np.pi                       # In order to make a degree angle interval [0 ; 180]°
            ip_complex_pdf = ip_complex_rad_pdf/np.trapz(ip_complex_rad_pdf,angles_comp)

            if arg=='pdf':
                plt.figure()
                plt.plot(angles_comp,ip_complex_pdf,color='green',label='PDF of $i_p$')
                plt.xlabel('Degree (°)')
                plt.ylabel('PDF - Spin axis \n $i_p$ = '+ (str(round(angles[np.argmax(ip_pdf)],2))) + '°')
                plt.legend()
        
            elif arg=='distribution':
                plt.figure()
                y, x, _ = plt.hist(self.ip, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
                ip_max =  x[np.where(y == y.max())][0]
                ip_err    =  np.std(ip_max)
                plt.title('Distribution - Spin axis \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
                plt.xlabel('Degree (°)')
            else : 
                plt.figure()
                y, x, _ = plt.hist(self.io, bins=bins, density=True, color='#50C878',label='Distribution of $i_p$')
                ip_max =  x[np.where(y == y.max())][0]
                ip_err    =  np.std(ip_max)
                plt.title('Spin axis  \n $i_p$ = '+ (str(round(ip_max,2)))+ '$\pm$'+ str(ip_err) + '°')
                plt.xlabel('Degree (°)')
                plt.plot(angles_comp,ip_complex_pdf,color='green',label='PDF of $i_p$')
                plt.xlabel('Degree (°)')
                plt.legend()

        else :
            return "You must choose a method (easy/complex)"     
        
    #def get_true_obliquity(self,arg,method):

    #def get_projected_obliquity(self,arg,method):

    #def plot_all(self):




        

    