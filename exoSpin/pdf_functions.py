'''
ExoSpin - PDF important functions


@authors : I. Abdoulwahab & P. Palma-Bifani & G. Chauvin & A. Simonnin

'''

# ------------------------------------------------------------------------------------------------------------------------------------------------
## Imports

import matplotlib.pyplot as plt

import numpy as np

import astropy.units as u
import astropy.constants as c

from scipy import interpolate
from scipy.signal import savgol_filter
from scipy.stats import gaussian_kde
from scipy.stats import uniform



# ------------------------------------------------------------------------------------------------------------------------------------------------
## Important functions 

def kde(data):
    """
    From your data, the function returns a 1D Kernel Density Estimation.

    Args:
        data (numpy.ndarray): A 1D array that contains the data that will be estimated.

    Returns:
        (scipy.stats.gaussian_kde): A 1D Kernel Density Estimation object.

    Raises:
        ValueError: If data are not a 1D array of if data is empty.
    """

    if data.ndim != 1:
        raise ValueError("The data must be in a 1D array.")

    return gaussian_kde(data,bw_method='scott')

def pdf(kde, domain):
    """
    Evaluate the Kernel Density Estimation (KDE) over a specified interval and return a normalized PDF.

    Args:
        kde (scipy.stats.gaussian_kde): A 1D Kernel Density Estimation object.
        domain (numpy.ndarray): 1D array representing the domain where the KDE will be evaluated.

    Returns:
        (numpy.ndarray): 1D array of KDE values evaluated at each point in the interval and normalized.
    
    Raises:
        ValueError: If interval is not a 1D array.
    """

    if domain.ndim != 1:
        raise ValueError("The domain must be in a 1D array.")
    
    pdf = kde(domain)
    pdf /= np.trapz(pdf,domain)

    return pdf


def ip_complex_pdf(v_kde,vsini_kde,v_range,n):
    """
    Evaluate the companion spin axis PDF by using M. Bryan et al. 2020 method.

    Args:
        v_kde (scipy.stats.gaussian_kde): A 1D KDE of companion velocity.
        vsini_kde (scipy.stats.gaussian_kde): A 1D KDE of companion rotational velocity.
        v_range (numpy.ndarray): 1D array representing the domain where the velocities will be evaluated.
        n (int): Number of evaluated points.

    Returns:
        (numpy.ndarray): 1D array representing the PDF of companion spin axis.
    
    Raises:
        ValueError: The number of evaluted points must be greater than 1 and v_range must be a 1D array.
    """

    if v_range.ndim != 1:
        raise ValueError("The velocity domain must be in a 1D array.")

    if n <= 1:
        raise ValueError("The number of evaluted points must be greater than 1")

    angles_rad = np.linspace(0,np.pi,n)                                              
    cos_ip_pdf = np.zeros_like(angles_rad)
    ### Integral calculation
    for k, cos_k in enumerate (np.cos(angles_rad)):
        int_dv = v_kde(v_range)*vsini_kde(v_range*np.sqrt(1-cos_k*cos_k))
        cos_ip_pdf[k] = np.trapz(int_dv,v_range)
    ### Normalization of cos_ip PDF
    cos_ip_pdf /= np.trapz(cos_ip_pdf,angles_rad)
    ### PDF of ip
    ip_pdf = cos_ip_pdf*np.abs(np.sin(angles_rad)) 
    ### Normalization of ip
    angles = angles_rad*180/np.pi                                                           
    ip_pdf /= np.trapz(ip_pdf,angles)
    return ip_pdf

def proj_obli_complex_pdf(io_kde,ip_pdf,n):
    """
    Evaluate the companion projected obliquity PDF.

    Args:
        io_kde (scipy.stats.gaussian_kde): A 1D KDE of orbital inclination.
        ip_pdf (numpy.ndarray): 1D array of KDE evaluated values of the companion spin axis
        n (int): Number of evaluated points

    Returns:
        (numpy.ndarray): 1D array representing the PDF of companion spin axis.
    
    Raises:
        ValueError: The number of evaluted points must be greater than 1
    """

    if n <= 1:
        raise ValueError("The number of evaluted points must be greater than 1")

    Lio = io_kde
    angles_rad = np.linspace(0,np.pi,n)
    proj_obli_pdf = np.zeros_like(angles_rad)
    ### Integral calculation
    for k, ang_k in enumerate (angles_rad):
        int_ = ip_pdf*(Lio(angles_rad-ang_k)+Lio(angles_rad+ang_k))
        proj_obli_pdf[k] = np.trapz(int_,angles_rad)
    # Normalization 
    angles = angles_rad*180/np.pi
    proj_obli_pdf /= np.trapz(proj_obli_pdf,angles)
    return proj_obli_pdf

def true_obli_complex_pdf(io_pdf,ip_pdf,lambda_pdf,nb):
    """
    Evaluate the companion true obliquity PDF.

    Args:
        io_pdf (numpy.ndarray): 1D array of KDE evaluated values of the orbital inclination.
        ip_pdf (numpy.ndarray): 1D array of KDE evaluated values of the companion spin axis.
        omega_o_pdf (numpy.ndarray): 1D array of KDE evaluated values of the star inclination.
        nb (int): Number of evaluated points.

    Returns:
        (numpy.ndarray): 1D array representing the PDF of companion true obliquity.
    
    Raises:
        ValueError: The number of evaluted points must be greater than 1.
    """
    

    if nb <= 1:
        raise ValueError("The number of evaluted points must be greater than 1")

    bins = 200
    angles = np.linspace(0,180,nb)
    angles_ip = np.linspace(0,180,100)

    io_samp = np.random.choice(a=angles, p=io_pdf/np.sum(io_pdf),size=nb)
    ip_samp = np.random.choice(a=angles_ip, p=ip_pdf/np.sum(ip_pdf),size=nb)
    lambda_samp = np.random.choice(a=angles, p=lambda_pdf/np.sum(lambda_pdf),size=nb)

    psi_samp = np.zeros_like(io_samp)

    for i in range (nb):    
        arg_1 = np.cos(np.deg2rad(io_samp[i]))*np.cos(np.deg2rad(ip_samp[i]))
        arg_2 = np.sin(np.deg2rad(io_samp[i]))*np.sin(np.deg2rad(ip_samp[i]))*np.cos(np.deg2rad(lambda_samp[i]))
        psi_samp[i] = np.arccos(arg_1+arg_2)

    psi_samp = np.rad2deg(psi_samp)

    psi_kde = gaussian_kde(psi_samp,bw_method='scott')
    psi_pdf = psi_kde(angles)
    
    return psi_pdf