'''
ExoSpin run script - PDF important functions


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


# ------------------------------------------------------------------------------------------------------------------------------------------------
## Important functions 

def kde(data):
    return gaussian_kde(data,bw_method='scott')

def pdf(data_kde,interval):
    return data_kde(interval)

def cos_pdf(pdf,domain):
    sin_domain = np.abs(np.sin(domain))

    for i in range (len(sin_domain)):
        if sin_domain[i]==0:
            sin_domain[i]=sin_domain[i+1]

    cos_pdf = pdf/sin_domain
    cos_pdf /= np.trapz(cos_pdf,domain)
    print(np.trapz(cos_pdf,domain))
    return cos_pdf

def sin_pdf(pdf,domain):
    cos_domain = np.abs(np.cos(domain))

    for i in range (len(cos_domain)):
        if cos_domain[i]==0:
            cos_domain[i]=sin_domain[i+1]

    sin_pdf = pdf/cos_domain
    sin_pdf /= np.trapz(sin_pdf,domain)
    return sin_pdf

def arccos_pdf(pdf,domain):
    arccos_domain = np.sqrt(1-domain*domain)
    arccos_pdf = pdf*arccos_domain
    arccos_pdf /= np.trapz(arccos_pdf,domain)
    return arccos_pdf

def product_pdf(pdf_1 , pdf_2 ,pdf_domain):
    prd_domain = pdf_domain
    new_pdf_1 = interp1d(pdf_domain,pdf_1,kind='linear', fill_value="extrapolate")
    new_pdf_2 = interp1d(pdf_domain,pdf_2,kind='linear', fill_value="extrapolate")
    prd_pdf = np.zeros_like(pdf_domain)
    for i, y in enumerate (prd_domain):
        int = new_pdf_1(y/pdf_domain)*new_pdf_2(pdf_domain)/np.abs(pdf_domain)
        prd_pdf[i]= np.trapz(int, pdf_domain)
    prd_pdf /= np.trapz(prd_pdf,pdf_domain)   
    return prd_pdf

def sum_pdf(pdf_1 , pdf_2 ,pdf_domain):
    sum_domain = pdf_domain
    new_pdf_1 = interp1d(pdf_domain,pdf_1,kind='linear', fill_value="extrapolate")
    new_pdf_2 = interp1d(pdf_domain,pdf_2,kind='linear', fill_value="extrapolate")
    sum_pdf = np.zeros_like(pdf_domain)
    for i, y in enumerate (sum_domain):
        int = new_pdf_1(pdf_domain)*new_pdf_2(y-pdf_domain)
        sum_pdf[i]= np.trapz(int, pdf_domain)
    return sum_pdf

def ip_complex_pdf(v_kde,vsini_kde,v_range,n_points):
    angles_rad = np.linspace(0,np.pi,n_points)                                              
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


def proj_obli_complex_pdf(io_kde,ip_pdf,n_points):
    Lio = io_kde
    angles_rad = np.linspace(0,np.pi,n_points)
    proj_obli_pdf = np.zeros_like(angles_rad)
    ### Integral calculation
    for k, ang_k in enumerate (angles_rad):
        int_ = ip_pdf*(Lio(angles_rad-ang_k)+Lio(angles_rad+ang_k))
        proj_obli_pdf[k] = np.trapz(int_,angles_rad)
    # Normalization 
    angles = angles_rad*180/np.pi
    proj_obli_pdf /= np.trapz(proj_obli_pdf,angles)
    return proj_obli_pdf

def true_obli_complex_pdf(io_kde,ip_pdf,n_points):
    Lio = io_kde
    angles_rad = np.linspace(0,np.pi,n_points)
    proj_obli_pdf = np.zeros_like(angles_rad)
    ### Integral calculation
    for k, ang_k in enumerate (angles_rad):
        int_ = ip_pdf*(Lio(angles_rad-ang_k)+Lio(angles_rad+ang_k))
        proj_obli_pdf[k] = np.trapz(int_,angles_rad)
    # Normalization 
    angles = angles_rad*180/np.pi
    proj_obli_pdf /= np.trapz(proj_obli_pdf,angles)
    return proj_obli_pdf