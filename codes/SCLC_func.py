###################################################
############### Useful function ###################
###################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_settings_screen
from scipy import stats,optimize,constants
import subprocess,shutil,os,glob
from itertools import repeat
from lmfit import Model
import warnings
import sys
# Don't show warnings
warnings.filterwarnings("ignore")
## Physics constants
q = constants.value(u'elementary charge')
eps_0 = constants.value(u'electric constant')
kb = constants.value(u'Boltzmann constant in eV/K')

#############################################################
################ SCLC analysis functions ####################
#############################################################

def MottGurney(V,mu,eps_r,Vbi,L):
    """Mott-Gurney equation typically used to fit single carrier devices JVs to get the mobility.

    J = (9/8)*eps_0*eps_r*mu*((V-Vbi)**2)/(L**3)

    Refs:
    Mott N F and Gurney R W 1940 Electronic Processes in Ionic Crystals (Oxford: Oxford University Press)
    Nice paper to read on MG ==> J. Phys.: Condens. Matter 30 (2018) 105901

    Parameters
    ----------
    V : 1-D sequence of floats
        Array containing the voltages (V)
    mu : float
        Mobility (m^2/(Vs))
    eps_r : float
        Relative dielectric constant
    Vbi : float
        Built-in voltage (V)
    L : float
        Thickness (m)

    Returns
    -------
    1-D sequence of floats
        Array containing the current-density (A m^-2)
    """    
    return (9/8)*eps_0*eps_r*mu*((V-Vbi)**2)/(L**3)


def fit_MottGurney(V,J,mu,eps_r,Vbi,L,var2fit=['mu']):
    """[summary]

    Parameters
    ----------
    V : 1-D sequence of floats
        Array containing the voltages (V)
    J : 1-D sequence of floats
        Array containing the current-density (A m^-2)
    mu : float
        Mobility (m^2/(Vs))
    eps_r : float
        Relative dielectric constant
    Vbi : float
        Built-in voltage (V)
    L : float
        Thickness (m)
    var2fit : list, optional
        Chose variable to fit (can be multiple), by default ['mu']

    Returns
    -------
    list
        list of the fitted values in the same order than var2fit
    """    
    fixed_var = ['V','mu','eps_r','Vbi','L']
    for i in var2fit:
        if i in fixed_var:
            fixed_var.remove(i)
    
    Mod_MottGurney = Model(MottGurney, independent_vars=fixed_var)
    result = Mod_MottGurney.fit(J, V=V,mu=mu,eps_r=eps_r,Vbi=Vbi,L=L)
   
    best_fit = []
    for i in var2fit:
        best_fit.append(result.params[i].value)
    return best_fit

def deriv(x,y):
    """Get the derivative of input data dy/dx

    Parameters
    ----------
    x : 1-D sequence of floats
        x-axis data
    y : 1-D sequence of floats
        y-axis data

    Returns
    -------
    1-D sequence of floats
        dy/dx
    """    
    dy = np.diff(y)/np.diff(x)
    dum = (y[-1] - y[-2])/(x[-1] - x[-2]) 
    dy = np.append(dy,[dum])
    return list(dy)

def log_slope(x,y):
    """Get the slope of the logarithmic data

    Parameters
    ----------
    x : 1-D sequence of floats
        x-axis data
    y : 1-D sequence of floats
        y-axis data

    Returns
    -------
    1-D sequence of floats
        slope of the logarithmic data
    """ 
    # Convert data into np.arrays
    x = np.asarray(x)
    y = np.asarray(y)   
    # remove V = 0V efore getting log
    if 0 in x:
        idx_0 = x.index(0)
        x = x.pop(idx_0)
        y = y.pop(idx_0)
    log_xf = np.log(x)
    log_yf = np.log(y)
    # Take derivative of data
    slopef = deriv(log_xf,log_yf)
    return np.asarray(slopef)

def calc_vnet_with_ions(ions,traps,L,eps_r):
    """Calculate Vnet for SCLC measurement in the presence of Ions
    as defined in equation (4) in ACS Energy Lett. 2021, 6, 3, 1087–1094
    https://doi.org/10.1021/acsenergylett.0c02599

    Vnet = q * n_net * L**2 /( 2 * eps_0 * eps_r) = q * (traps-ions) * L**2 /( 2 * eps_0 * eps_r)

    Parameters
    ----------
    ions : float
        Ion density (m^-3)
    traps : float
        Trap density (m^-3)
    L : float
        Thickness (m)
    eps_r : float
        Relative dielectric constant

    Returns
    -------
    float
        returns vnet as defined in 
    """    
    vnet = q * (traps-ions) * L**2 /( 2 * eps_0 * eps_r)
    return vnet

def calc_net_charge(Vnet,L,eps_r):
    """Calculate the net charge for SCLC measurement
    as defined in equation (4) in ACS Energy Lett. 2021, 6, 3, 1087–1094
    https://doi.org/10.1021/acsenergylett.0c02599

    Vnet = q * n_net * L**2 /( 2 * eps_0 * eps_r) 

    Note that if there are not ions or dopant n_net = traps

    Parameters
    ----------
    Vnet : float
        Vnet see ref
    L : float
        Thickness (m)
    eps_r : float
        Relative dielectric constant

    Returns
    -------
    float
        net-charge in the bulk
    """    
    net_charge= Vnet * ( 2 * eps_0 * eps_r)/ (q * L**2)
    return net_charge # m^-3

def SCLC_get_data_plot(volt,curr):
    # Get loglog graph slopes and JV without 0V
    slopesf = log_slope(volt,curr)
    V_slopef = np.asarray(volt)
    J_slopef = np.asarray(curr)

    # Get max slope if max slopes > 2.2
    idx_maxf = np.argmax(slopesf)
    max_slopesf = slopesf[idx_maxf]
    
    if max_slopesf>2.0:
        get_tangentf = 1
    else: get_tangentf = 0


    # Calculate tangent and crossing point
    tang_val_V1f,tang_val_V2f,tang_val_V3f,V1f,J1f,V2f,J2f = [],[],[],[],[],[],[]
    if get_tangentf == 1:
            tang_val_V1f = np.exp(np.log(J_slopef[0])+1*(np.log(V_slopef)-np.log(V_slopef[0])))
            tang_val_V2f = np.exp(np.log(J_slopef[idx_maxf])+max_slopesf*(np.log(V_slopef)-np.log(V_slopef[idx_maxf])))
            tang_val_V3f = np.exp(np.log(J_slopef[-1])+2*(np.log(V_slopef)-np.log(V_slopef[-1])))
            V1f = np.exp((1/(1-max_slopesf)*np.log((J_slopef[idx_maxf]*V_slopef[0])/(J_slopef[0]*(V_slopef[idx_maxf]**max_slopesf)))))
            J1f = np.exp(np.log(J_slopef[0])+1*(np.log(V1f)-np.log(V_slopef[0])))
            V2f = np.exp((1/(2-max_slopesf)*np.log((J_slopef[idx_maxf]*(V_slopef[-1]**2))/(J_slopef[-1]*(V_slopef[idx_maxf]**max_slopesf)))))
            J2f = np.exp(np.log(J_slopef[-1])+2*(np.log(V2f)-np.log(V_slopef[-1])))

    # print(calc_net_charge(V1,thick))
    # print(calc_net_charge(V2,thick))

    return V_slopef,J_slopef,slopesf,get_tangentf,idx_maxf,max_slopesf,tang_val_V1f,tang_val_V2f,tang_val_V3f,V1f,J1f,V2f,J2f


