######################################################################
####################### Useful for Impedance #########################
######################################################################
import numpy as np
import pandas as pd
from scipy import stats,optimize,constants
def preprocess_Impedance_data(df,f):
    """Preprocess ZimT output data to remove the non stabilized
    part of the AC signal.

    Parameters
    ----------
    df : DataFrame
        Dataframe contening the tj_file output from ZimT

    f : float
        frequency of the measurement in Hz

    Returns
    -------
    df : DataFrame
        Dataframe contening the tj_file output from ZimT with
        the first period removed to get the stabilized AC signal.
    
    """
    df['t'] = df['t']-1/f
    df = df[df.t>0]
    return df

def sin_func(t,ampli,f,phi,offset):
    """Define Sine wave function,
    f(t) = ampli*np.sin(2*np.pi*freq*t + phi) + offset

    Parameters
    ----------
    t : 1-D sequence of floats
        contains the time vector

    ampli : float
        amplitude of the sine perturbation

    f : float
        frequency of the signal in Hz

    phi : float
        phase shift (in rad)

    offset : float
        offset of the signal

    Returns
    -------
    1-D sequence of floats
        Output the sine wave.
    """    
    return ampli*np.sin(2*np.pi*f*t + phi) + offset

def fit_sin_func(t,data,f):
    """Fit sin_func to the signal

    Parameters
    ----------
    t : 1-D sequence of floats
        contains the time data

    data : 1-D sequence of floats
        contains the data to fit

    f : float
        frequency of the measurement in Hz

    Returns
    -------
    ampli : float
        fitted amplitude of the sine perturbation

    freq : float
        fitted frequency of the signal in Hz

    phi : float
        fitted phase shift (in rad)

    offset : float
        fitted offset of the signal
    """    
    guess_offset = np.mean(data)
    guess_ampli = max(data) - guess_offset
    guess_freq = f
    guess_phi = 0

    ampli, freq, phi, offset = optimize.curve_fit(sin_func,t,data ,[guess_ampli, guess_freq, guess_phi, guess_offset],maxfev=int(1e5))[0]
    
    return ampli, freq, phi, offset 

def get_complex_impedance(df,f):
    """Estimate complex impedance from Impedance transient signal

    Parameters
    ----------
    df : DataFrame
        Dataframe containing the tj_file output from ZimT

    f : float
        frequency of the measurement in Hz

    Returns
    -------
    complex,
        Complex impedance
    """    
    # ampli_Va, freq_Va, phi_Va, offset_Va = fit_sin_func(np.asarray(df['t']),np.asarray(df['Va']),f)
    ampli_Va, freq_Va, phi_Va, offset_Va = fit_sin_func(np.asarray(df['t']),np.asarray(df['Vext']),f)
    ampli_Vdev, freq_Vdev, phi_Vdev, offset_Vdev = fit_sin_func(np.asarray(df['t']),np.asarray(df['Vext']),f)
    ampli_Jdev, freq_Jdev, phi_Jdev, offset_Jdev = fit_sin_func(np.asarray(df['t']),np.asarray(df['Jext']),f)
    phi = phi_Jdev - (phi_Vdev - phi_Va)

    return (ampli_Vdev/ampli_Jdev)*complex(np.cos(-phi),np.sin(-phi)) 

