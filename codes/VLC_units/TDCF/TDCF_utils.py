######################################################################
###################### Numerical eq for TDCF #########################
######################################################################


def TDCF_fit(x,k2,nBG):
    """Equation to fit for TDCF experiment
    dn/dt = -k2 * (ncol**2 + 2* ncol *nBG)
    see equation 4 in doi: 10.1063/1.5129037

    Parameters
    ----------
    x : 1-D sequence of floats
        Array contaning the ncol values.

    k2 : float 
        Bimolecular/radiative recombination constant.
        R = k2*n*p
    
    nBG : float
        Background carrier density.

    Returns
    -------
    y : 1-D sequence of floats
        Array of dn/dt
    """
    y = -k2*(x**2 + 2*nBG*x)

    return y