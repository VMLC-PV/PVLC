###################################################
############### Useful function ###################
###################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_settings_screen
from scipy import stats,optimize,constants, fftpack
import subprocess,shutil,os,glob
from itertools import repeat
import warnings
import sys
from pathlib import Path
# Don't show warnings
warnings.filterwarnings("ignore")
## Physics constants
q = constants.value(u'elementary charge')
eps_0 = constants.value(u'electric constant')
kb = constants.value(u'Boltzmann constant in eV/K')

def sci_notation(number, sig_fig=2):
    """Make proper scientific notation for graphs

    Parameters
    ----------
    number : float
        Number to put in scientific notation.

    sig_fig : int, optional
        Number of significant digits (Defaults = 2).

    Returns
    -------
    output : str
        String containing the number in scientific notation
    """
    if sig_fig != -1:
        if number == 0:
            output = '0'
        else:
            ret_string = "{0:.{1:d}e}".format(number, sig_fig)
            a,b = ret_string.split("e")
            if int(b) > 0:
                b = int(b) #removed leading "+" and strips leading zeros too.
                c = ''
            else: 
                b = abs(int(b))
                c = u"\u207B" # superscript minus sign
            SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
            b = str(b).translate(SUP)
            output =a + ' x 10' + c + b
    else:
        if number == 0:
            output = '0'
        else: 
            ret_string = "{0:.{1:d}e}".format(number, 0)
            a,b = ret_string.split("e")
            b = int(b) #removed leading "+" and strips leading zeros too.
            SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
            #SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
            b = str(b).translate(SUP)
            output = '10' + b    
    return output

def monoExp(t, k, A, B):
    """ Monoexponential decay function
    f(t) = A * np.exp(-k * t ) + B

    Parameters
    ----------
    t : 1-D sequence of floats
        time
    k : float
        decay rate
    A : float
        initial quantity
    B : float
        offset

    Returns
    -------
    1-D sequence of floats
        f(t)
    """
    return A * np.exp(-k * t ) + B

def get_Jsc(Volt,Curr):
    """Get the short-circuit current (Jsc) from solar cell JV-curve by interpolating the current at 0 V

    Parameters
    ----------
    Volt : 1-D sequence of floats
        Array containing the voltages.

    Curr : 1-D sequence of floats
        Array containing the current-densities.

    Returns
    -------
    Jsc : float
        Short-circuit current value
    """
    Jsc_dumb = np.interp(0, Volt, Curr)
    return Jsc_dumb

def get_Voc(Volt,Curr):
    """Get the Open-circuit voltage (Voc) from solar cell JV-curve by interpolating the Voltage when the current is 0

    Parameters
    ----------
    Volt : 1-D sequence of floats
        Array containing the voltages.

    Curr : 1-D sequence of floats
        Array containing the current-densities.

    Returns
    -------
    Voc : float
        Open-circuit voltage value
    """
    Voc_dumb = np.interp(0, Curr, Volt)
    return Voc_dumb

def get_FF(Volt,Curr):
    """Get the fill factor (FF) from solar cell JV-curve by calculating the maximum power point

    Parameters
    ----------
    Volt : 1-D sequence of floats
        Array containing the voltages.

    Curr : 1-D sequence of floats
        Array containing the current-densities.

    Returns
    -------
    FF : float
        Fill factor value
    """
    power = []
    Volt_oc = get_Voc(Volt,Curr)
    Curr_sc = get_Jsc(Volt,Curr)
    for i,j in zip(Volt,Curr):
        if (i < Volt_oc and j > Curr_sc):
            power.append(i*j)
    power_max = min(power)
    FF_dumb = power_max/(Volt_oc*Curr_sc)
    return FF_dumb

def get_PCE(Volt,Curr,suns):
    """Get the power conversion efficiency (PCE) from solar cell JV-curve

    Parameters
    ----------
    Volt : 1-D sequence of floats
        Array containing the voltages.

    Curr : 1-D sequence of floats
        Array containing the current-densities.

    Returns
    -------
    PCE : float
        Power conversion efficiency value.
    """
    Voc_dumb = get_Voc(Volt,Curr)
    Jsc_dumb = get_Jsc(Volt, Curr)
    FF_dumb = get_FF(Volt, Curr)
    PCE_dumb = Voc_dumb*Jsc_dumb*FF_dumb/suns
    return PCE_dumb  

def get_ideality_factor(suns,Vocs,T=295):
    """Returns ideality factor from suns-Voc data linear fit of Voc = (nIF/Vt)*log(suns) + intercept

    Parameters
    ----------
    suns : 1-D sequence of floats
        Array containing the intensity in sun.

    Vocs : 1-D sequence of floats
        Array containing the open-circuit voltages.
    
    T : float optional
        Temperature in Kelvin (Default = 295 K).

    Returns
    -------
    nIF : float
        Ideality factor value.
    intercept : float
        Intercept of the regression line.
    rvalue : float
        Correlation coefficient.
    pvalue : float
        Two-sided p-value for a hypothesis test whose null hypothesis is
        that the slope is zero, using Wald Test with t-distribution of
        the test statistic.
    stderr : float
        Standard error of the estimated gradient.
    """
    Vt = kb*T
    suns = np.log(suns)
    slope_d, intercept_d, r_value_d, p_value_d, std_err_d = stats.linregress(suns,Vocs)
    nIF = slope_d/Vt
    return nIF,intercept_d, r_value_d**2, p_value_d, std_err_d

def get_alpha_factor(suns,Jscs):
    """Returns alpha from suns-Jsc data linear fit of log(Jsc) = alpha*log(suns) + b

    Parameters
    ----------
    suns : 1-D sequence of floats
        Array containing the intensity in sun.

    Vocs : 1-D sequence of floats
        Array containing the open-circuit voltages.
    
    Returns
    -------
    alpha : float
        Alpha value.

    intercept : float
        Intercept of the regression line.

    rvalue : float
        Correlation coefficient.

    pvalue : float
        Two-sided p-value for a hypothesis test whose null hypothesis is
        that the slope is zero, using Wald Test with t-distribution of
        the test statistic.

    stderr : float
        Standard error of the estimated gradient.
    """
    suns = np.log(suns)
    Jscs = np.log(Jscs)
    alpha, intercept_d, r_value_d, p_value_d, std_err_d = stats.linregress(suns,Jscs)
    return alpha,intercept_d, r_value_d**2, p_value_d, std_err_d


######################################################################
############### Function to run simulation code ######################
######################################################################

def run_SIMsalabim_light_dep(str2run,G_frac,System):
    """Run SIMsalabim in the current folder for specified light intensity.

    Parameters
    ----------
    st2run : str
        String to run with SIMsalabim.

    Gfrac : float 
        Gfrac argument for SIMsalabim, light intensity.
    
    System : String,
        String to specify the operating system.
        The options are: 'Linux' or 'Windows'

    Returns
    -------
    
    """
    if System == 'Windows':
        subprocess.call(r'SIMsalabim.exe -ScPars_file scPars'+str(G_frac)+'.dat -Gfrac '+str(G_frac)+' '+'-JV_file JV_'+str(G_frac)+'.dat'+' '+str2run)
    elif System == 'Linux':
        subprocess.check_call(('./SIMsalabim -ScPars_file scPars'+str(G_frac)+'.dat -Gfrac '+str(G_frac)+' '+'-JV_file JV_'+str(G_frac)+'.dat'+' '+str2run).split())
    else: print('Wrong system input')

def run_SIMsalabim(str2run,System,path=''):
    """Run SIMsalabim in the folder specified by 'path'.

    Parameters
    ----------
    st2run : str
        String to run with SIMsalabim.
    
    System : str
        String to specify the operating system.
        The options are: 'Linux' or 'Windows'
    
    path : str
        path to the folder containing SIMsalabim 
        ('./SIMsalabim' in Linux and 'SIMsalabim.exe' in Windows ).

    Returns
    -------
    
    """
    curr_dir = os.getcwd()
    os.chdir(path)
    FNULL = open(os.devnull, 'w')
    # stdout=FNULL, stderr=subprocess.STDOUT is used to turn of the display of SIMsalabim output on the terminal
    if System == 'Windows':
        subprocess.call(path+r'SIMsalabim.exe ' + str2run,stdout=FNULL, stderr=subprocess.STDOUT)
    elif System == 'Linux':
        subprocess.check_call(('./SIMsalabim ' + str2run).split(),stdout=FNULL, stderr=subprocess.STDOUT)
    else: print('Wrong system input')
    os.chdir(curr_dir)
   
def run_zimt(str2run,System,path=''):
    """Run zimt in the folder specified by 'path'.

    Parameters
    ----------
    st2run : str
        String to run with SIMsalabim.
    
    System : str
        String to specify the operating system.
        The options are: 'Linux' or 'Windows'
    
    path : str
        path to the folder containing zimt 
        ('./zimt' in Linux and 'zimt.exe' in Windows ).

    Returns
    -------
    
    """
    curr_dir = os.getcwd()
    print(path)
    os.chdir(path)
    FNULL = open(os.devnull, 'w')
    # stdout=FNULL, stderr=subprocess.STDOUT is used to turn of the display of zimt output on the terminal
    if System == 'Windows':
        subprocess.call(path+r'zimt.exe ' + str2run,stdout=FNULL, stderr=subprocess.STDOUT)
    elif System == 'Linux':
        subprocess.check_call(('./zimt ' + str2run).split(),stdout=FNULL, stderr=subprocess.STDOUT)
    else: print('Wrong system input')
    os.chdir(curr_dir)

def run_code(name_prog,path,str2run,System):
    """Run program 'name_prog' in the folder specified by 'path'.

    Parameters
    ----------
    name_prog : str
        name of the program to run.

    st2run : str
        String to run for the name_prog.
    
    System : str
        String to specify the operating system.
        The options are: 'Linux' or 'Windows'
    
    path : str
        path to the folder containing zimt 
        ('./zimt' in Linux and 'zimt.exe' in Windows ).

    Returns
    -------
    
    """
    curr_dir = os.getcwd()
    os.chdir(curr_dir+path)
    if System == 'Windows':
        print(curr_dir+path+name_prog+'.exe ' + str2run)
        subprocess.call(curr_dir+path+name_prog+'.exe ' + str2run)
    elif System == 'Linux':
        subprocess.check_call(('./'+name_prog+' ' + str2run).split())
    else: print('Wrong system input')
    os.chdir(curr_dir)



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


######################################################################
####################### Useful for Impedance #########################
######################################################################

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


######################################################################
#################### Plot function SIMsalabim ########################
######################################################################

def make_df_JV(path_JV_file):
    """Output panda dataframe containing JV_file

    Parameters
    ----------
    path_JV_file : str
        path to file containing the JV_file output from SIMsalabim
    
    Returns
    -------
    DataFrame
        panda dataframe containing JV_file
    """    
    #names = ['Vext','Vint','Jext','Jint','P','recLan','recSRH','Jbimo','JSRH_bulk','JSRH_LI','JSRH_RI','Jph','Jn_l','Jp_l','Jn_r','Jp_r']
    df = pd.read_csv(path_JV_file,delim_whitespace=True)

    return df

def make_df_Var(path_Var_file):
    """Output panda dataframe containing Var_file

    Parameters
    ----------
    path_JV_file : str
        path to file containing the Var_file output from SIMsalabim
    
    Returns
    -------
    DataFrame
        panda dataframe containing JV_file
    """   
    #names = ['x','V','n','p','Evac','Ec','Ev','phin','phip','ntrap','ptrap','nid','pid','nion','pion','mun','mup','rec','dp','Gm','Jn','Jp']
    df = pd.read_csv(path_Var_file,delim_whitespace=True)

    return df

def SIMsalabim_nrj_diag(num_fig,data_var,th_TL_left,th_TL_right,pic_save_name,save_yes):
    """"Make energy diagram plot from Var_file output SIMsalabim

    Parameters
    ----------
    num_fig : int
        figure number where to plot the energy diagram
    data_var : DataFrame
        Panda Dataframe containing the Var_file output from SIMsalabim (see function "make_df_Var")
    th_TL_left : float
        Thickness of the left transport layer
    th_TL_right : float
        Thickness of the right transport layer
    pic_save_name : str
        name of the file where the figure is saved
    save_yes : bool
        If True, save energy diagram as an image with the  file name defined by "pic_save_name" 
    """    
    
    line_thick = 3
    # Color control for the different layers
    color_nTL = '#c7e6a3'
    color_pTL = '#8cb5f0'
    color_pero = '#ba0000'
    color_electrode ='#999999'

    # Plotting
    plt.figure(num_fig)
    ax_nrj_diag = plt.axes()
    # ax_nrj_diag.plot('x','Evac',data=data_var,label = r'E$_{vac}$',linestyle='-',linewidth=2,color = 'k')
    ax_nrj_diag.plot('x','Ec',data=data_var,label = r'E$_{c}$',linestyle='-', linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','Ev',data=data_var,label = r'E$_{v}$',linestyle='-', linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','phin',data=data_var,label = r'E$_{fn}$',linestyle='--',linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','phip',data=data_var,label = r'E$_{fp}$',linestyle='--',linewidth=line_thick,color = 'k')
    TL_left = data_var[data_var['x']<th_TL_left]
    TL_right = data_var[data_var['x']>max(data_var['x'])-th_TL_right]
    AL = data_var[data_var['x']<max(data_var['x'])-th_TL_right]
    AL = AL[AL['x']>th_TL_left]
    ax_nrj_diag.fill_between(TL_left['x'],TL_left['Ec'],y2=0,color=color_nTL)
    ax_nrj_diag.fill_between(TL_left['x'],TL_left['Ev'],y2=-8,color=color_nTL)
    ax_nrj_diag.fill_between(TL_right['x'],TL_right['Ec'],y2=0,color=color_pTL)
    ax_nrj_diag.fill_between(TL_right['x'],TL_right['Ev'],y2=-8,color=color_pTL)
    ax_nrj_diag.fill_between(AL['x'],AL['Ec'],y2=0,color=color_pero)
    ax_nrj_diag.fill_between(AL['x'],AL['Ev'],y2=-8,color=color_pero)
    ax_nrj_diag.plot([-10,0],[min(data_var['phin']),min(data_var['phin'])],color='k')
    ax_nrj_diag.plot([max(data_var['x']),max(data_var['x'])+10],[max(data_var['phip']),max(data_var['phip'])],color='k')
    ax_nrj_diag.fill_between([-10,0],[min(data_var['phin']),min(data_var['phin'])],y2=-8,color=color_electrode)
    ax_nrj_diag.fill_between([max(data_var['x']),max(data_var['x'])+10],[max(data_var['phip']),max(data_var['phip'])],y2=-8,color=color_electrode)
    # plt.axhline(y=max(data_var['phip']), color='k', linestyle='-')

    # Hide axis and spines
    ax_nrj_diag.get_xaxis().set_visible(False)
    ax_nrj_diag.get_yaxis().set_visible(False)
    for sides in ['right','left','top','bottom']:
        ax_nrj_diag.spines[sides].set_visible(False)

    # Legend
    plt.legend(loc='center',frameon=False,ncol = 2, bbox_to_anchor=(0.52,1.02),fontsize = 40)
    plt.tight_layout()

    # Save file
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)

def SIMsalabim_JVs_plot(num_fig,data_JV,x='Vext',y=['Jext'],xlimits=[],ylimits=[],plot_type=0,labels='',colors='b',line_type = ['-'],mark='',legend=True,plot_jvexp=False,data_JVexp=pd.DataFrame(),save_yes=False,pic_save_name='JV.jpg'):
    """ Make JV_plot for SIMsalabim  
    
    Parameters
    ----------
    num_fig : int
        number of the fig to plot JV
    data_JV : DataFrame
        Panda DataFrame containing JV_file
    x : str, optional
        xaxis data  (default = 'Vext'), by default 'Vext'
    y : list of str, optional
        yaxis data can be multiple like ['Jext','Jbimo']  (default = ['Jext']), by default ['Jext']
    xlimits : list, optional
        x axis limits if = [] it lets python chose limits, by default []
    ylimits : list, optional
        y axis limits if = [] it lets python chose limits, by default []
    plot_type : int, optional
        type of plot 1 = logx, 2 = logy, 3 = loglog else linlin (default = linlin), by default 0
    labels : str, optional
        label of the JV, by default ''
    colors : str, optional
        color for the JV line, by default 'b'
    line_type : list, optional
        type of line for simulated data plot
        size line_type need to be = size(y), by default ['-']
    mark : str, optional
        type of Marker for the JV, by default ''
    legend : bool, optional
        Display legend or not, by default True
    plot_jvexp : bool, optional
        plot an experimental JV or not, by default False
    data_JVexp : [type], optional
        Panda DataFrame containing experimental JV_file with 'V' the voltage and 'J' the current, by default pd.DataFrame()
    save_yes : bool, optional
        If True, save JV as an image with the  file name defined by "pic_save_name", by default False
    pic_save_name : str, optional
        name of the file where the figure is saved, by default 'JV.jpg'
    """    
        
    if len(y) > len(line_type):
        print('Invalid line_type list, we meed len(y) == len(line_type)')
        print('We will use default line type instead')
        line_type = []
        for counter, value in enumerate(y):
            line_type.append('-')

    plt.figure(num_fig)
    ax_JVs_plot = plt.axes()
    for i,line in zip(y,line_type):
        if plot_type == 1:
            ax_JVs_plot.semilogx(data_JV['Vext'],data_JV[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            if plot_jvexp:
                ax_JVs_plot.semilogx(data_JVexp['V'],data_JVexp['J']/10,'o',markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        
        elif plot_type == 2:
            ax_JVs_plot.semilogy(data_JV['Vext'],data_JV[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)   
            if plot_jvexp:
                ax_JVs_plot.semilogy(data_JVexp['V'],data_JVexp['J']/10,'o',markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 3:
            ax_JVs_plot.loglog(data_JV['Vext'],data_JV[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            if plot_jvexp:
                ax_JVs_plot.loglog(data_JVexp['V'],data_JVexp['J']/10,'o',markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        else:
            ax_JVs_plot.plot(data_JV['Vext'],data_JV[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            if plot_jvexp:
                ax_JVs_plot.plot(data_JVexp['V'],data_JVexp['J']/10,'o',markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
        
    
    # legend
    if legend == True:
        plt.legend(loc='best',frameon=False,fontsize = 30)
    if xlimits != []:
        plt.xlim(xlimits)
    if ylimits != []:
        plt.ylim(ylimits)
    plt.grid(b=True,which='both')
    plt.xlabel('Applied Voltage [V]')
    plt.ylabel('Current Density [mA cm$^{-2}$]')
    plt.tight_layout()
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)

def SIMsalabim_dens_plot(num_fig,data_Var,Vext='nan',y=['n','p'],xlimits=[],ylimits=[],plot_type=0,labels='',colors='b',line_type = ['-'],legend=True,save_yes=False,pic_save_name='density.jpg'):
    """Make Var_plot for SIMsalabim 

    Parameters
    ----------
    num_fig : int
        number of the fig to plot JV
    data_JV : DataFrame
        Panda DataFrame containing JV_file
    Vext : float
        float to define the voltage at which the densities will be plotted if Vext='nan' then take Vext as max(Vext), if Vext does not exist we plot the closest voltage, default 'nan'
    y : list of str, optional
        yaxis data can be multiple like ['n','p'], by default ['n','p']
    xlimits : list, optional
        x axis limits if = [] it lets python chose limits, by default []
    ylimits : list, optional
        y axis limits if = [] it lets python chose limits, by default []
    plot_type : int, optional
        type of plot 1 = logx, 2 = logy, 3 = loglog else linlin (default = linlin), by default 0
    labels : str, optional
        label of the line, by default ''
    colors : str, optional
        color for the line, by default 'b'
    line_type : list, optional
        type of line for the plot
        size line_type need to be = size(y), by default ['-']
    legend : bool, optional
        Display legend or not, by default True
    save_yes : bool, optional
        If True, save density plot as an image with the  file name defined by "pic_save_name", by default False
    pic_save_name : str, optional
        name of the file where the figure is saved, by default 'density.jpg'
    """    
    
    if len(y) > len(line_type):
        print('Invalid line_type list, we meed len(y) == len(line_type)')
        print('We will use default line type instead')
        line_type = []
        for counter, value in enumerate(y):
            line_type.append('-')

    if Vext == 'nan':
        Vext = max(data_Var['Vext'])
        print('Vext was not specified so Vext = {:.2f} V was plotted'.format(Vext))
    
    data_Var = data_Var[abs(data_Var.Vext -Vext) == min(abs(data_Var.Vext -Vext))]
    data_Var = data_Var.reset_index(drop=True)
    if min(abs(data_Var.Vext -Vext)) != 0:
        print('Vext = {:.2f} V was not found so {:.2f} was plotted'.format(Vext,data_Var['Vext'][0]))
    
    # print(data_Var)
    plt.figure(num_fig)
    ax_Vars_plot = plt.axes()
    for i,line in zip(y,line_type):
        if plot_type == 1:
            ax_Vars_plot.semilogx(data_Var['x'],data_Var[i]/1e6,color=colors,label=labels,linestyle=line)
                       
        
        elif plot_type == 2:
            ax_Vars_plot.semilogy(data_Var['x'],data_Var[i]/1e6,color=colors,label=labels,linestyle=line)   
            
        elif plot_type == 3:
            ax_Vars_plot.loglog(data_Var['x'],data_Var[i]/1e6,color=colors,label=labels,linestyle=line)
            
        else:
            ax_Vars_plot.plot(data_Var['x'],data_Var[i]/1e6,color=colors,label=labels,linestyle=line)
            
    
    # legend
    if legend == True:
        # plt.legend(loc='best',frameon=False,fontsize = 30)
        from matplotlib.legend import Legend
        leg = Legend(ax_Vars_plot, ['-','--'], ['line C', 'line D'],
             loc='lower right', frameon=False)
        ax_Vars_plot.add_artist(leg);
    if xlimits != []:
        plt.xlim(xlimits)
    if ylimits != []:
        plt.ylim(ylimits)
    plt.grid(b=True,which='both')
    plt.xlabel('x [nm]')
    plt.ylabel('Density [cm$^{-3}$]')
    plt.tight_layout()
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)


#####################################################################
###################### Plot function ZimT ###########################
#####################################################################

def zimt_tj_plot(num_fig,data_tj,x='t',y=['Jext'],xlimits=[],ylimits=[],plot_type=0,labels='',colors='b',line_type = ['-'],mark='',legend=True,save_yes=False,pic_save_name='transient.jpg'):
    """ Make tj_file transient plot for ZimT  
    Default time on the x axis in $\mu$s
    
    Parameters
    ----------
    num_fig : int
        number of the fig to plot tj
    data_tj : DataFrame
        Panda DataFrame containing tj_file
    x : str, optional
        xaxis data, by default 't'
    y : list of str, optional
        yaxis data can be multiple like ['Jext','Jncat'], by default ['Jext']
    xlimits : list, optional
        x axis limits if = [] it lets python chose limits, by default []
    ylimits : list, optional
        y axis limits if = [] it lets python chose limits, by default []
    plot_type : int, optional
        type of plot 1 = logx, 2 = logy, 3 = loglog else linlin (default = linlin), by default 0
    labels : str, optional
        label of the tj, by default ''
    colors : str, optional
        color for the JV line, by default 'b'
    line_type : list, optional
        type of line used for the plot
        size line_type need to be = size(y), by default ['-']
    mark : str, optional
        type of Marker used for the plot, by default ''
    legend : bool, optional
        Display legend or not, by default True
    pic_save_name : str, optional
        name of the file where the figure is saved, by default 'transient.jpg'
    """  

    if len(y) > len(line_type):
        print('Invalid line_type list, we meed len(y) == len(line_type)')
        print('We will use default line type instead')
        line_type = []
        for counter, value in enumerate(y):
            line_type.append('-')

    plt.figure(num_fig)
    ax_JVs_plot = plt.axes()
    for i,line in zip(y,line_type):
        if plot_type == 1:
            ax_JVs_plot.semilogx(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 2:
            ax_JVs_plot.semilogy(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 3:
            ax_JVs_plot.loglog(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        else:
            ax_JVs_plot.plot(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)

    # legend
    if legend == True:
        plt.legend(loc='best',frameon=False,fontsize = 30)
    
    if xlimits != []:
        plt.xlim(xlimits)
    if ylimits != []:
        plt.ylim(ylimits)

    plt.grid(b=True,which='both')
    plt.xlabel('Time [s]')
    plt.ylabel('Current Density [mA cm$^{-2}$]')
    # plt.tight_layout()
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)

def zimt_tj_JV_plot(num_fig,data_tj,x='Vext',y=['Jext'],xlimits=[],ylimits=[],plot_type=0,labels='',colors='b',line_type = ['-'],mark='',legend=True,save_yes='False',pic_save_name='transient_JV.jpg'):
    """ Make tj_file transient current-voltage curve plot for ZimT  
    Default Voltage on the x axis 
    
    Parameters
    ----------
    num_fig : int
        number of the fig to plot tj
    data_tj : DataFrame
        Panda DataFrame containing tj_file
    x : str, optional
        xaxis data, by default 'Va'
    y : list of str, optional
        yaxis data can be multiple like ['Jext','Jncat'], by default ['Jext']
    xlimits : list, optional
        x axis limits if = [] it lets python chose limits, by default []
    ylimits : list, optional
        y axis limits if = [] it lets python chose limits, by default []
    plot_type : int, optional
        type of plot 1 = logx, 2 = logy, 3 = loglog else linlin (default = linlin), by default 0
    labels : str, optional
        label of the tj, by default ''
    colors : str, optional
        color for the JV line, by default 'b'
    line_type : list, optional
        type of line used for the plot
        size line_type need to be = size(y), by default ['-']
    mark : str, optional
        type of Marker used for the plot, by default ''
    legend : bool, optional
        Display legend or not, by default True
    save_yes : bool, optional
        If True, save JV as an image with the  file name defined by "pic_save_name", by default False
    pic_save_name : str, optional
        name of the file where the figure is saved, by default 'transient_JV.jpg'
    """  

    if len(y) > len(line_type):
        print('Invalid line_type list, we meed len(y) == len(line_type)')
        print('We will use default line type instead')
        line_type = []
        for counter, value in enumerate(y):
            line_type.append('-')

    plt.figure(num_fig)
    ax_JVs_plot = plt.axes()
    for i,line in zip(y,line_type):
        if plot_type == 1:
            ax_JVs_plot.semilogx(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 2:
            ax_JVs_plot.semilogy(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 3:
            ax_JVs_plot.loglog(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        else:
            ax_JVs_plot.plot(data_tj[x],data_tj[i]/10,color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)

    # legend
    if legend == True:
        plt.legend(loc='best',frameon=False,fontsize = 30)
    if xlimits != []:
        plt.xlim(xlimits)
    if ylimits != []:
        plt.ylim(ylimits)
    plt.grid(b=True,which='both')
    plt.xlabel('Voltage [V]')
    plt.ylabel('Current Density [mA cm$^{-2}$]')
    plt.tight_layout()
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)


def zimt_Voltage_transient_plot(num_fig,data_tj,x='t',y=['Vext'],xlimits=[],ylimits=[],plot_type=0,labels='',colors='b',line_type = ['-'],mark='',legend=True,save_yes=False,pic_save_name='transient_volt.jpg'):
    """ Make tj_file transientvoltage curve plot for ZimT  
    Default time on the x axis in $\mu$s
    
    Parameters
    ----------
    num_fig : int
        number of the fig to plot tj
    data_tj : DataFrame
        Panda DataFrame containing tj_file
    x : str, optional
        xaxis data, by default 't'
    y : list of str, optional
        yaxis data can be multiple like ['Vext','Va'], by default ['Vext']
    xlimits : list, optional
        x axis limits if = [] it lets python chose limits, by default []
    ylimits : list, optional
        y axis limits if = [] it lets python chose limits, by default []
    plot_type : int, optional
        type of plot 1 = logx, 2 = logy, 3 = loglog else linlin (default = linlin), by default 0
    labels : str, optional
        label of the tj, by default ''
    colors : str, optional
        color for the JV line, by default 'b'
    line_type : list, optional
        type of line used for the plot
        size line_type need to be = size(y), by default ['-']
    mark : str, optional
        type of Marker used for the plot, by default ''
    legend : bool, optional
        Display legend or not, by default True
    save_yes : bool, optional
        If True, save JV as an image with the  file name defined by "pic_save_name", by default False
    pic_save_name : str, optional
        name of the file where the figure is saved, by default 'transient_volt.jpg'
    """  

    if len(y) != len(line_type):
        print('Invalid line_type list, we meed len(y) == len(line_type)')
        print('We will use default line type instead')
        line_type = []
        for counter, value in enumerate(y):
            line_type.append('-')

    plt.figure(num_fig)
    ax_JVs_plot = plt.axes()
    for i,line in zip(y,line_type):
        if plot_type == 1:
            ax_JVs_plot.semilogx(data_tj[x],data_tj[i],color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 2:
            ax_JVs_plot.semilogy(data_tj[x],data_tj[i],color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        elif plot_type == 3:
            ax_JVs_plot.loglog(data_tj[x],data_tj[i],color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)
            
        else:
            ax_JVs_plot.plot(data_tj[x],data_tj[i],color=colors,label=labels,linestyle=line,marker=mark,markeredgecolor=colors,markersize=10,markerfacecolor='None',markeredgewidth = 3)

    # legend
    if legend == True:
        plt.legend(loc='best',frameon=False,fontsize = 30)
    if xlimits != []:
        plt.xlim(xlimits)
    if ylimits != []:
        plt.ylim(ylimits)
    plt.grid(b=True,which='both')
    plt.xlabel('Time [s]')
    plt.ylabel('Volatge [V]')
    plt.tight_layout()
    if save_yes:
        plt.savefig(pic_save_name,dpi=300,transparent=True)



######################################################################
#################### Useful for multiprocessing ######################
######################################################################


def starmap_with_kwargs(pool, fn, args_iter, kwargs_iter):
    args_for_starmap = zip(repeat(fn), args_iter, kwargs_iter)
    return pool.starmap(apply_args_and_kwargs, args_for_starmap)

def apply_args_and_kwargs(fn, args, kwargs):
    return fn(*args, **kwargs)

######################################################################
####################### Cleaning output files ########################
######################################################################

def clean_up_output(filename_start,path):
    """Delete output files from the simulation

    Parameters
    ----------
    filename_start : str
        string containing the begining of the filename to delete
    path : str
        path to the directory where we clean the output
    """ 
    for fname in os.listdir(path):
        if fname.startswith(filename_start):
            os.remove(os.path.join(path,fname))

def store_output_in_folder(filenames,folder_name,path):
    """Move output files from the simulation into new folder

    Parameters
    ----------
    filenames : list of str
        list of string containing the name of the files to move
    folder_name : str
        name of the folder where we store the output files
    path : str
        directory of the folder_name (creates one if it does not already exist)
    """    

    # Create directory if it does not exist
    if not os.path.exists(Path(path,folder_name)):
        os.makedirs(Path(path,folder_name))
    # move file into the new folder
    for f in filenames:
        if os.path.exists(Path(path,f)):
            os.replace(Path(path,f),Path(path,folder_name,f))


def calcZC(tV, tJ, freq, equalCircuit):

    tV_fft = fftpack.fft(tV)[1:int(len(tV)/2)]

    tJ_fft = fftpack.fft(tJ)[1:int(len(tV)/2)]
    Z = max(tV_fft)/max(tJ_fft)
   

    if equalCircuit == 'RCp':
        C = 1/(2*np.pi*freq) * (1/Z).imag 
    elif equalCircuit == 'RCs' :
        C = -1/(2*np.pi*freq*Z.imag) 
  

    # plt.stem(abs(tV_fft))
    # plt.show()

    # plt.stem(abs(tJ_fft))
    # plt.show()

    return Z, C


