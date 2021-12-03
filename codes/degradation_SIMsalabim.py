#########################################################
########## Run and plot SIMsalabim output ###############
#########################################################
# by Vincent M. Le Corre
# Package import
import os,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# Import homemade package by VLC
from VLC_units.plots.SimSS_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.useful_functions.aux_func import *



def run_plot_SIMsalabim(fixed_str = None, parameters = None, path2SIMsalabim = None, run_simu = True, plot_JVs = True, plot_effs = True ,  verbose = True):
    """Run and plot SIMsalabim outputs.

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    parameters : dict, optional
        Dictionary with the input for the SIMsalabim function, we will simulate all permutations, by default None.
    path2SIMsalabim : str, optional
        Path to Simss, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_JVs : bool, optional
        Make JV plot?, by default True
    plot_effs : bool, optional
        Make efficiency versus time plot, by default True
    verbose : bool, optional
        Verbose?, by default True
    
    Returns
    -------
    JV_files,Var_files,scPars_files : list,list
        Lists of JV and variable and performance filenames
    """
    ## General Inputs
    warnings.filterwarnings("ignore")           # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                 # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    do_multiprocessing = True                      # Use multiprocessing
    if system == 'Windows':                     # cannot easily do multiprocessing in Windows
            max_jobs = 1
            do_multiprocessing = False
            try:                                # kill all running jobs to avoid conflicts
                os.system('taskkill.exe /F /IM simss.exe')
            except:
                pass


    curr_dir = os.getcwd()              # Current working directory
    if path2SIMsalabim  is None:
        path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabim_v425/SimSS')                    # Path to SimSS
                                                                                            # Rerun simulation

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ''  # add any fixed string to the simulation command

    # Parameters to vary
    if parameters is None:
        time = np.geomspace(1e-1,1000,13)
        time = np.insert(time,0,0)
        a = MonoExpDecay(time, 1, 5e-8, 1e-10)
        b = StretchedExp(time, 1, 3/12, 4.5e-8, 3e-8)
        mob = a + b
        lang = MonoExpInc(time, 10, 0.001, 0.002,)

        parameters = []
        parameters.append({'name':'mun_0','values':list(mob)})
        parameters.append({'name':'mup_0','values':list(mob)})
        parameters.append({'name':'-kdirect','values':list(lang)})
    else:
        for i in range(len(parameters)): # get time from parameters and then delete the line
            if parameters[i]['name'] == 'time':
                time = np.array(parameters[i]['values'])
                del parameters[i]

    ParFileDic = ReadParameterFile(f"{path2SIMsalabim}/device_parameters.txt") # Read device parameters

    # Initialize     
    str_lst,labels,JV_files,Var_files,scPars_files,code_name_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (10, 8)
    num_fig = 0
    # JV_file plots
    if plot_JVs:
        num_fig += 1
        num_JV_plot = num_fig
        f_JVs = plt.figure(num_JV_plot,figsize=size_fig)
    
    plot_effs = True # Make efficiency plot
    if plot_effs:
        num_fig = num_fig + 1
        num_eff_plot = num_fig
        f_effs = plt.figure(num_eff_plot,figsize=size_fig)
    
    for param in parameters: # initalize lists of values and names
        val.append(param['values'])
        nam.append(param['name'])

    
    idx = 0
    
    for idx, _ in enumerate(parameters[0]['values']): # loop over values
    # for i in list(itertools.product(*val)): # iterate over all combinations of parameters
        str_line = ''
        lab = ''
        JV_name = 'JV'
        Var_name = 'Var'
        scPars_name = 'scPars'
        for param in parameters: # loop over parameters
            str_line = str_line +'-'+param['name']+' {:.2e} '.format(param['values'][idx])
            lab = lab+param['name']+' {:.2e} '.format(param['values'][idx])
            JV_name = JV_name +'_'+param['name']+'_{:.2e}'.format(param['values'][idx])
            Var_name = Var_name +'_'+param['name']+'_{:.2e}'.format(param['values'][idx])
            scPars_name = scPars_name +'_'+param['name']+'_{:.2e}'.format(param['values'][idx])
        str_lst.append(fixed_str+ ' ' +str_line+ '-JV_file '+JV_name+ '.dat -Var_file '+Var_name+'.dat -scPars_file '+scPars_name+'.dat')# -ExpJV '+JVexp_lst[idx])
        JV_files.append(os.path.join(path2SIMsalabim , str(JV_name+ '.dat')))
        Var_files.append(os.path.join(path2SIMsalabim , str(Var_name+ '.dat')))

        # Check if we output scPars file, i.e. if we have a solar cell
        Gehp = float(ChosePar('Gehp',GetParFromStr(fixed_str+ ' ' +str_line),ParFileDic))
        Gfrac = float(ChosePar('Gfrac',GetParFromStr(fixed_str+ ' ' +str_line),ParFileDic))
        if Gfrac * Gehp > 0:
            scPars_files.append(os.path.join(path2SIMsalabim , str(scPars_name+ '.dat')))
        labels.append(lab)
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)
        idx += 1

    # print(str_lst)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots

    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)



    idx = 0
    Voc,Jsc,FF,PCE = [],[],[],[]
    for Simu_str,lbl,JV_file_name,Var_file_name,scPars_file_name in zip(str_lst,labels,JV_files,Var_files,scPars_files):

        # Get performance data
        data_scpars = pd.read_csv(scPars_file_name,delim_whitespace=True) # Load scPars_file
        Jsc.append(float(data_scpars['Jsc']))
        Voc.append(float(data_scpars['Voc']))
        FF.append(float(data_scpars['FF']))
        PCE.append(float((abs(data_scpars['Jsc'])/10)*data_scpars['Voc']*data_scpars['FF']))

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file

            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,xlimits=[-0.5,1.1],ylimits=[-25,5],pic_save_name=os.path.join(path2SIMsalabim , str('JV'+ext_save_pic)))
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=2,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp)
    
        idx += 1
    
    # Efficiency plots
    if plot_effs:
        plt.figure(num_eff_plot)
        x_var = time
        x_label = 'Time [s]'
        Jsc,Voc,FF,PCE = np.asarray(Jsc),np.asarray(Voc),np.asarray(FF),np.asarray(PCE)

        plt.subplot(221)
        plt.semilogx(x_var,PCE/max(PCE),linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('Norm. PCE')
        plt.ylim([0.2,1.2])
        plt.grid(b=True,which='both')
        plt.subplot(222)
        plt.semilogx(x_var,FF,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('FF')
        # plt.ylim([0,1])
        plt.grid(b=True,which='both')
        plt.subplot(223)
        plt.semilogx(x_var,Voc,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('V$_{OC}$ [V]')
        # plt.ylim([0.7,1.3])
        plt.grid(b=True,which='both')
        plt.subplot(224)
        plt.semilogx(x_var,abs(np.asarray(Jsc)),linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('J$_{SC}$ mA cm$^{-2}$')
        # plt.ylim([10,30])
        # plt.legend(loc='best', title="Scan direction:",fontsize=35)
        plt.tight_layout()
        plt.grid(b=True,which='both')

    return JV_files,Var_files,scPars_files
  


if __name__ == '__main__':
    
    run_plot_SIMsalabim()
    plt.show()