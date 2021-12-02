###################################################
###### Simulate JV sweep using ZimT ##########
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants
from time import time
# Import homemade package by VLC
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *
from VLC_units.cleanup_folder.clean_folder import *
from VLC_units.useful_functions.aux_func import *

# Main Program
def JV_sweep(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = True, plot_tjs = True, move_ouput_2_folder = True, Store_folder = 'JV_sweep',clean_output = False,verbose = True):  
    """Run single JV sweep simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_JV_sweepfunction (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'JV_sweep'
    clean_output : bool, optional
        Clean up output?, by default False
    verbose : bool, optional
        Verbose?, by default True
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
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass

    curr_dir = os.getcwd()              # Current working directory
    if path2ZimT is None:
        path2ZimT = os.path.join(os.getcwd(),'Simulation_program/SIMsalabim_v425/ZimT')                  # Path to ZimT in curr_dir

 
    ## Physics constants
    q = constants.value(u'elementary charge')

    # TPC input
    # see zimt_TPC in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No JV_sweep input dictionary given, using default values')
        Vmin = -0.1                                                 # Minimum voltage (unit: V)
        Vmax = 1.1                                                  # Maximum voltage (unit: V)
        scans = [0.1]                                               # Scan speed (unit: V/s)
        Gens = [1.28e28]                                             # Average generation rate (unit: m^-3 s^-1) 
        steps = 100                                                 # Number of voltage step
    else:
        if 'Vmin' in input_dic:
            Vmin = input_dic['Vmin']
        if 'Vmax' in input_dic:
            Vmax = input_dic['Vmax']
        if 'scans' in input_dic:
            scans = input_dic['scans']
        if 'Gens' in input_dic:
            Gens = input_dic['Gens']
        if 'steps' in input_dic:
            steps = input_dic['steps']

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = '-accDens 0.05 '  # add any fixed string to the simulation command

    # Initialize 
    str_lst,code_name_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    size_fig = (16, 12)

    # plot JV curve
    if plot_tjs:
        num_fig_JV = 0
        f_JVs = plt.figure(num_fig_JV,figsize=size_fig) 

    colors = cm.viridis((np.linspace(0,1,max(len(Gens),4)+1)) ) # Color range for plotting
       
    if run_simu:
    # Generate tVG files and str_lst
        for scan in scans:
            for Gen in Gens:
                zimt_JV_sweep(Vmin,Vmax,scan,Gen,steps,tVG_name=os.path.join(path2ZimT,'tVG_JV_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen)))
                str_lst.append(fixed_str+' -tVG_file tVG_JV_scan_{:.2e}_G_{:.2e}.txt -tj_file tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen,scan,Gen))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_JV_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                tj_lst.append('tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen))

        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose=verbose)
            
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    
    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    
    ########################################################
    ################## Plotting ############################
    ########################################################
    if plot_tjs:
        # Plot JV curve
        idx = 0
        for scan in scans:
            for Gen in Gens:
                    data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen)),delim_whitespace=True)
                    zimt_tj_JV_plot(num_fig_JV,data_tj,x='Vext',y=['Jext'],xlimits = [Vmin-0.1,Vmax+0.1],ylimits=[-22,2] ,colors=colors[idx],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=os.path.join(path2ZimT,Store_folder,'JV_transient.jpg'))
                    idx += 1          

    
    # Make plots
    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        if verbose:
            print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':

    JV_sweep()
    
    
    

