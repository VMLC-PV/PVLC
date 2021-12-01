###################################################
########## Simulate CELIV using ZimT ###############
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

# Main Program
def CELIV(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = True, plot_tjs = True, move_ouput_2_folder = True, Store_folder = 'CELIV',clean_output = False,verbose = True):
    """Run CELIV simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_CELIV function (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'CELIV'
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

    ## CELIV input
    # see zimt_CELIV in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No CELIV input dictionary given, using default values')
        Gens = [0,1e22]                 # Number of generated charge from the gaussian laser pulse 
        slopes =[-1/1e-6]               # Slope of the voltage ramp [V/s]
        Voffsets =[0]                   # Voltage offset [V]
        tdelay = 0                      # Delay between the laser pulse and the begining of the ramp [s]
        tpulse = 5e-8                   # Middle of the gaussian laser pulse [s]
        tmin = 1e-8                     # First time step after 0 [s]
        tmax = 3e-6                     # Final time step [s]
        tstep = 1e-9                    # Time step for the linear regime before the ramp [s]
        steps = 100                     # Number of steps
    else:
        if 'Gens' in input_dic.keys():
            Gens = input_dic['Gens']
        if 'slopes' in input_dic.keys():
            slopes = input_dic['slopes']
        if 'Voffsets' in input_dic.keys():
            Voffsets = input_dic['Voffsets']
        if 'tdelay' in input_dic.keys():    
            tdelay = input_dic['tdelay']
        if 'tpulse' in input_dic.keys():
            tpulse = input_dic['tpulse']
        if 'tmin' in input_dic.keys():
            tmin = input_dic['tmin']
        if 'tmax' in input_dic.keys():
            tmax = input_dic['tmax']
        if 'tstep' in input_dic.keys():
            tstep = input_dic['tstep']
        if 'steps' in input_dic.keys():
            steps = input_dic['steps'] 
    
    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ''  # add any fixed string to the simulation command

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs = 0
    colors = cm.viridis((np.linspace(0,1,max(len(slopes),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

    # Initialize 
    code_name_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    if run_simu: # Run simulation  
        # Generate tVG files and str_lst
        for Gen in Gens:
            for slope in slopes:
                for Voffset in Voffsets:
                    zimt_CELIV(tmin,tmax,Voffset,slope,Gen,tpulse,tstep,tdelay,width_pulse = 6e-9,time_exp=True,steps=steps,tVG_name=os.path.join(path2ZimT,'tVG_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.txt'.format(Gen,slope,Voffset)))
                    str_lst.append(fixed_str+' -tVG_file tVG_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.txt -tj_file tj_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.dat'.format(Gen,slope,Voffset,Gen,slope,Voffset))
                    code_name_lst.append('zimt')
                    path_lst.append(path2ZimT)
                    tVG_lst.append('tVG_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.txt'.format(Gen,slope,Voffset))
                    tj_lst.append('tj_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.dat'.format(Gen,slope,Voffset))
        
   
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
    ################## JVs_file ############################
    ########################################################
    if plot_tjs:
        for Gen in Gens: # Loop over Gens
            idx = 0
            for slope in slopes: # Loop over slopes
                for Voffset in Voffsets: # Loop over Voffsets
                    data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_CELIV_G_{:.2e}_slope_{:.2f}_Voff_{:.2f}.dat'.format(Gen,slope,Voffset)),delim_whitespace=True)
                    zimt_tj_plot(num_fig_tjs,data_tj,colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = os.path.join(path2ZimT,Store_folder,'transient.jpg'))
                idx += 1 

    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))

    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':
    CELIV()
    
    
    

