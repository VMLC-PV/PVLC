###################################################
########## Simulate TPV using ZimT ################
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from scipy.integrate import simps
from scipy.optimize import curve_fit
from scipy import constants
from time import time
from itertools import repeat
# Don't show warnings
warnings.filterwarnings("ignore")
# Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,zimt_Voltage_transient_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_TPV

# Main Program
def TPV():    
    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
            max_jobs = 1
            slash = '/'
            try:
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass
    else:
        slash = '/'

    curr_dir = os.getcwd()              # Current working directory
    path2ZimT = 'Simulation_program/DDSuite_v400/ZimT'+slash                      # Path to ZimT in curr_dir

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    run_simu = True                                        # Rerun simu?
    plot_tjs = True                                        # make plot ?
    plot_output = False
    move_ouput_2_folder = True
    Store_folder = 'TPV'+slash
    clean_output = False
    L = 140e-9                                                  # Device thickness (m)
    Gens = [1e27]                                               # Max generation rate for the gaussian laser pulse
    G0s = [1e26]#,1e25]
    tpulse = 5e-8

    # Initialize 
    sys_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    colors = cm.viridis((np.linspace(0,1,max(len(Gens),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

    if run_simu:           
    
        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for G0 in G0s:
                zimt_TPV(1e-8,1,Gen,G0,1e-8,tpulse,time_exp =True,tVG_name=curr_dir+slash+path2ZimT+'tVG_TPV_G_{:.1e}_G0_{:.2f}.txt'.format(Gen,G0))
                str_lst.append('-L '+str(L)+' -tVG_file tVG_TPV_G_{:.1e}_G0_{:.2f}.txt -tj_file tj_TPV_G_{:.1e}_G0_{:.2f}.dat'.format(Gen,G0,Gen,G0))
                sys_lst.append(system)
                path_lst.append(curr_dir+slash+path2ZimT)
                tVG_lst.append('tVG_TPV_G_{:.1e}_G0_{:.2f}.txt'.format(Gen,G0))
                tj_lst.append('tj_TPV_G_{:.1e}_G0_{:.2f}.dat'.format(Gen,G0))
        
 
        # Run ZimT
        # str_lst = str_lst[::-1] # reverse list order to start with longest delays
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_zimt,list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)
        p.close()
        p.join()
         
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,curr_dir+slash+path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,curr_dir+slash+path2ZimT)

    ########################################################
    ################## JVs_file ############################
    ########################################################
    if plot_tjs:
        for G0 in G0s:
            idx = 0
            for Gen in Gens:
                data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_TPV_G_{:.1e}_G0_{:.2f}.dat'.format(Gen,G0),delim_whitespace=True)
                zimt_Voltage_transient_plot(num_fig_tjs,data_tj,y=['Vext'],xlimits=[-1e-4,1e-3],colors=colors[idx],plot_type=0,save_yes=True,pic_save_name = curr_dir+slash+path2ZimT+Store_folder+'transient.jpg') 
                idx = idx + 1

        
    # if plot_output:
            # Make plots
    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)

    

    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':
    TPV()