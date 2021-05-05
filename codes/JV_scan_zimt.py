###################################################
###### Simulate JV hysteresis using ZimT ##########
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
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,zimt_tj_JV_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_JV_sweep

# Main Program
def JV_scan(str2run=[''],Gens=[0],num_fig=0,path2ZimT='',Store_folder=''):
    """Run simple (one direction) JV scan

    Parameters
    ----------
    str2run : 1-D sequence of strings
        Input string to run with ZimT,  by default ['']
    
    Gens : 1-D sequence of floats, optional
        Max generation rate (unit: m^-3 s^-1), by default [0]

    num_fig : int, optional
        Starting figure number for the first plot, for the other figure the figure number will be incremented as num_fig + 1, by default 0
        
    path2ZimT : str, optional
        Path to folder containing ./ZimT or zimt.exe in current directory, by default = ''

    Store_folder : str, optional
        Path to folder where the output of the simulation is stored in path2ZimT directory, by default = ''
    """     

    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = 1                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
            max_jobs = 1
            slash = '/'
            try:
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass
    else:
        slash = '/'

    curr_dir = os.getcwd()                                      # Current working directory
    path2ZimT = path2ZimT+slash                                 # Path to ZimT in curr_dir
    Store_folder = Store_folder+slash                           # Path to folder containing output data in curr_dir
    run_simu = False                                            # Rerun simu?  
    move_ouput_2_folder = True                                  # Move (True) output of simulation to Store_folder
    clean_output = False                                        # Clean output after simulation (delete files)

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    Vmin = -0.1                                                 # Minimum voltage (unit: V)
    Vmax = 1.1                                                  # Maximum voltage (unit: V)
    scans = [0.1]                                               # Scan speed (unit: V/s)
    steps = 100                                                 # Number of voltage step

    # Figures control
    size_fig = (16, 12)
    num_fig_JV= 0

    # plot JV cyrve
    plot_JV = True                                        # make plot ?
    if plot_JV:
        num_fig_JV = num_fig
        f_JVs = plt.figure(num_fig_JV,figsize=size_fig) 

    colors = cm.viridis((np.linspace(0,1,max(len(str2run),4)+1)) ) # Color range for plotting
   

    # Initialize 
    str_lst,sys_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    
    if run_simu:
    # Generate tVG files and str_lst
        for strs in str2run:
            for scan in scans:
                for Gen in Gens:
                    zimt_JV_sweep(Vmin,Vmax,scan,Gen,steps,tVG_name=curr_dir+slash+path2ZimT+'tVG_JV_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                    str_lst.append(strs+' -tVG_file tVG_JV_scan_{:.2e}_G_{:.2e}.txt -tj_file tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen,scan,Gen))
                    sys_lst.append(system)
                    path_lst.append(curr_dir+slash+path2ZimT)
                    tVG_lst.append('tVG_JV_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                    tj_lst.append('tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen))

        print(str_lst)
        # Run ZimT
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
    ################## Plotting ############################
    ########################################################
    if plot_JV:
        for strs in str2run:
            idx = 0
            for scan in scans:
                for Gen in Gens:
                        data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_JV_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen),delim_whitespace=True)
                        zimt_tj_JV_plot(num_fig_JV,data_tj,x='Vdev',y=['Jdev'],xlimits = [Vmin-0.1,Vmax+0.1],ylimits=[-22,2] ,colors=colors[idx],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=curr_dir+slash+path2ZimT+Store_folder+'JV_transient.jpg')
                        idx = idx + 1          

    
    # Make plots
    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)

    

    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':

    JV_scan(str2run=['-L 120e-9'],Gens=[0,1e28],num_fig=0,path2ZimT = 'Simulation_program/ZimT043_BETA',Store_folder='test_output_JV_scan')

    
    
    

