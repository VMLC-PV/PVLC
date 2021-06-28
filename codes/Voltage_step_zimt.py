###################################################
########## Simulate BACE using ZimT ###############
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
from scipy import integrate
from scipy.optimize import curve_fit
from scipy import constants
from time import time
from itertools import repeat
# Don't show warnings
warnings.filterwarnings("ignore")
# Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_voltage_step

# Main Program
def Voltage_step(L,num_fig=0,str2run='',path2ZimT='',Store_folder=''): 
    """Run voltage step simulation

    Parameters
    ----------
    L : float
        Device total thickness (unit: m)
    
    num_fig : int, optional
        Starting figure number for the first plot, for the other figure the figure number will be incremented as num_fig + 1, by default 0

    str2run : str, optional
        Input string to run with ZimT,  by default ''

    path2ZimT : str, optional
        Path to folder containing ./ZimT or zimt.exe in current directory, by default = ''

    Store_folder : str, optional
        Path to folder where the output of the simulation is stored in path2ZimT directory, by default = ''

    """       
    ## General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = 10                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
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
    run_simu = True                                            # Rerun simu?  
    move_ouput_2_folder = True                                  # Move (True) output of simulation to Store_folder
    clean_output = False                                        # Clean output after simulation (delete files)
    make_fit = False                                            # Make fit dn/dt

    ## Physics constants
    q = constants.value(u'elementary charge')
    eps_0 = constants.value(u'electric constant')

    ## Simulation input for zimt_voltage_step (see tVG_gen.py)   
    tmin = 1e-8                                                  # First time step after 0 (unit: s)
    tmax = 5e-6                                                     # Final time step (unit: s)
    Gens = [0e28]                                                # Constant generation rate for the entire duration of the pulse (unit m^-3 s^-1)
    Vstarts = [1.3]                                             # Initial voltage (unit: V)
    Vfinals =[0]                                                 # Final voltage (unit: V)
    steps = 100                                                  # Number of time step

    

    # Figures control
    size_fig = (16, 12)
    colors = cm.viridis((np.linspace(0,1,max(len(Vfinals),4)+1)) ) # Color range for plotting
    save_fig = True

    # plot current transient
    plot_Curr = True
    if plot_Curr:
        num_fig_Curr = num_fig
        num_fig = num_fig + 1
        f_Curr = plt.figure(num_fig_Curr,figsize=size_fig)
    
    # plot extracted charges (integrate the current)
    plot_Charge = True
    if plot_Curr and plot_Curr:
        num_fig_Charge = num_fig
        num_fig = num_fig + 1
        f_Charge = plt.figure(num_fig_Charge,figsize=size_fig)

   # Initialize 
    sys_lst,str_lst,charge,path_lst,tj_lst,tVG_lst = [],[],[],[],[],[]
    idx = 0
    start = time()

    if run_simu:
        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for Vstart in Vstarts:
                for Vfinal in Vfinals:
                    zimt_voltage_step(tmin,tmax,Vstart,Vfinal,Gen,steps,trf = 10e-9,time_exp =False,tVG_name=curr_dir+slash+path2ZimT+'tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt'.format(Gen,Vstart,Vfinal))
                    str_lst.append('-L '+str(L)+' -tVG_file tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt -tj_file tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal,Gen,Vstart,Vfinal))
                    sys_lst.append(system)
                    path_lst.append(curr_dir+slash+path2ZimT)
                    tVG_lst.append('tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt'.format(Gen,Vstart,Vfinal))
                    tj_lst.append('tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal))       
        
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
    if plot_Curr:
        for Gen in Gens:
            for Vstart in Vstarts:
                charge_dumb = []
                idx = 0
                for Vfinal in Vfinals:
                    # Load Data 
                    data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal),delim_whitespace=True)

                    # Make Plot
                    zimt_tj_plot(num_fig_Curr,data_tj,y=['Jdev'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vstart {:.2f} Vfinal {:.2f}'.format(Vstart,Vfinal),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = curr_dir+slash+path2ZimT+Store_folder+'Current_transient.jpg') 

                    # Calculate number of extracted charges
                    if plot_Charge:
                        extrat_charge = integrate.cumtrapz(data_tj['Jdev'], data_tj['t'], initial=0) 
                        data_tj['Qext']= extrat_charge/(q*L)
                        zimt_tj_plot(num_fig_Charge,data_tj,y=['Qext'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vstart {:.2f} Vfinal {:.2f}'.format(Vstart,Vfinal),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = curr_dir+slash+path2ZimT+Store_folder+'Charge_transient.jpg') 


                    idx = idx + 1
             


    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':
    Voltage_step(L=340e-9,num_fig=0,path2ZimT = 'Simulation_program/ZimT',Store_folder='test_output_Vstep')
    
    
    

