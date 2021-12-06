###########################################################
########## Simulate Voltage step using ZimT ###############
###########################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from time import time
from scipy import constants,interpolate
from scipy.integrate import simps
# Import homemade package by VLC
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *
from VLC_units.cleanup_folder.clean_folder import *
from VLC_units.useful_functions.aux_func import *

# Main Program
def Voltage_step_zimt(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = True, plot_tjs = True, plot_Charge = True, move_ouput_2_folder = True, Store_folder = 'Voltage_step',clean_output = False,verbose = True):
    """Run voltage step simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_voltage_step function (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    plot_Charge : bool, optional
        calculate the extracted charge (integrating the current) and make extracted charge plot?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'Voltage_step'
    clean_output : bool, optional
        Clean up output?, by default False
    verbose : bool, optional
        Verbose?, by default True
    """  
    ## General Inputs
    warnings.filterwarnings("ignore")           # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                 # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    do_multiprocessing = False                      # Use multiprocessing
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
    
    # Physics constants
    q = constants.value(u'elementary charge')

    ## Voltage step Inputs
    # see zimt_voltage_step in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No voltage step input dictionary given, using default values')
        tmin = 1e-8                                                  # First time step after 0 (unit: s)
        tmax = 5e-5                                                  # Final time step (unit: s)
        Gens = [1.28e28]                                                # Constant generation rate for the entire duration of the pulse (unit m^-3 s^-1)
        Vstarts = [0.95]                                              # Initial voltage (unit: V)
        Vfinals =[0]                                                 # Final voltage (unit: V)
        steps = 100                                                  # Number of time step
    else :
        if 'tmin' not in input_dic.keys():
            tmin = input_dic['tmin']
        if 'tmax' not in input_dic.keys():
            tmax = input_dic['tmax']
        if 'Gens' not in input_dic.keys():
            Gens = input_dic['Gens']
        if 'Vstarts' not in input_dic.keys():
            Vstarts = input_dic['Vstarts']
        if 'Vfinals' not in input_dic.keys():
            Vfinals = input_dic['Vfinals']
        if 'steps' not in input_dic.keys():
            steps = input_dic['steps']

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ''  # add any fixed string to the simulation command

    # Get thicknesses (needed later)
    ParFileDic = ReadParameterFile(f"{path2ZimT}/device_parameters.txt") # Read device parameters
    L = float(ChosePar('L',GetParFromStr(fixed_str),ParFileDic))
    L_LTL = float(ChosePar('L_LTL',GetParFromStr(fixed_str),ParFileDic))
    L_RTL = float(ChosePar('L_RTL',GetParFromStr(fixed_str),ParFileDic))

    # Initialize 
    code_name_lst,str_lst,charge,path_lst,tj_lst,tVG_lst = [],[],[],[],[],[]
    idx = 0
    start = time()
    num_fig = 0

    # Figures control
    size_fig = (16, 12)
    colors = cm.viridis((np.linspace(0,1,max(len(Vfinals),4)+1)) ) # Color range for plotting
    save_fig = True

    # plot current transient
    if plot_tjs:
        num_fig_Curr = num_fig
        num_fig = num_fig + 1
        f_Curr = plt.figure(num_fig_Curr,figsize=size_fig)
    
    # plot extracted charges (integrate the current)
    if plot_tjs and plot_Charge:
        num_fig_Charge = num_fig
        num_fig = num_fig + 1
        f_Charge = plt.figure(num_fig_Charge,figsize=size_fig)

   

    if run_simu:
        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for Vstart in Vstarts:
                for Vfinal in Vfinals:
                    zimt_voltage_step(tmin,tmax,Vstart,Vfinal,Gen,steps,trf = 10e-9,time_exp =True,tVG_name=os.path.join(path2ZimT,'tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt'.format(Gen,Vstart,Vfinal)))
                    str_lst.append(fixed_str+' -tVG_file tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt -tj_file tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal,Gen,Vstart,Vfinal))
                    code_name_lst.append('zimt')
                    path_lst.append(path2ZimT)
                    tVG_lst.append('tVG_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.txt'.format(Gen,Vstart,Vfinal))
                    tj_lst.append('tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal))       
        
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
        for Gen in Gens:
            for Vstart in Vstarts:
                charge_dumb = []
                idx = 0
                for Vfinal in Vfinals:
                    # Load Data 
                    data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_Vstep_G_{:.1e}_Vstart_{:.2f}_Vfinal{:.2f}.dat'.format(Gen,Vstart,Vfinal)),delim_whitespace=True)

                    # Make Plot
                    zimt_tj_plot(num_fig_Curr,data_tj,y=['Jext'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vstart {:.2f} Vfinal {:.2f}'.format(Vstart,Vfinal),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = os.path.join(path2ZimT,Store_folder,'Current_transient.jpg')) 

                    # Calculate number of extracted charges
                    if plot_Charge:
                        extrat_charge = integrate.cumtrapz(data_tj['Jext'], data_tj['t'], initial=0) 
                        data_tj['Qext']= extrat_charge/(q*(L-L_LTL-L_RTL))
                        zimt_tj_plot(num_fig_Charge,data_tj,y=['Qext'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vstart {:.2f} Vfinal {:.2f}'.format(Vstart,Vfinal),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = os.path.join(path2ZimT,Store_folder,'Charge_transient.jpg')) 
                    idx += 1
             

    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':
    Voltage_step_zimt()
    plt.show()
    
    
    

