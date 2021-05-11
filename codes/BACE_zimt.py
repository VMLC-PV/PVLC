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
from scipy.optimize import curve_fit
from scipy import constants
from time import time
from itertools import repeat
# Don't show warnings
warnings.filterwarnings("ignore")
# Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_BACE

# Main Program
def main():    
    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                         # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
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
    path2ZimT = 'Simulation_program/DDSuite_v400/ZimT'+slash    # Path to ZimT in curr_dir

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    run_simu = True                                         # Rerun simu?
    plot_tjs = True                                             # make plot ?
    plot_output = True
    move_ouput_2_folder = True
    Store_folder = 'BACE'+slash
    clean_output = False
    L = 100e-9                                                  # Device thickness (m)
    Gens = [1.44e28]                                               # Max generation rate for the gaussian laser pulse
    Vpres = [0.77]
    Vextrs =[0,-1,-2,-3,-4]

    # Initialize 
    sys_lst,str_lst,nextr,path_lst,tj_lst,tVG_lst = [],[],[],[],[],[]
    idx = 0
    start = time()
    p = multiprocessing.Pool(max_jobs)

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    colors = cm.viridis((np.linspace(0,1,max(len(Vextrs),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

   
        
    # Generate tVG files and str_lst for Dark simulation
    for Vpre in Vpres:
        for Vextr in Vextrs:
            zimt_BACE(1e-8,1e-6,0e28,Vpre,Vextr,1e-10,time_exp=True,steps=100,tVG_name=curr_dir+slash+path2ZimT+'tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Vpre,Vextr))
            str_lst.append('-FailureMode 1 -L '+str(L)+' -tVG_file tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt -tj_file tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr,Vpre,Vextr))
            sys_lst.append(system)
            path_lst.append(curr_dir+slash+path2ZimT)
            tVG_lst.append('tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Vpre,Vextr))
            tj_lst.append('tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr))
        

    # Generate tVG files and str_lst for light pulse simulation
    for Gen in Gens:
        for Vpre in Vpres:
            for Vextr in Vextrs:
                zimt_BACE(1e-8,1e-6,Gen,Vpre,Vextr,1e-10,time_exp=True,steps=100,tVG_name=curr_dir+slash+path2ZimT+'tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Gen,Vpre,Vextr))
                str_lst.append('-FailureMode 1 -L '+str(L)+' -tVG_file tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt -tj_file tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr,Gen,Vpre,Vextr))
                sys_lst.append(system)
                path_lst.append(curr_dir+slash+path2ZimT)
                tVG_lst.append('tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Gen,Vpre,Vextr))
                tj_lst.append('tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr))
        

    if run_simu:       
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
    ################## Plotting ############################
    ########################################################
    if plot_tjs:
        for Gen in Gens:
            for Vpre in Vpres:
                nextr_dumb = []
                idx = 0
                for Vextr in Vextrs: 
                    data_tjdark = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr),delim_whitespace=True)
                    data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr),delim_whitespace=True)
                    data_tj['JBACE'] = abs(data_tj['Jext']-data_tjdark['Jext'])
                    zimt_tj_plot(num_fig_tjs,data_tj,y=['JBACE'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vpre {:.2f} Vext {:.2f}'.format(Vpre,Vextr),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = curr_dir+slash+path2ZimT+Store_folder+'transient.jpg') 
                    nextr.append(abs(simps(data_tj['JBACE']/(q*L),x=data_tj['t'])))
                    idx = idx + 1
    print(nextr)         

    
    if plot_output:
            # Make plots
            f_Q = plt.figure(10,figsize=size_fig)
        
            idx = 0
            for i in range(len(Gens)):
                for j in range(len(Vpres)):
                    plt.plot(Vextrs,nextr,'o',markersize=10)
                    idx = idx + 1

            plt.grid(b=True,which='both')
            plt.xlabel('Extraction voltage [V]')
            plt.ylabel('Carrier Density [m$^{-3}$]')
            plt.tight_layout()
            plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'extrac_charges.jpg',dpi=600,transparent=True)
    
    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
if __name__ == '__main__':
    main()
    
    
    

