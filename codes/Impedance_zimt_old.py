###################################################
########## Simulate Impedance using ZimT ##########
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,warnings,cmath
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
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,Store_output_in_folder,clean_up_output,preprocess_Impedance_data,get_complex_impedance,fit_sin_func,sin_func
from tVG_gen import zimt_impedance

# Main Program
def main(): 
    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = 100                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
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
    path2ZimT = 'ZimT043_BETA'+slash                      # Path to ZimT in curr_dir

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    run_simu = True                                        # Rerun simu?
    plot_tjs = False                                        # make plot ?
    plot_output = True
    move_ouput_2_folder = True
    Store_folder = 'Impedance'+slash
    clean_output = False                                         # Make output plots?
    make_fit = False                                            # make fit dn/dt
    L = 140e-9                                                  # Device thickness (m)
    freqs = np.geomspace(1e-4,1e7,num=50,endpoint=True)                         # frequencies to simulate (Hz)
    Vapps = [0]                                                 # Applied voltage (V)
    Vamp = 0.01                                                 # Amplitude voltage perturbation
    Gen = 0                                                     # Average generation rate 
    # Initialize 
    str_lst,sys_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()
    
    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    colors = cm.viridis((np.linspace(0,1,max(len(freqs),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)
    f_impe = plt.figure(num_fig_tjs+1,figsize=size_fig)

    
  
    if run_simu:
        
        # Generate tVG files 
        for freq in freqs:
            for Vapp in Vapps:
                zimt_impedance(Vapp,Vamp,freq,Gen,steps=100,tVG_name=curr_dir+slash+path2ZimT+'tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq))
                str_lst.append('-L '+str(L)+' -tVG_file tVG_{:.2f}V_f_{:.1e}Hz.txt -tj_file tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq,Vapp,freq))
                sys_lst.append(system)
                path_lst.append(curr_dir+slash+path2ZimT)
                tVG_lst.append('tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq))
                tj_lst.append('tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq))
        
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

    for freq in freqs:
        for Vapp in Vapps:
            data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq),delim_whitespace=True)
            data_tj['Jdev_norm'] = data_tj['Jdev']/max(data_tj['Jdev'])*10
            data_tj['Vdev_norm'] = data_tj['Vdev']/max(data_tj['Vdev'])*10
            data_tj['t'] = data_tj['t']/max(data_tj['t'])*10
            if plot_tjs:
                zimt_tj_plot(num_fig_tjs,data_tj,y=['Jdev_norm'],ylimits= [-1.1,1.1],labels='Norm. J ' +sci_notation(freq,sig_fig=0)+' Hz',colors=colors[idx],plot_type=0,save_yes=True,legend=False)
                zimt_tj_plot(num_fig_tjs,data_tj,y=['Vdev_norm'],ylimits= [-1.1,1.1],labels='Norm. V',line_type=['--'],colors=colors[idx],plot_type=0,save_yes=True,legend=False)
            data_tj['Z'] = data_tj['Vdev']/data_tj['Jdev']
            idx = idx+1
    
    
    Zs,ReZ,ImZ,Zmag = [],[],[],[]
    C,R,phase = [],[],[]
    for freq in freqs:
        for Vapp in Vapps:
            data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq),delim_whitespace=True)
            data_tj = preprocess_Impedance_data(data_tj,freq)
            comp = get_complex_impedance(data_tj,freq)
            Zs.append(comp)
            ReZ.append(comp.real)
            ImZ.append(comp.imag)
            Zmag.append(abs(comp))
            phase.append(cmath.phase(comp))
            C.append(1/(-2*np.pi*abs(comp)*np.sin(cmath.phase(comp))))
            


    Zs = np.asarray(Zs)

    if plot_output:
        plt.figure(num_fig_tjs+1)
        plt.scatter(ReZ,ImZ)
    # plt.semilogx(freqs,ImZ)
    # plt.semilogx(freqs,C)
        plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'Impedance.jpg',dpi=600,transparent=True)

    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)

    

    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
if __name__ == '__main__':
    main()

   