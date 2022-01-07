###################################################
###### Simulate JV hysteresis using ZimT ##########
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
def JV_Hyst(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = True, plot_tjs = True, move_ouput_2_folder = True, Store_folder = 'JV_Hyst',clean_output = False,verbose = True):  
    """Run single JV sweep simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_JV_Hyst function (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'JV_Hyst'
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
        path2ZimT = os.path.join(os.getcwd(),'Simulation_program/SIMsalabim_v427/ZimT')                  # Path to ZimT in curr_dir

 
    ## Physics constants
    q = constants.value(u'elementary charge')

    # TPC input
    # see zimt_Hyst in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No JV_Hyst input dictionary given, using default values')
        Vstart = 1.17                                               # Start voltage (unit: V)
        Vfinal = 0                                                  # Final voltage (unit: V)
        scans = np.geomspace(1e-3, 1e4,num=10)                                               # Scan speed (unit: V/s)
        Gens = [0e27]                                             # Average generation rate (unit: m^-3 s^-1) 
        steps = 100                                                 # Number of voltage step
    else:
        if 'Vstart' in input_dic:
            Vstart = input_dic['Vstart']
        if 'Vfinal' in input_dic:
            Vfinal = input_dic['Vfinal']
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
        fixed_str = ' -accDens 0.05'  # add any fixed string to the simulation command

    # Initialize 
    str_lst,code_name_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Initialize 
    code_name_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    size_fig = (16, 12)
    colors = cm.viridis((np.linspace(0,1,max(len(scans),4)+1)) ) # Color range for plotting
    f_JVs = plt.figure(0,figsize=size_fig)
    f_effs = plt.figure(1,figsize=size_fig)

    if run_simu:
    # Generate tVG files and str_lst 
        for scan in scans:
            for Gen in Gens:
                zimt_JV_double_sweep(Vstart,Vfinal,scan,Gen,steps,tVG_name=os.path.join(path2ZimT,'tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen)))
                str_lst.append(fixed_str+' -tVG_file tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt -tj_file tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen,scan,Gen))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                tj_lst.append('tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen))

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
        idx = 0
        Voc_dir1,Jsc_dir1,FF_dir1,PCE_dir1 = [],[],[],[]
        Voc_dir2,Jsc_dir2,FF_dir2,PCE_dir2 = [],[],[],[]
        for scan in scans:
            for Gen in Gens:
                    data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen)),delim_whitespace=True)
                    
                    # Seperate scan directions
                    # abs(data_Var.Vext -Vext) == min(abs(data_Var.Vext -Vext))
                    switch = data_tj.index[abs(data_tj['Vext']-Vfinal)==min(abs(data_tj['Vext']-Vfinal))].tolist() # get index switch scan direction
                    scan_dir1 = data_tj.head(switch[0]+1)
                    scan_dir2 = data_tj.tail(switch[0]+1)
                    if Vstart > Vfinal: # reorder dataframe so that the perf can be calculated
                        scan_dir1 = scan_dir1[::-1].reset_index()
                    else:
                        scan_dir2 = scan_dir2[::-1].reset_index()

                    # Get perfs
                    Jsc_dir1.append(get_Jsc(scan_dir1['Vext'],scan_dir1['Jext'])/10)
                    Voc_dir1.append(get_Voc(scan_dir1['Vext'],scan_dir1['Jext']))
                    FF_dir1.append(get_FF(scan_dir1['Vext'],scan_dir1['Jext']))
                    PCE_dir1.append(get_PCE(scan_dir1['Vext'],scan_dir1['Jext'])/10)
                    Jsc_dir2.append(get_Jsc(scan_dir2['Vext'],scan_dir2['Jext'])/10)
                    Voc_dir2.append(get_Voc(scan_dir2['Vext'],scan_dir2['Jext']))
                    FF_dir2.append(get_FF(scan_dir2['Vext'],scan_dir2['Jext']))
                    PCE_dir2.append(get_PCE(scan_dir2['Vext'],scan_dir2['Jext'])/10)
                    

                    zimt_tj_JV_plot(0,scan_dir1,x='Vext',y=['Jext'],xlimits = [-0.1,1.25],ylimits=[-22,2] ,colors=colors[idx],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=os.path.join(path2ZimT,Store_folder,'Fast_hyst_JV.jpg'),legend=False)
                    zimt_tj_JV_plot(0,scan_dir2,x='Vext',y=['Jext'],xlimits = [-0.1,1.25],ylimits=[-22,2] ,colors=colors[idx],line_type=['--'],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=os.path.join(path2ZimT,Store_folder,'Fast_hyst_JV.jpg'),legend=False)
                    idx = idx + 1 

        plt.figure(1)
        plt.subplot(221)
        plt.semilogx(scans,PCE_dir2,linestyle='-',marker='None',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,PCE_dir1,linestyle='-',marker='None',markersize=10,markerfacecolor='w',label='Reverse')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('PCE [%]')
        plt.ylim([10,25])
        plt.grid(b=True,which='both')
        plt.subplot(222)
        plt.semilogx(scans,FF_dir2,linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,FF_dir1,linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('FF')
        plt.ylim([0.5,1])
        plt.grid(b=True,which='both')
        plt.subplot(223)
        plt.semilogx(scans,Voc_dir2,linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,Voc_dir1,linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('V$_{OC}$ [V]')
        plt.ylim([0.9,1.3])
        plt.grid(b=True,which='both')
        plt.subplot(224)
        plt.semilogx(scans,abs(np.asarray(Jsc_dir2)),linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,abs(np.asarray(Jsc_dir1)),linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('J$_{SC}$ mA cm$^{-2}$')
        plt.ylim([18,26])
        # plt.legend(loc='best', title="Scan direction:",fontsize=35)
        plt.tight_layout()
        plt.grid(b=True,which='both')
        plt.savefig(os.path.join(path2ZimT,Store_folder,'Fast_hyst_eff.jpg'),dpi=100,transparent=True)

        print(PCE_dir1[0]-PCE_dir1[-1])
        

    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))

    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
if __name__ == '__main__':
    JV_Hyst()
    plt.show()
    
    
    

