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
from pathlib import Path
from itertools import repeat
# Don't show warnings
warnings.filterwarnings("ignore")
# Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,zimt_tj_JV_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_JV_sweep

# Main Program
def JV_scan(str2run=[''],Gens=[0],num_fig=0,path2ZimT='',Store_folder='',hasTL=True,PlotName='',MaxVol=1.2):
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

    hasTL : bool
        Changes SCAPS data reading for devices with/without transport layers, by default True

    PlotName: str
        Name of the saved plot JPG, by default = ''

    MaxVol: int, optional
        Changes ZimT's voltage range. Impacts simulation's speed at the cost of data range, by default 1.2
    """     

    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = 14                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
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
    run_simu = True                                         # Rerun simu?
    move_ouput_2_folder = False                                 # Move (True) output of simulation to Store_folder
    clean_output = False                                        # Clean output after simulation (delete files)

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    Vmin = 0                                                # Minimum voltage (unit: V)
    Vmax = MaxVol                                              # Maximum voltage (unit: V)
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

    colors = plt.cm.viridis(np.linspace(0,1,max(len(str2run),4)+1))  # Color range for plotting
   

    # Initialize 
    str_lst,sys_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = -1
    start = time()

    
    if run_simu:
    # Generate tVG files and str_lst
        idx=-1
        for strs in str2run:
            idx=idx+1
            for scan in scans:
                for Gen in Gens:
                    zimt_JV_sweep(Vmin,Vmax,scan,Gen,steps,tVG_name=curr_dir+slash+path2ZimT+'tVG_JV_scan_{:.2e}_G_{:.2e}_{}.txt'.format(scan,Gen,idx))
                    str_lst.append(strs+' -tVG_file tVG_JV_scan_{:.2e}_G_{:.2e}_{}.txt -tj_file tj_JV_scan_{:.2e}_G_{:.2e}_{}.dat'.format(scan,Gen,idx,scan,Gen,idx))
                    sys_lst.append(system)
                    path_lst.append(curr_dir+slash+path2ZimT)
                    tVG_lst.append('tVG_JV_scan_{:.2e}_G_{:.2e}_{}.txt'.format(scan,Gen,idx))
                    tj_lst.append('tj_JV_scan_{:.2e}_G_{:.2e}_{}.dat'.format(scan,Gen,idx))


        # Run ZimT
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_zimt,list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)
        p.close()
        p.join()
        print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds
    else:
        print('/!\ Simulations were not ran again, to rerun simulations set "run_simu = True"')
         
    

    
    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,curr_dir+slash+path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,curr_dir+slash+path2ZimT)
    else:
        Store_folder = ''

    
    ########################################################
    ################## Plotting ############################
    ########################################################
    if plot_JV:
        idx=-1                                                                # index for file names and plot colors
        idc=-1
        for strs,lbl,exp_name in zip(str2run,labels,JVscaps_lst):
            idx = idx + 1
            idc = idc + 1
            ##Import of SCAPS data ##
            if has_TL:
                names_scaps = ['v(V)', 'jtot(mA/cm2)', 'j_total_rec(mA/cm2)', 'j_total_gen(mA/cm2)', '   jbulk(mA/cm2)', 'jifr(mA/cm2)',
                               'jminor_left(mA/cm2)', 'jminor_right(mA/cm2)', 'j_SRH(mA/cm2)', 'j_Radiative(mA/cm2)', 'j_Auger(mA/cm2)']
                data_JVscaps = pd.read_csv(exp_name, names=names_scaps, engine="python", header=None, delim_whitespace=True,
                                       usecols=[0, 1], skiprows=30, skipfooter=10)
            else:
                names_scaps = ['v(V)', 'jtot(mA/cm2)', 'j_total_rec(mA/cm2)', 'j_total_gen(mA/cm2)', '   jbulk(mA/cm2)', 'jifr(mA/cm2)',
                               'jminor_left(mA/cm2)', 'jminor_right(mA/cm2)', 'j_SRH(mA/cm2)', 'j_Radiative(mA/cm2)', 'j_Auger(mA/cm2)']
                data_JVscaps = pd.read_csv(exp_name, names=names_scaps, engine="python", header=None,
                                       delim_whitespace=True, usecols=[0, 1], skiprows=27, skipfooter=10)


            ## Compare SIMsalabim and scaps JVs ##
            plt.figure(0)                           #loads plot
            plt.plot(data_JVscaps['v(V)'], data_JVscaps['jtot(mA/cm2)'], 'o',                   #plots SCAPS data into ZimT plot
                    markeredgecolor=colors[idx], markersize=10, markerfacecolor='None', markeredgewidth=3)
            plt.legend(prop={"size": 20})


            for scan in scans:
                for Gen in Gens:
                        data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_JV_scan_{:.2e}_G_{:.2e}_{}.dat'.format(scan,Gen,idc),delim_whitespace=True)
                        zimt_tj_JV_plot(num_fig_JV,data_tj,x='Vext',y=['Jext'],xlimits = [Vmin-0.1,max(1.4,Vmax+0.1)],ylimits=[-27,2] ,colors=colors[idx],labels=lbl,pic_save_name=curr_dir+slash+path2ZimT+Store_folder+PlotName)


            
    # Make plots
    # plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)

    

    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

if __name__ == '__main__':


    curr_dir = os.getcwd()  # Current working directory
    slash = '/'
    path2ZimT = 'Simulation_program/DDSuite_v400/ZimT'+slash    # Path to SIMsalabim in curr_dir
    ext_save_pic = '.jpg'


    ## Simulation types ##
    MIM_configuration = False
    Pin_nip_configuration = False
    Traps_configuration = False
    Interface_traps_configuration = True


    ## First setup for MIM ##
    if MIM_configuration:
        print('\n')
        print('Start the MIM configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/MIM Configuration')
        str_lst = ['-Nc 5e24 -L 300e-9 -Var_file Var_MIM_ref.dat',
                   '-Nc 5e24 -kdirect 1e-17 -Var_file Var_MIM_1.dat',
                   '-Nc 5e24 -mun_0 1e-7 -mup_0 1e-7 -Var_file Var_MIM_2.dat',
                   '-Nc 5e24 -mup_0 1e-7 -Var_file Var_MIM_3.dat',
                   '-Nc 5e24 -W_L 3.8 -W_R 5 -Var_file Var_MIM_5.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 5']
        JVscaps_lst = [path_scaps / 'testref.iv', path_scaps / 'test1.iv', path_scaps / 'test2.iv',
                       path_scaps / 'test3.iv', path_scaps / 'test5.iv']
        JV_plot_filename = 'JV_MIM_configuration_ZimT'+ext_save_pic
        has_TL = False
        JV_scan(str2run=str_lst,Gens=[4.3e27],num_fig=0,path2ZimT = path2ZimT,Store_folder='',PlotName= JV_plot_filename,MaxVol=1.32)

        # add test number 4
        str_lst = ['-Nc 5e24 -L 100E-9 -Gmax 1.3e28 -Var_file Var_MIM_4.dat']
        labels = ['Test number 4']
        JVscaps_lst = [path_scaps / 'test4.iv']
        JV_scan(str2run=str_lst,Gens=[1.3e28],num_fig=0,path2ZimT = path2ZimT,Store_folder='',PlotName= JV_plot_filename,MaxVol=1.4)
        # add test number 6
        str_lst = ['-Nc 5e24 -VB 4.8 -W_R 4.8 -Var_file Var_MIM_6.dat']
        labels = ['Test number 6']
        JVscaps_lst = [path_scaps / 'test6.iv']
        JV_scan(str2run=str_lst,Gens=[4.28e27],num_fig=0,path2ZimT = path2ZimT,Store_folder='',PlotName= JV_plot_filename,MaxVol=0.82)

        plt.show()



    ## Pin-nip simulation ##
    if Pin_nip_configuration:
        print('\n')
        print('Start the pin-nip configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Pin-nip Configuration')
        str_lst = ['-Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_TL_ref.dat',
                   '-Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -Var_file Var_TL_1.dat',
                   '-Nc 1e24 -L 400e-9 -L_LTL 50e-9 -L_RTL 50e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -Var_file Var_TL_2.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_TL_3.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -Var_file Var_TL_4.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -L 400e-9 -L_LTL 50e-9 -L_RTL 50e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -Var_file Var_TL_5.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2',
                   'Test number 3', 'Test number 4', 'Test number 5']
        JVscaps_lst = [path_scaps / 'TLtestref.iv', path_scaps / 'TLtest1.iv', path_scaps / 'TLtest2.iv',
                       path_scaps / 'TLtest3.iv', path_scaps / 'TLtest4.iv', path_scaps / 'TLtest5.iv']
        JV_plot_filename = 'JV_Pin_nip_configuration_ZimT' + ext_save_pic
        has_TL = True
        JV_scan(str2run=str_lst, Gens=[4.2e27], num_fig=0, path2ZimT=path2ZimT,
                Store_folder='',PlotName=JV_plot_filename,MaxVol=1.45)
        plt.show()


    ## Traps simulation ##
    if Traps_configuration:
        print('\n')
        print('Start the bulk traps configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Traps Configuration')
        str_lst = ['-Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_ref.dat -Tr_type_B 0',
                   '-Bulk_tr 1e20 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_1.dat -Tr_type_B 0',
                   '-Tr_type_B 1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_2.dat',
                   '-Tr_type_B 1 -Bulk_tr 1e20 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_3.dat',
                   '-Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_4.dat',
                   '-Tr_type_B -1 -Bulk_tr 1e20 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_6.dat',
                   '-Etrap 3.7 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_7.dat',
                   '-Etrap 3.9 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_8.dat',
                   '-Etrap 5.1 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_9.dat',
                   '-Etrap 5.1 -Tr_type_B -1 -Bulk_tr 1e22 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_10.dat',
                   '-Etrap 5.1 -Tr_type_B 1 -Bulk_tr 1e22 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_11.dat', ]
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 4', 'Test number 5',
                  'Test number 6', 'Test number 7', 'Test number 8', 'Test number 9', 'Test number 10', 'Test number 11', ]
        JV_files = ['JV_Traps_ref.dat','JV_Traps_1.dat','JV_Traps_2.dat','JV_Traps_3.dat','JV_Traps_4.dat','JV_Traps_6.dat','JV_Traps_7.dat','JV_Traps_8.dat','JV_Traps_9.dat','JV_Traps_10.dat','JV_Traps_11.dat', ]
        Var_files = ['Var_Traps_ref.dat','Var_Traps_1.dat','Var_Traps_2.dat','Var_Traps_3.dat','Var_Traps_4.dat','Var_Traps_6.dat','Var_Traps_7.dat','Var_Traps_8.dat','Var_Traps_9.dat','Var_Traps_10.dat','Var_Traps_11.dat', ]
        JVscaps_lst = [path_scaps / 'traptestref.iv', path_scaps / 'traptest1.iv', path_scaps / 'traptest2.iv',
                     path_scaps / 'traptest3.iv',path_scaps / 'traptest4.iv',
                     path_scaps / 'traptest6.iv',path_scaps / 'traptest7.iv', path_scaps / 'traptest8.iv',
                     path_scaps / 'traptest9.iv', path_scaps / 'traptest10.iv', path_scaps / 'traptest11.iv']
        JV_plot_filename = 'JV_Traps_configuration_ZimT'+ext_save_pic 
        has_TL = True
        JV_scan(str2run=str_lst, Gens=[4.2e27], num_fig=0, path2ZimT=path2ZimT,
                Store_folder='',PlotName=JV_plot_filename,MaxVol=1.35)
        
        str_lst = ['-Cn 1e-12 -Cp 1e-12 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_Traps_5.dat' ]
        labels = ['Test number 5']
        JV_files = ['JV_Traps_5.dat' ]
        Var_files = ['Var_Traps_5.dat' ]
        JVscaps_lst = [path_scaps / 'traptest5.iv']
        JV_scan(str2run=str_lst, Gens=[4.68e27], num_fig=0, path2ZimT=path2ZimT,
                Store_folder='',PlotName=JV_plot_filename,MaxVol=1.35)


        plt.show()

    # Interface traps simulation
    if Interface_traps_configuration:
        print('\n')
        print('Start the interface traps configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Interface Traps Configuration')
        Store_Folder = Path(curr_dir + slash + 'Test simulation/Interface Traps Configuration')
        str_lst = ['-St_L 1e14 -Tr_type_L 0 -St_R 1e14 -Tr_type_R 0 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_ref.dat',
                   '-St_L 1e13 -Tr_type_L 0 -St_R 1e13 -Tr_type_R 0 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_1.dat',
                   '-St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_2.dat',
                   '-St_L 1e14 -Tr_type_L -1 -St_R 1e12 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_3.dat',
                   '-St_L 1e14 -Tr_type_L 1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_4.dat',
                   '-Cn 1e-12 -Cp 1e-12 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_5.dat',
                   '-St_L 1e14 -Tr_type_L 1 -St_R 1e14 -Tr_type_R -1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_6.dat',
                   '-Etrap 3.7 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_7.dat',
                   '-Etrap 3.9 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_8.dat',
                   '-Etrap 5.1 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_9.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -St_L 1e14 -Tr_type_L 0 -St_R 1e14 -Tr_type_R 0 -Nc 1e24 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_10.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -St_L 1e13 -Tr_type_L 0 -St_R 1e13 -Tr_type_R 0 -Nc 1e24 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -Var_file Var_InterfaceTraps_11.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 4', 'Test number 5',
                  'Test number 6', 'Test number 7', 'Test number 8', 'Test number 9', 'Test number 10',
                  'Test number 11', ]
        JV_files = ['JV_InterfaceTraps_ref.dat','JV_InterfaceTraps_1.dat','JV_InterfaceTraps_2.dat','JV_InterfaceTraps_3.dat','JV_InterfaceTraps_4.dat','JV_InterfaceTraps_5.dat','JV_InterfaceTraps_6.dat','JV_InterfaceTraps_7.dat','JV_InterfaceTraps_8.dat','JV_InterfaceTraps_9.dat','JV_InterfaceTraps_10.dat','JV_InterfaceTraps_11.dat', ]
        Var_files = ['Var_InterfaceTraps_ref.dat','Var_InterfaceTraps_1.dat','Var_InterfaceTraps_2.dat','Var_InterfaceTraps_3.dat','Var_InterfaceTraps_4.dat','Var_InterfaceTraps_5.dat','Var_InterfaceTraps_6.dat','Var_InterfaceTraps_7.dat','Var_InterfaceTraps_8.dat','Var_InterfaceTraps_9.dat','Var_InterfaceTraps_10.dat','Var_InterfaceTraps_11.dat', ]
        JVscaps_lst = [path_scaps / 'iftestref.iv', path_scaps / 'iftest1.iv', path_scaps / 'iftest2.iv',
                     path_scaps / 'iftest3.iv',path_scaps / 'iftest4.iv', path_scaps / 'iftest5.iv',
                     path_scaps / 'iftest6.iv',path_scaps / 'iftest7.iv', path_scaps / 'iftest8.iv',
                     path_scaps / 'iftest9.iv', path_scaps / 'iftest10.iv', path_scaps / 'iftest11.iv']
        NPexp_lst = [path_scaps / 'npiftestref.eb', path_scaps / 'npiftest1.eb', path_scaps / 'npiftest2.eb',
                     path_scaps / 'npiftest3.eb',path_scaps / 'npiftest4.eb', path_scaps / 'npiftest5.eb',
                     path_scaps / 'npiftest6.eb',path_scaps / 'npiftest7.eb', path_scaps / 'npiftest8.eb',
                     path_scaps / 'npiftest9.eb', path_scaps / 'npiftest10.eb', path_scaps / 'npiftest11.eb']
        JV_plot_filename = 'JV_Interface_traps_configuration_ZimT'+ext_save_pic 
        has_TL = True
        JV_scan(str2run=str_lst, Gens=[4.2e27], num_fig=0, path2ZimT=path2ZimT,
                Store_folder='',PlotName=JV_plot_filename,MaxVol=1.35)    