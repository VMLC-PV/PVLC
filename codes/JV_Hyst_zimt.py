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
from VLC_useful_func import * #sci_notation,run_zimt,zimt_tj_plot,zimt_tj_JV_plot,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_JV_double_sweep
q = constants.value(u'elementary charge')
# Main Program
def main():    
    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                       # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
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
    path2ZimT = 'Simulation_program/DDSuite_v403/ZimT'+slash    # Path to ZimT in curr_dir

    # Physics constants
    q = constants.value(u'elementary charge')

    # Simulation input
    run_simu = True                                         # Rerun simu?
    plot_tjs = True                                        # make plot ?
    plot_output = False
    move_ouput_2_folder = True
    Store_folder = 'JV_hysteresis'+slash
    clean_output = False
    L = 440e-9                                                  # Device thickness (m)
    Gens = [215/(q*400e-9)]                                               # Max generation rate for the gaussian laser pulse
    print(Gens)
    # scans = np.geomspace(1e-3,5e3,num=10)
    
    Vstart = 1.22
    Vfinal = 0

    # scans_exp = [0.01,0.06,0.2,1,2,10,20,100,200,1000,2000]
    # rev_exp = [18.8,19,19.07694,19.15601,19.45013,19.72537,19.81781,19.809,20.20233,20.35866,20.1]
    # for_exp = [18.8,18.8,18.2,17.71846,17.601,17.25488,17.69805,18.618,19.22702,19.85594,19.98172]

    # scans_exp = np.asarray([0.00569,0.0285,0.057,0.286,0.572,2.86,5.72,28.6,57.2,286,572])*4
    # for_exp = [18.59206,17.85167,17.67393,17.33705,17.1059,16.58675,16.54971,17.85843,18.6046,19.73997,19.78995]
    # rev_exp = [18.31429,18.32017,18.61573,19.11934,19.3399,19.70609,19.86326,20.03486,20.06419,20.23221,20.44356]
    # data_exp = pd.read_csv(curr_dir+slash+'summary_Fast_hysteresis_eff.dat',sep=',')
    # print(data_exp)
    # scans_exp = data_exp['speed_mean']*4
    # for_exp = data_exp['PCE_for_mean']
    # rev_exp = data_exp['PCE_rev_mean']
    data_exp = pd.read_csv(curr_dir+slash+'405_pxA_effs.dat',sep=',')
    print(data_exp)
    scans_exp = data_exp['speed_mean']*4
    scans = data_exp['speed_mean']*4
    for_exp = data_exp['PCE_for_mean']
    rev_exp = data_exp['PCE_rev_mean']
    # Initialize 
    sys_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
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
                zimt_JV_double_sweep(Vstart,Vfinal,scan,Gen,100,tVG_name=curr_dir+slash+path2ZimT+'tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                str_lst.append('-L '+str(L)+' -NP 2000 -tVG_file tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt -tj_file tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen,scan,Gen))
                sys_lst.append(system)
                path_lst.append(curr_dir+slash+path2ZimT)
                tVG_lst.append('tVG_JV_Hyst_scan_{:.2e}_G_{:.2e}.txt'.format(scan,Gen))
                tj_lst.append('tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen))

        # print(str_lst)
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
        idx = 0
        Voc_dir1,Jsc_dir1,FF_dir1,PCE_dir1 = [],[],[],[]
        Voc_dir2,Jsc_dir2,FF_dir2,PCE_dir2 = [],[],[],[]
        for scan in scans:
            for Gen in Gens:
                    data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_JV_Hyst_scan_{:.2e}_G_{:.2e}.dat'.format(scan,Gen),delim_whitespace=True)
                    
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
                    

                    zimt_tj_JV_plot(0,scan_dir1,x='Vext',y=['Jext'],xlimits = [-0.1,1.25],ylimits=[-22,2] ,colors=colors[idx],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=curr_dir+slash+path2ZimT+Store_folder+'Fast_hyst_JV.jpg',legend=False)
                    zimt_tj_JV_plot(0,scan_dir2,x='Vext',y=['Jext'],xlimits = [-0.1,1.25],ylimits=[-22,2] ,colors=colors[idx],line_type=['--'],labels=sci_notation(scan, sig_fig=1)+' V s$^{-1}$',pic_save_name=curr_dir+slash+path2ZimT+Store_folder+'Fast_hyst_JV.jpg',legend=False)
                    idx = idx + 1 

        plt.figure(1)
        plt.subplot(221)
        plt.semilogx(scans,PCE_dir2,linestyle='-',marker='None',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,PCE_dir1,linestyle='-',marker='None',markersize=10,markerfacecolor='w',label='Reverse')
        plt.semilogx(scans_exp,for_exp,linestyle='None',marker='o',markersize=10,markerfacecolor='w',label='Forward exp')
        plt.semilogx(scans_exp,rev_exp,linestyle='None',marker='^',markersize=10,markerfacecolor='w',label='Reverse exp')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('PCE [%]')
        plt.ylim([10,25])
        plt.grid(b=True,which='both')
        plt.subplot(222)
        plt.semilogx(scans,FF_dir2,linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,FF_dir1,linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.semilogx(scans_exp,data_exp['FF_for_mean'],linestyle='None',marker='o',markersize=10,markerfacecolor='w',label='Forward exp')
        plt.semilogx(scans_exp,data_exp['FF_rev_mean'],linestyle='None',marker='^',markersize=10,markerfacecolor='w',label='Reverse exp')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('FF')
        plt.ylim([0.5,1])
        plt.grid(b=True,which='both')
        plt.subplot(223)
        plt.semilogx(scans,Voc_dir2,linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,Voc_dir1,linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.semilogx(scans_exp,data_exp['Voc_for_mean'],linestyle='None',marker='o',markersize=10,markerfacecolor='w',label='Forward exp')
        plt.semilogx(scans_exp,data_exp['Voc_rev_mean'],linestyle='None',marker='^',markersize=10,markerfacecolor='w',label='Reverse exp')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('V$_{OC}$ [V]')
        plt.ylim([0.9,1.3])
        plt.grid(b=True,which='both')
        plt.subplot(224)
        plt.semilogx(scans,abs(np.asarray(Jsc_dir2)),linestyle='-',marker='o',markersize=10,markerfacecolor='w',label='Forward')
        plt.semilogx(scans,abs(np.asarray(Jsc_dir1)),linestyle='-',marker='^',markersize=10,markerfacecolor='w',label='Reverse')
        plt.semilogx(scans_exp,data_exp['Jsc_for_mean'],linestyle='None',marker='o',markersize=10,markerfacecolor='w',label='Forward exp')
        plt.semilogx(scans_exp,data_exp['Jsc_rev_mean'],linestyle='None',marker='^',markersize=10,markerfacecolor='w',label='Reverse exp')
        plt.xlabel('Scan speed [V s$^{-1}$]')
        plt.ylabel('J$_{SC}$ mA cm$^{-2}$')
        plt.ylim([18,26])
        # plt.legend(loc='best', title="Scan direction:",fontsize=35)
        plt.tight_layout()
        plt.grid(b=True,which='both')
        plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'Fast_hyst_eff.jpg',dpi=100,transparent=True)

        print(PCE_dir1[0]-PCE_dir1[-1])
        




    
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
    main()
    
    
    

