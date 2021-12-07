#####################################################
####### Test SIMsalabim and ZimT physics ############
#####################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,subprocess,shutil,math,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from time import time
from scipy import stats
from scipy import constants
from pathlib import Path
import scipy.optimize
# Package by VLC
from VLC_units.plots.SimSS_plots import *
import VLC_units.plots.plot_settings_screen
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *
from VLC_units.cleanup_folder.clean_folder import *
from VLC_units.useful_functions.aux_func import *
from VLC_units.impedance.impedance_func import *
from VLC_units.SCLC.SCLC_func import *


## General Inputs
warnings.filterwarnings("ignore")           # Don't show warnings
system = platform.system()                  # Operating system
max_jobs = os.cpu_count()-2                 # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
do_multiprocessing = True                      # Use multiprocessing
if system == 'Windows':                     # cannot easily do multiprocessing in Windows
        max_jobs = 1
        do_multiprocessing = False
        try:                                # kill all running jobs to avoid conflicts
            os.system('taskkill.exe /F /IM simss.exe')
        except:
            pass

# Initialize
curr_dir = os.getcwd()
idx = 0
lines = ['-', '--', '-.', ':']

# Figures control
size_fig = (14, 6)
fig_idx = 0

# Inputs
curr_dir = os.getcwd()                      # Current working directory
path2SIMsalabim = os.path.join(os.getcwd(),'Simulation_program/SIMsalabim_v425/SimSS')    # Path to SIMsalabim in curr_dir
path2ZimT = os.path.join(os.getcwd(),'Simulation_program/SIMsalabim_v425/ZimT')                     # Path to ZimT in curr_dir
ext_save_pic = '.jpg'

# copy and save current device_parameters.txt file so it does not get lost and replace it with the one needed for the test
os.rename(os.path.join(path2SIMsalabim,'device_parameters.txt'), os.path.join(path2SIMsalabim,'device_parameters_saved.txt'))
shutil.copyfile(os.path.join(os.getcwd(),'Test simulation','device_parameters_test_phys_SIMsalabim.txt'), os.path.join(path2SIMsalabim,'device_parameters.txt'))
os.rename(os.path.join(path2ZimT,'device_parameters.txt'), os.path.join(path2ZimT,'device_parameters_saved.txt'))
shutil.copyfile(os.path.join(os.getcwd(),'Test simulation','device_parameters_test_phys_zimt.txt'), os.path.join(path2ZimT,'device_parameters.txt'))

shutil.copyfile(os.path.join(os.getcwd(),'Test simulation','blue.txt'), os.path.join(path2SIMsalabim,'blue.txt'))
shutil.copyfile(os.path.join(os.getcwd(),'Test simulation','red.txt'), os.path.join(path2SIMsalabim,'red.txt'))





# Chose test to run:
Test_gen_profile = True
Test_SCLC_MottGurney = True
Test_SCLC_Traps = True
Test_TPV = True
Test_TPC = True
Test_RCtime = True
Test_Impedance = True


### Test generation profile
if Test_gen_profile:
    # Test whether the generation profile is correclty inputed in SIMsalabim
    print('\n')
    print('Start the Test generation profile:')
    L = 300e-9 # m
    L_LTL = 20e-9 # m
    L_RTL = 20e-9 # m
    TLsAbsorbs  = [0,1] # only 0 or 1
    Store_Folder = os.path.join(curr_dir,'Test simulation')
    gen_profiles = ['None','blue.txt','red.txt']
    labels = ['No profile', 'Blue profile', 'Red profile']
    JV_files = ['JV_no_profile.dat','JV_blue.dat','JV_red.dat']
    Var_files = ['Var_no_profile.dat','Var_blue.dat','Var_red.dat']
    code_name_lst, path_lst = [], []

    f_gen_pro,ax = plt.subplots(1,3,num =fig_idx, figsize=size_fig)
    fig_idx = fig_idx + 1
    subplot_num = 0
    for TLsAbsorb in TLsAbsorbs:
        # f_gen_pro = plt.figure(fig_idx, figsize=size_fig)
                
        str_lst = ['-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -Gen_profile  '+gen_profiles[0],
        '-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -Gen_profile  '+gen_profiles[1],'-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -grad 10 -Gen_profile  '+gen_profiles[2]]
        # print(str_lst)
        # Simulation input
        run_simu = True     
        start = time()
        count = 0
        for i in str_lst:
            str_lst[count] = str_lst[count] +' -JV_file '+ JV_files[count] +' -Var_file '+ Var_files[count]
            code_name_lst.append('SimSS')
            path_lst.append(path2SIMsalabim)
            count = count + 1
        
        # Color range for plotting
        colors = plt.cm.viridis(np.linspace(0, 1, max(len(str_lst), 4) + 1))
        
        # Run simulation
        if run_simu:
            if do_multiprocessing:
                run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
            else:
                for i in range(len(str_lst)):
                    run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)
        print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds
        print('Now plotting the results...') # Time in seconds
        
        # Plot
        for gen,jv,var,color,lbl in zip(gen_profiles,JV_files,Var_files,colors,labels):
            if gen != 'None':
                data_var = make_df_Var(os.path.join(path2SIMsalabim, var) )  
                data_gen = pd.read_csv(os.path.join(path2SIMsalabim , gen),delim_whitespace=True,names=['x','Gehp'])
                if TLsAbsorb == 0: # correct profile to account for the part that is cut by SIMsalabim
                    data_var =  data_var[data_var['x'] >= L_LTL]
                    data_var =  data_var[data_var['x'] <= L-L_RTL]
                    data_gen['x'] = data_gen['x']/max(data_gen['x'])
                    data_gen = data_gen[data_gen['x'] >= L_LTL/L]
                    data_gen = data_gen[data_gen['x'] <= (L-L_RTL)/L]
                    data_gen['x'] = (data_gen['x']-min(data_gen['x']))/max(data_gen['x']-min(data_gen['x']))

                ax[subplot_num].plot(data_gen['x']/max(data_gen['x']),data_gen['Gehp']/data_gen['Gehp'].mean(),color=color,linestyle='None',marker='x',markeredgecolor=color,markersize=10,markerfacecolor='None',markeredgewidth = 1,markevery=5)
                ax[subplot_num].plot((data_var['x']-min(data_var['x']))/max(data_var['x']-min(data_var['x'])),data_var['Gehp']/data_var['Gehp'].mean(),color=color,label=lbl,linestyle='-')
            else:
                data_var = make_df_Var(os.path.join(path2SIMsalabim, var))
                if TLsAbsorb == 0: # correct profile to account for the part that is cut by SIMsalabim
                    data_var =  data_var[data_var['x'] >= L_LTL]
                    data_var =  data_var[data_var['x'] <= L-L_RTL]

                ax[subplot_num].plot(data_var['x']/max(data_var['x']),data_var['Gehp']/data_var['Gehp'].mean(),color=color,label=lbl,linestyle='-')
            data_jv = make_df_JV(os.path.join(path2SIMsalabim, jv))
            ax[2].plot(data_jv['Vext'],data_jv['Jphoto']/10,color=color,linestyle=lines[subplot_num],label='TLAbsorb = {:.0f}'.format(TLsAbsorb))
        subplot_num = subplot_num + 1
        if TLsAbsorb == 0:
            ax[0].set_title('TLAbsorb = 0')
            ax[0].set_xlabel('Norm. position')
            ax[0].set_ylabel('Norm. Generation rate')
            ax[0].legend(loc='best',frameon=False)
        else:
            ax[1].set_title('TLAbsorb = 1')
            ax[1].set_xlabel('Norm. position')
            ax[1].set_ylabel('Norm. Generation rate')
        legend_elements = [
            Line2D([0], [0], color=colors[0],linestyle= lines[0], label='TLAbsorb = 0'),
            Line2D([0], [0],  color=colors[0],linestyle= lines[0], label='TLAbsorb = 1')]
        ax[2].legend(handles= legend_elements ,loc='best',frameon=False)
        ax[2].set_xlabel('Voltage [V]')
        ax[2].set_ylabel('Jphoto [mA cm$^{-2}$]')
        
    plt.tight_layout()
    plt.savefig(os.path.join(Store_Folder, 'test_gen_pro.jpg'),dpi=100,transparent=True)
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds


    

### Test textbook SCLC
# Test fit mobility with Mott-Gurney equation
if Test_SCLC_MottGurney:
    # from SCLC_func import fit_MottGurney,MottGurney,log_slope
    # Test whether the generation profile is correclty inputed in SIMsalabim
    print('\n')
    print('Start the Test fit Mott-Gurney to SCLC:')
    L = 300e-9 # m
    L_LTL = 0e-9 # m
    L_RTL = 0e-9 # m
    W_L = 3.9 # eV
    W_R = 3.9 # eV
    eps_r = 3.5
    Store_Folder = os.path.join(curr_dir,'Test simulation')
    mun_0s = [1e-8,1e-7,1e-4]
    labels, JV_files, Var_files, str_lst = [],[],[],[]
    code_name_lst, path_lst = [], []
    for i in mun_0s:
        labels.append('{:.1e}'.format(i))
        JV_files.append('JV_mu_{:.1e}.dat'.format(i))
        Var_files.append('Var_mu_{:.1e}.dat'.format(i))
        str_lst.append( '-W_L '+str(W_L)+' -W_R '+str(W_R)+' -eps_r '+str(eps_r)+' -L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -Rseries 0 -Gehp 0 -Vmin 0.001 -Vmax 25 -Vacc -0.1 -Vdistribution 2 -NJV 50 -until_Voc 0 -tolJ 5e-6 -mun_0 {:.1e} -mup_0 {:.1e} -JV_file JV_mu_{:.1e}.dat -Var_file Var_mu_{:.1e}.dat'.format(i,i,i,i))
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)

    f_gen_pro,ax = plt.subplots(1,3,num =fig_idx, figsize=size_fig)
    fig_idx = fig_idx + 1
    subplot_num = 0

    # Color range for plotting
    colors = plt.cm.viridis(np.linspace(0, 1, max(len(str_lst), 4) + 1))
    # Simulation input
    run_simu = True     
    start = time()
    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds
    print('Now plotting the results...') # Time in seconds

    # Plot
    fitted_mu = []
    for mun_0,jv,var,color,lbl in zip(mun_0s,JV_files,Var_files,colors,labels):
        data_jv = make_df_JV(os.path.join(path2SIMsalabim, jv))
        ax[0].loglog(data_jv['Vext'],data_jv['Jext']/10,color=color,label=lbl)
        result = fit_MottGurney(data_jv['Vext'],data_jv['Jext'],mu=1e-8,eps_r=eps_r,Vbi=W_L-W_R,L=L-L_LTL-L_RTL)
        mu_fit = result[0]
        fitted_mu.append(mu_fit)
        ax[0].loglog(data_jv['Vext'],MottGurney(data_jv['Vext'],mu_fit,eps_r,W_L-W_R,L)/10,color=color,linestyle='None',marker='x',markeredgecolor=color,markersize=10,markerfacecolor='None',markeredgewidth = 1,markevery=5)
        # slope = np.diff(np.log(data_jv['Jext']))/np.diff(np.log(data_jv['Vext']))
        slope = log_slope(np.asarray(data_jv['Vext']),np.asarray(data_jv['Jext']))
        ax[1].semilogx(data_jv['Vext'],slope,color=color,label=lbl)
    ax[2].loglog([min(mun_0s)/2,max(mun_0s)*5],[min(mun_0s)/2,max(mun_0s)*5],'k-')
    ax[2].loglog(mun_0s,fitted_mu,'o')
    
    # ax[0].set_ylim([1e-3,1e4])
    ax[0].set_xlabel('Voltage [V]')
    ax[0].set_ylabel('J [mA cm$^{-2}$]')
    ax[0].legend(loc='best',frameon=False)
    ax[1].set_xlabel('Voltage [V]')
    ax[1].set_ylabel('Slope')
    ax[2].set_xlabel('True value')
    ax[2].set_ylabel('Fitted Value')
    plt.tight_layout()
    plt.savefig(os.path.join(Store_Folder, 'test_SCLC_MottGurney.jpg'),dpi=100,transparent=True)
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

# test to get trap density from Vtfl
if Test_SCLC_Traps:
    # from SCLC_func import* 
    # Test whether the generation profile is correclty inputed in SIMsalabim
    print('\n')
    print('Start the Test SCLC with traps:')
    L = 300e-9 # m
    L_LTL = 0e-9 # m
    L_RTL = 0e-9 # m
    W_L = 3.9 # eV
    W_R = 3.9 # eV
    eps_r = 3.5 # relative dielectric constant
    Nc = 1e24
    mun_0 = 1e-4 # m^2/Vs, zero field mobility
    Tr_type_B = -1 # Trap type of bulk and grain boundary traps: -1: acceptor, 0: neutral, 1: donor
    Bulk_trs = [5e21,6e21,7e21,8e21]
    Store_Folder = os.path.join(curr_dir,'Test simulation')
    
    labels, JV_files, Var_files, str_lst = [],[],[],[]
    code_name_lst, path_lst = [], []
    for i in Bulk_trs:
        labels.append('{:.1e}'.format(i))
        JV_files.append('JV_Bulk_tr_{:.1e}.dat'.format(i))
        Var_files.append('Var_Bulk_tr_{:.1e}.dat'.format(i))
        str_lst.append( '-W_L '+str(W_L)+' -W_R '+str(W_R)+' -eps_r '+str(eps_r)+' -L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -mun_0 '+str(mun_0)+' -mup_0 '+str(mun_0)+' -Tr_type_B '+str(Tr_type_B)+' -Nc '+str(Nc)+' -Rseries 0 -Gehp 0 -Vmin 1e-3 -Vmax 30 -Vdistribution 2 -NJV 100 -tolJ 1e-4 -NP 1000 -Vacc -0.1 -until_Voc 0 -Bulk_tr {:.1e} -JV_file JV_Bulk_tr_{:.1e}.dat -Var_file Var_Bulk_tr_{:.1e}.dat'.format(i,i,i))
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)
    
    # print(str_lst)
    f_gen_pro,ax = plt.subplots(1,3,num =fig_idx, figsize=size_fig)
    fig_idx = fig_idx + 1
    subplot_num = 0

    # Color range for plotting
    colors = plt.cm.viridis(np.linspace(0, 1, max(len(str_lst), 4) + 1))
    # Simulation input
    run_simu = True     
    start = time()
    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds
    print('Now plotting the results...') # Time in seconds

    # Plot
    fitted_Bulk_tr = []
    for Bulk_tr,jv,var,color,lbl in zip(Bulk_trs,JV_files,Var_files,colors,labels):
        data_jv = make_df_JV(os.path.join(path2SIMsalabim, jv))
        ax[0].loglog(data_jv['Vext'],data_jv['Jext']/10,color=color,label=lbl)

        V_slope,J_slope,slopes,get_tangent,idx_max,max_slopes,tang_val_V1,tang_val_V2,tang_val_V3,V1,J1,V2,J2 = SCLC_get_data_plot(data_jv['Vext'],data_jv['Jext'])
        if get_tangent == 1:
            Nt_vtfl1 = calc_net_charge(V1,L,eps_r)
            Nt_vtfl2 = calc_net_charge(V2,L,eps_r)    
            Nt_vtfl_inf = calc_net_charge(V_slope[idx_max],L,eps_r)
        fitted_Bulk_tr.append(Nt_vtfl2)
        if get_tangent == 1:
            ax[0].loglog(V_slope,tang_val_V1/10,'--',color=color)
            ax[0].loglog(V_slope,tang_val_V2/10,'--',color=color)
            ax[0].loglog(V_slope,tang_val_V3/10,'--',color=color)
            ax[0].loglog(V_slope[idx_max],J_slope[idx_max]/10,'b^',markersize=5)
            ax[0].axvline(x=calc_vnet_with_ions(0,Bulk_tr,L,eps_r), color=color, linestyle='-')

        ax[0].loglog(V1,J1/10,'rs',markersize=5)
        ax[0].loglog(V2,J2/10,'mo',markersize=5)
        slope = log_slope(np.asarray(data_jv['Vext']),np.asarray(data_jv['Jext']))
        ax[1].semilogx(data_jv['Vext'],slope,color=color,label=lbl)
    ax[2].loglog([min(Bulk_trs)/2,max(Bulk_trs)*5],[min(Bulk_trs)/2,max(Bulk_trs)*5],'k-')
    ax[2].loglog(Bulk_trs,fitted_Bulk_tr,'o')
    
    ax[0].set_xlabel('Voltage [V]')
    ax[0].set_ylabel('J [mA cm$^{-2}$]')
    ax[0].legend(loc='best',frameon=False)
    ax[1].set_xlabel('Voltage [V]')
    ax[1].set_ylabel('Slope')
    ax[2].set_xlabel('True value')
    ax[2].set_ylabel('Fitted Value')
    plt.tight_layout()
    plt.savefig(os.path.join(Store_Folder, 'test_SCLC_Traps.jpg'),dpi=100,transparent=True)
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
    


# test TPV compare Voc from ZimT and SIMsalabim
if Test_TPV:
    print('\n')
    print('Start the Test TPV compare Voc from ZimT and SIMsalabim:')
    # from tVG_gen import zimt_light_decay
    # Simulation input
    run_simu = True                                        # Rerun simu?
    plot_tjs = True                                        # make plot ?
    plot_output = False
    move_ouput_2_folder = True
    Store_folder = 'TPV'
    clean_output = False
    L = 140e-9                                                  # Device thickness (m)
    Gens = [1e28]                                               # Max generation rate for the gaussian laser pulse
    G0s = [1e27]#,1e25]
    tpulse = 5e-8

    # Initialize 
    code_name_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    
    colors = cm.viridis((np.linspace(0,1,max(len(Gens),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(fig_idx,figsize=size_fig)
    
    if run_simu:           
    
        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for G0 in G0s:
                zimt_light_decay(1e-8,5e-5,Gen,G0,'oc',steps=100,trf = 20e-9,time_exp =True,tVG_name=os.path.join(path2ZimT,'tVG_TPV_G_{:.1e}_G0_{:.1e}.txt'.format(Gen,G0)))
                str_lst.append('-L '+str(L)+' -tVG_file tVG_TPV_G_{:.1e}_G0_{:.1e}.txt -tj_file tj_TPV_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0,Gen,G0))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_TPV_G_{:.1e}_G0_{:.1e}.txt'.format(Gen,G0))
                tj_lst.append('tj_TPV_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0))
        

        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose = True)

        print('ZimT calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    # Check Voc with SIMsalabim
    print('Calculate Steady-state:') # Time in seconds
    str_lst,JV_files,code_name_lst,path_lst,scPars_files = [],[],[],[],[]
    Gehps = Gens + G0s

    for Gen in Gehps:    
        str_lst.append('-L '+str(L)+' -Gehp {:.1e} -scPars_file sc_G_{:.1e}.dat -JV_file JV_TPV_G_{:.1e}.dat'.format(Gen,Gen,Gen))
        JV_files.append('JV_TPV_G_{:.1e}.dat'.format(Gen))
        scPars_files.append('sc_G_{:.1e}.dat'.format(Gen))
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)

    
    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    
    ########################################################
    ################## tjs_file ############################
    ########################################################
    if plot_tjs:
        print('Now plotting the results...') 
        plt.figure(fig_idx)
        # colors=['k']
        for scPars_file in scPars_files:
            data_sc = pd.read_csv(os.path.join(path2SIMsalabim,scPars_file),delim_whitespace=True)
            plt.axhline(y=float(data_sc['Voc']),color="black", linestyle="--")

        for G0 in G0s:
            idx = 0
            for Gen in Gens:
                data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_TPV_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0)),delim_whitespace=True)
                zimt_Voltage_transient_plot(fig_idx,data_tj,y=['Vext'],xlimits=[],colors=colors[idx],plot_type=0,save_yes=True,pic_save_name = os.path.join(curr_dir,'Test simulation','test_TPV.jpg')) 
                idx = idx + 1
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
    fig_idx = fig_idx + 1
    

# test TPC compare Jsc from ZimT and SIMsalabim
if Test_TPC:
    print('\n')
    print('Start the Test TPC compare Jsc from ZimT and SIMsalabim:')
    # from tVG_gen import zimt_light_decay
    # Simulation input
    run_simu = True                                        # Rerun simu?
    plot_tjs = True                                        # make plot ?
    plot_output = False
    move_ouput_2_folder = True
    Store_folder = 'TPC'
    clean_output = False
    L = 140e-9                                                  # Device thickness (m)
    Gens = [1e28]                                               # Max generation rate for the gaussian laser pulse
    G0s = [1e27]#,1e25]
    tpulse = 5e-8

    # Initialize 
    code_name_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    colors = cm.viridis((np.linspace(0,1,max(len(Gens),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(fig_idx,figsize=size_fig)

    if run_simu:           
    
        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for G0 in G0s:
                zimt_light_decay(1e-8,1e-5,Gen,G0,0,100,trf = 20e-9,time_exp =True,tVG_name=os.path.join(path2ZimT,'tVG_TPC_G_{:.1e}_G0_{:.1e}.txt'.format(Gen,G0)))
                str_lst.append('-L '+str(L)+' -tVG_file tVG_TPC_G_{:.1e}_G0_{:.1e}.txt -tj_file tj_TPC_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0,Gen,G0))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_TPC_G_{:.1e}_G0_{:.1e}.txt'.format(Gen,G0))
                tj_lst.append('tj_TPC_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0))
        
 
        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose = True)

        print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    # Check Voc with SIMsalabim
    print('Calculate Steady-state:') # Time in seconds
    str_lst,JV_files,code_name_lst,path_lst,scPars_files = [],[],[],[],[]
    Gehps = Gens + G0s

    for Gen in Gehps:    
        str_lst.append('-L '+str(L)+' -Gehp {:.1e} -scPars_file sc_G_{:.1e}.dat -JV_file JV_TPC_G_{:.1e}.dat'.format(Gen,Gen,Gen))
        JV_files.append('JV_TPC_G_{:.1e}.dat'.format(Gen))
        scPars_files.append('sc_G_{:.1e}.dat'.format(Gen))
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)

    
    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    
    ########################################################
    ################## tjs_file ############################
    ########################################################
    if plot_tjs:
        print('Now plotting the results...')
        plt.figure(fig_idx)
        # colors=['k']
        for scPars_file in scPars_files:
            data_sc = pd.read_csv(os.path.join(path2SIMsalabim,scPars_file),delim_whitespace=True)
            plt.axhline(y=float(data_sc['Jsc']/10),color="black", linestyle="--")

        for G0 in G0s:
            idx = 0
            for Gen in Gens:
                data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_TPC_G_{:.1e}_G0_{:.1e}.dat'.format(Gen,G0)),delim_whitespace=True)
                zimt_tj_plot(fig_idx,data_tj,y=['Jext'],colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name = os.path.join(curr_dir,'Test simulation','test_TPC.jpg'))
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
    fig_idx = fig_idx + 1




# test RC-time decay
if Test_RCtime:
    print('\n')
    print('Start the Test RC-time decay:')
    # Simulation input
    run_simu = True                                        # Rerun simu?
    plot_tjs = True                                        # make plot ?
    plot_output = False
    move_ouput_2_folder = True
    Store_folder = 'RC_decay'
    clean_output = False
    L = 100e-9 # m
    L_LTL = 0e-9 # m
    L_RTL = 0e-9 # m
    CB = 4 # eV
    VB = 6 # eV
    W_L = 5 # eV
    W_R = 5 # eV
    eps_r = 3.5 # relative dielectric constant
    mu = 1e-7 # m^2/Vs, zero field mobility
    Rseries = 3e-4
    Cgeo = (eps_r*eps_0/L) 
    Vstarts = [0]                                               
    Vfinals = [1]
    tpulse = 5e-8

    # Initialize 
    code_name_lst,str_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    idx = 0
    start = time()

    # Figures control
    colors = cm.viridis((np.linspace(0,1,max(len(Vstarts),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(fig_idx,figsize=size_fig)

    if run_simu:           
        # from tVG_gen import zimt_voltage_step
        # Generate tVG files and str_lst for light pulse simulation
        for Vstart in Vstarts:
            for Vfinal in Vfinals:
                zimt_voltage_step(1e-8,1e-6,Vstart,Vfinal,0,steps=100,trf = 10e-9,time_exp =True,tVG_name=os.path.join(path2ZimT,'tVG_RCtime_Vstart_{:.1f}_Vend_{:.1f}.txt'.format(Vstart,Vfinal)))
                str_lst.append('-W_L '+str(W_L)+' -W_R '+str(W_R)+' -eps_r '+str(eps_r)+' -L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -mun_0 '+str(mu)+' -mup_0 '+str(mu)+' -CB '+str(CB)+' -VB '+str(VB)+' -Rseries '+str(Rseries)+' -Gehp 0 -tVG_file tVG_RCtime_Vstart_{:.1f}_Vend_{:.1f}.txt -tj_file tj_RCtime_Vstart_{:.1f}_Vend_{:.1f}.dat'.format(Vstart,Vfinal,Vstart,Vfinal))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_RCtime_Vstart_{:.1f}_Vend_{:.1f}.txt'.format(Vstart,Vfinal))
                tj_lst.append('tj_RCtime_Vstart_{:.1f}_Vend_{:.1f}.dat'.format(Vstart,Vfinal))
        
 
        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose = True)

        print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    
    ########################################################
    ################## tjs_file ############################
    ########################################################
    if plot_tjs:        
        for Vfinal in Vfinals:
            idx = 0
            for Vstart in Vstarts:
                data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_RCtime_Vstart_{:.1f}_Vend_{:.1f}.dat'.format(Vstart,Vfinal)),delim_whitespace=True)

                data2fit = data_tj #.copy()
                t_pic = data2fit['t'][data_tj['Jext'].idxmax()]
                data2fit['t'] = data2fit['t']-t_pic
                data2fit = data2fit[data2fit.t >= 0]
                # perform the fit
                p0 = (Rseries*Cgeo, max(data2fit['Jext']), min(data2fit['Jext'])) # start with values near those we expect
                params, cv = scipy.optimize.curve_fit(MonoExpDecay,data2fit['t'], data2fit['Jext'], p0)
                k, A, B = params

                plt.plot(data2fit['t'], MonoExpDecay(data2fit['t'], k,A,B)/10, '--', label='RC-fit = {:.2e}'.format(k))
                
                zimt_tj_plot(fig_idx,data_tj,y=['Jext'],colors=colors[idx],plot_type=0,save_yes=True,legend=True,labels='RC-input = {:.2e}'.format(Rseries*Cgeo),pic_save_name = os.path.join(curr_dir,'Test simulation','test_RC_decay.jpg'))

# test Impedance simple capacitor
if Test_Impedance:
    print('\n')
    print('Start the Test Impedance:')
    import cmath
    # from tVG_gen import zimt_impedance
    from scipy import interpolate
    from impedance.models.circuits import circuits as impcir
    from impedance import preprocessing as imppre
    from impedance import visualization as impvis
    curr_dir = os.getcwd()                           # Current working directory

    ## Physics constants
    q = constants.value(u'elementary charge')
    eps_0 = constants.epsilon_0

    ## Simulation input
    run_simu = True                                           # Rerun simu? (bool)
    plot_tjs = False                                           # make tJ and tV plots? (bool)
    plot_output = 1                                        # make Nyquist and Bode plots? (bool)
    move_ouput_2_folder = True                             # (bool)
    Store_folder = 'Impedance'
    clean_output = False                                   # Make output plots? (bool)
    calc_capacitance = 1                                   # calculate capacitance? (bool)
    nFcm2 = 1                                              # Capacitance unit (bool)

    L = 100e-9                                             # Device thickness (m)
    eps = 4                                              # active layer permittivity, same as in device_parameters.txt
    L_LTL = 0e-9 # m
    L_RTL = 0e-9 # m
    CB = 4 # eV
    VB = 6 # eV
    W_L = 5 # eV
    W_R = 5 # eV
    mu = 1e-7 # m^2/Vs, zero field mobility
    Rseries = 0e-4 # ohm m^2
    C_geo = (eps*eps_0/L)                             # geometric capacitance (Ohm)
    freqs1 = np.geomspace(1e2,1e4,num=15,endpoint=False)
    freqs2 = np.geomspace(1e4,1e9,num=20,endpoint=True)
    freqs = np.append(freqs1, freqs2)                      # frequencies to simulate (Hz)
    freqs_interp = np.geomspace(1e2,1e7,num=2000)          # same but more dense range of frequencies (for interpolation)
    # print(freqs)
    Vapps = [0]      # Applied voltage (V)
    Vamp = 0.01                       # Amplitude voltage perturbation
    Gen = 0e27                        # Average generation rate 
    sun = 1                      # generation rate at 1 sun

    ## Figure control
    plottitle = 'Test Impedance' #'ZimT: Organic Solar Cell  {:.0f}nm  {:.1f}sun'.format(L*1e9, Gen/sun)  # full title of plots
    savetitle = 'Test_Impedance'    # start of filename of plots
    colors = cm.viridis((np.linspace(0, 1, max(len(freqs), 4)+1)))  # Freq colors (virid-ish)
    colors1 = cm.winter(np.linspace(0, 1, max(len(Vapps), 4)+1))    # Vapp colors (blue-ish)
    colors1[1] = colors1[0]  # Color adjustment for dark colors: skip second color
    colors1 = colors1[1:]    # same
    color1 = colors1[0]      # Color of left axis


    ## Initialize 
    str_lst,code_name_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    start = time()
    print('Vapps:', str(Vapps))

    ###########################################################################
    ### Run ZimT ##############################################################
    ###########################################################################

    if run_simu:
        print('sim loop...')
        # Generate tVG files 
        for freq in freqs:
            for Vapp in Vapps:
                
                zimt_impedance(Vapp,Vamp,freq,Gen,steps=200,tVG_name=os.path.join(path2ZimT,'tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq)))
                str_lst.append('-W_L '+str(W_L)+' -W_R '+str(W_R)+' -eps_r '+str(eps)+' -L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -mun_0 '+str(mu)+' -mup_0 '+str(mu)+' -CB '+str(CB)+' -VB '+str(VB)+' -Rseries '+str(Rseries)+' -tVG_file '+'tVG_{:.2f}V_f_{:.1e}Hz.txt '.format(Vapp,freq)+'-tj_file '+'tj_{:.2f}V_f_{:.1e}Hz.dat '.format(Vapp,freq))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq))
                tj_lst.append('tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq))
        
        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose=True)
    
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in ,seconds

    
    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    ###########################################################################
    ### Calculate Complex Impedance ###########################################
    ###########################################################################

    Zs,ReZ,ImZ,Zmag = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    Cap,R,phase = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    fZ_file = [[] for f in Vapps]
    Jm, tm, Jmo, tmo = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    print('Z loop...')
    for idx in range(len(Vapps)):
        for freq in freqs:
            with open(os.path.join(path2ZimT,Store_folder ,'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapps[idx],freq)), 'r') as file:# open in read mode
                data_tj2 = pd.read_csv(file, delim_whitespace=True)
            Jmo[idx].append(data_tj2['Jext'])  # save original data for comparison
            tmo[idx].append(data_tj2['t'])
            data_tj2 = preprocess_Impedance_data(data_tj2,freq)
            Jm[idx].append(data_tj2['Jext'])  # save preprocessed data for comparison
            tm[idx].append(data_tj2['t'])
            comp = get_complex_impedance(data_tj2,freq)
            ## save
            Zs[idx].append(comp)
            ReZ[idx].append(comp.real)
            ImZ[idx].append(comp.imag)
            Zmag[idx].append(abs(comp))
            phase[idx].append(cmath.phase(comp))
            Cap[idx].append((1/comp).imag/(2*np.pi*freq))
        ## save f and Z in file for later use of impedance.py
        fZ = np.transpose([freqs, np.asarray(ReZ[idx]), np.asarray(ImZ[idx])])
        fZ = fZ[np.asarray(ImZ[idx])<0]  # only save realistic data: Im(Z)<0
        fZ_file[idx] = os.path.join(path2ZimT,Store_folder,'fZ_{:.2f}V.txt'.format(Vapps[idx]))
        np.savetxt(fZ_file[idx], fZ, delimiter=',')  # ',' needed for impedance.preprocessing

    ## convert lists to arrays to enable calculations
    Zs, ReZ, ImZ, Cap = np.asarray(Zs), np.asarray(ReZ), np.asarray(ImZ), np.asarray(Cap)

      
    ###########################################################################
    ### J/t, V/t graphs #######################################################
    ###########################################################################

    if plot_tjs:
        print('tJ loop...')
        lines = ['-', '--', ':', '-.', 'x-', 'x--', 'x:', 'x-.', 'o-', 'o--', 'o:',
                'o-.', '+-', '+--', '+-.', '+:']  # linestyles for up to 16 Vapps
        for idx in range(len(Vapps)):
            for i in range(len(freqs)):
                with open(os.path.join(path2ZimT,Store_folder,'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapps[idx],freqs[i])), 'r') as file:
                    data_tj = pd.read_csv(file, delim_whitespace=True)
                ## norm data to interval (-1, 1)
                data_tj['Jext_norm'] = data_tj['Jext']/max(data_tj['Jext'])*10
                data_tj['Vext_norm'] = data_tj['Vext']/max(data_tj['Vext'])*10
                data_tj['t_norm'] = data_tj['t']/max(data_tj['t'])*10
                if i%int((len(freqs)-1)/2) == 0:    # plot 3 lines per Vapp
                    ## plot normed J/t
                    zimt_tj_plot(100, data_tj, x='t_norm', y=['Jext_norm'],
                                ylimits= [-1.1,1.1],
                                labels=sci_notation(freqs[i],sig_fig=0)+' Hz, '+str(Vapps[idx])+' V',
                                colors=colors[i], line_type=[lines[idx]],
                                plot_type=0, save_yes=True, legend=True)
                    ## plot normed V/t
                    zimt_tj_plot(101, data_tj, x='t_norm', y=['Vext_norm'],
                                ylimits= [-1.1,1.1],
                                labels=sci_notation(freqs[i],sig_fig=0)+' Hz '+str(Vapps[idx])+' V',
                                colors=colors[i], line_type=[lines[idx]],
                                plot_type=0, save_yes=True, legend=True)

                ## check whether the sin fit works
                if i%int((len(freqs)-1)/2) == 0:  # test sin fit of 3 curves per Vapp
                    try:
                        Jf_amp, Jf_freq, Jf_phi, Jf_off = fit_sin_func(np.asarray(tm[idx][i]), np.asarray(Jm[idx][i]), freqs[i])
                        Jf = Jf_amp*np.sin(Jf_freq*2*np.pi*tm[idx][i] + Jf_phi) + Jf_off  # fitted J
                        plt.figure()
                        plt.plot(tmo[idx][i], Jmo[idx][i], 'o', label='original')
                        plt.plot(tm[idx][i], Jm[idx][i], 'o', label='preprocessed')
                        plt.plot(tm[idx][i], Jf, label='fitted')
                        plt.xlabel('t')
                        plt.ylabel('J')
                        plt.legend()
                        plt.title('freq {:.1e}Hz Vapp {:.1f}V'.format(freqs[i], Vapps[idx]))
                    except:
                        print('tJ fit plot failed.')

    ###########################################################################
    ### Plots #################################################################
    ###########################################################################

    if plot_output:
        print('plotting...')

        ## Bode plot C/f
        if nFcm2:
            Cap = Cap/(1e-9/1e-4)  # to get nF/cm^2
            C_geo = C_geo/(1e-9/1e-4)  # to get nF/cm^2

                    
        try:
            fig4, ax5 = plt.subplots(figsize=size_fig)
            ax5.set_xlabel('Frequency  f  [Hz]')
            ax5.set_ylabel('Capacitance  C  [F]')
            if nFcm2:
                ax5.set_ylabel(r'Capacitance  C  [nF/cm$^2$]')
                ax5.set_ylim(-5, 83)
                ax5.set_xlim(1e2, 1e9)
            else:
                ax5.set_yscale('log')
            for idx in range(len(Vapps)):
                ax5.semilogx(freqs, Cap[idx], 'o-', color=colors1[idx], label='{:.1f} V, ZimT'.format(Vapps[idx]))

            ax5.axhline(C_geo, c='k', label='geometric capacitance')
            plt.legend(fontsize='small')
            plt.title(plottitle)
            plt.savefig(os.path.join(curr_dir,'Test simulation','Test_Impedance_C_f.png'),
                        dpi=100,transparent=True, bbox_inches='tight')
        except:
            print('Bode plot C/f failed')

    plt.show()


    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',path2ZimT)
        clean_up_output('tVG',path2ZimT)
        print('Ouput data was deleted from '+path2ZimT)

    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds

    fig_idx = fig_idx + 1 

# reload the original device_parameters.txt file
os.rename(os.path.join(path2SIMsalabim,'device_parameters_saved.txt'), os.path.join(path2SIMsalabim,'device_parameters.txt'))
os.rename(os.path.join(path2ZimT,'device_parameters_saved.txt'), os.path.join(path2ZimT,'device_parameters.txt'))                



plt.show()