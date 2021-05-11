###################################################
######### Test SIMsalabim versus scaps ############
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,subprocess,shutil,math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from matplotlib.lines import Line2D
from time import time
import plot_settings_screen
from scipy import stats
from VLC_useful_func import sci_notation, run_SIMsalabim, SIMsalabim_nrj_diag, SIMsalabim_JVs_plot, make_df_JV, SIMsalabim_dens_plot, make_df_Var, make_df_JV
import warnings
from scipy import constants
from pathlib import Path
# Don't show warnings
warnings.filterwarnings("ignore")

# General Inputs
system = platform.system()                  # Operating system
# Max number of parallel simulations (for number of CPU use: os.cpu_count() )
max_jobs = os.cpu_count()-2
if system == 'Windows':             # cannot easily do multiprocessing in Windows
    max_jobs = 1
    slash = '/'
    try:
        os.system('taskkill.exe /F /IM SIMsalabim.exe')
    except:
        pass
else:
    slash = '/'

                                    # Rerun simu?

# Initialize
curr_dir = os.getcwd()
idx = 0
lines = ['-', '--', '-.', ':']
sys_lst, path_lst = [], []
# Figures control
size_fig = (14, 6)
fig_idx = 0

# Inputs
curr_dir = os.getcwd()                      # Current working directory
path2SIMsalabim = 'Simulation_program/DDSuite_v400/SIMsalabim'+slash    # Path to SIMsalabim in curr_dir
ext_save_pic = '.jpg'


### Test generation profile
Test_gen_profile = True
if Test_gen_profile:
    # Test whether the generation profile is correclty inputed in SIMsalabim
    print('\n')
    print('Start the Test generation profile:')
    L = 300e-9 # m
    L_LTL = 20e-9 # m
    L_RTL = 20e-9 # m
    TLsAbsorbs  = [0,1] # only 0 or 1
    Store_Folder = Path(curr_dir+slash+'Test simulation/')
    gen_profiles = ['None','blue.txt','red.txt']
    labels = ['No profile', 'Blue profile', 'Red profile']
    JV_files = ['JV_no_profile.dat','JV_blue.dat','JV_red.dat']
    Var_files = ['Var_no_profile.dat','Var_blue.dat','Var_red.dat']

    f_gen_pro,ax = plt.subplots(1,3,num =fig_idx, figsize=size_fig)
    fig_idx = fig_idx + 1
    subplot_num = 0
    for TLsAbsorb in TLsAbsorbs:
        # f_gen_pro = plt.figure(fig_idx, figsize=size_fig)
                
        str_lst = ['-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -Gen_profile  '+gen_profiles[0],
        '-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -Gen_profile  '+gen_profiles[1],'-L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -TLsAbsorb '+str(TLsAbsorb)+' -OutputRatio 100 -grad 10 -Gen_profile  '+gen_profiles[2]]
        # Simulation input
        run_simu = True     
        start = time()
        count = 0
        for i in str_lst:
            str_lst[count] = str_lst[count] +' -JV_file '+ JV_files[count] +' -Var_file '+ Var_files[count]
            sys_lst.append(system)
            path_lst.append(curr_dir+slash+path2SIMsalabim)
            count = count + 1
        
        # Color range for plotting
        colors = plt.cm.viridis(np.linspace(0, 1, max(len(str_lst), 4) + 1))
        
        # Run SIMsalabim
        if run_simu:
            p = multiprocessing.Pool(max_jobs)
            results = parmap.starmap(run_SIMsalabim, list(
                zip(str_lst, sys_lst, path_lst)), pm_pool=p, pm_processes=max_jobs, pm_pbar=True)
            p.close()
            p.join()
        
        # Plot
        for gen,jv,var,color,lbl in zip(gen_profiles,JV_files,Var_files,colors,labels):
            if gen != 'None':
                data_var = make_df_Var(Path(curr_dir+slash+path2SIMsalabim) / var)   
                data_gen = pd.read_csv(Path(curr_dir+slash+path2SIMsalabim) / gen,delim_whitespace=True,names=['x','Gehp'])
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
                data_var = make_df_Var(Path(curr_dir+slash+path2SIMsalabim) / var)
                if TLsAbsorb == 0: # correct profile to account for the part that is cut by SIMsalabim
                    data_var =  data_var[data_var['x'] >= L_LTL]
                    data_var =  data_var[data_var['x'] <= L-L_RTL]

                ax[subplot_num].plot(data_var['x']/max(data_var['x']),data_var['Gehp']/data_var['Gehp'].mean(),color=color,label=lbl,linestyle='-')
            data_jv = make_df_JV(Path(curr_dir+slash+path2SIMsalabim) / jv)
            ax[2].plot(data_jv['Vext'],data_jv['Jphoto']/10,color=color,linestyle=lines[subplot_num],label='TLAbsorb = {:.0f}'.format(TLsAbsorb))
            print(data_jv['Jphoto'].mean())
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
    plt.savefig(Path(Store_Folder) / 'test_gen_pro.jpg',dpi=100,transparent=True)

plt.show()
    
    
