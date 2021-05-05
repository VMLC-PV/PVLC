###################################################
########## Make plots from Var_file ###############
###################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plot_settings_screen
from scipy import stats
import subprocess,shutil,os
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,zimt_tj_JV_plot,run_code
from tVG_gen import zimt_light_decay,zimt_voltage_decay,zimt_JV_sweep,zimt_TPV
import warnings

# Don't show warnings
warnings.filterwarnings("ignore")
color_nTL = '#c7e6a3'
color_pTL = '#8cb5f0'
color_pero = '#ba0000'
color_electrode ='#999999'
ext_save_pic = '.jpg'
Vt = 0.0256 # at 25 C

# Inputs
sytem = 'Windows'
size_fig = (16, 12)
num_fig_tjs= 0
num_fig_JV = 1

# zimt_light_decay(1e-10,20e-6,1.1e28,1e28,0,1e-8,tVG_name='tVG1.txt')
# zimt_light_decay(1e-10,20e-6,1.1e28,1e28,0,1e-7,tVG_name='tVG2.txt')

# zimt_voltage_decay(1e-10,10e-6,0.5,0.1,0,1e-9,tVG_name='tVG1.txt')
# zimt_voltage_decay(1e-10,10e-6,0.5,0.1,0,1e-8,tVG_name='tVG2.txt')

# zimt_JV_sweep(0,1.4,0.1,4.2e27,0.5,tVG_name='tVG1.txt')
# zimt_JV_sweep(0,1.4,0.1,4.2e27,0.1,tVG_name='tVG2.txt')
zimt_TPV(1e-8,1,1e30,1e-8,5e-8,time_exp =True,tVG_name='tVG1.txt')

Simu_str1 = '-tVG_file tVG1.txt' 
# Simu_str2 = '-tVG_file tVG2.txt' 

str_lst = [Simu_str1]#,Simu_str2]

labels = [sci_notation(1e-8,sig_fig=1)]#,sci_notation(1e-7,sig_fig=1)]

# colors = cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1))
colors = cm.rainbow(np.linspace(0,1,max(len(str_lst),4)+1))


run_simu = True
# JV_file plots
plot_tjs = True


idx = 0
curr_dir = os.getcwd()
path2ZimT = ''
for Simu_str,lbl in zip(str_lst,labels):
    # Run simulation
    if run_simu:
        run_zimt(Simu_str,sytem)
        
    # Import files
    path = ''
    JV_names = ['t', 'Va', 'Vdev', 'G', 'Rec', 'Jdev', 'range', 'Qnet', 'Jncat', 'Jpan']
    JV_file_name = r'\tj.dat'
    data_tj = pd.read_csv(curr_dir+path2ZimT+JV_file_name,delim_whitespace=True)
    

    ########################################################
    ################## JVs_file ############################
    ########################################################
    if plot_tjs:
        f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)
        zimt_tj_plot(num_fig_tjs,data_tj,0,lbl,colors[idx]) # input = num_fig,data_JV,plot_type,pic_save_name,save_yes,plot_jvexp,jvexp_name
        f_JVs = plt.figure(num_fig_JV,figsize=size_fig)
        zimt_tj_JV_plot(num_fig_JV,data_tj,0,lbl,colors[idx])
        

    
    # plt.semilogy(data_JV['Vext'],(data_JV['Jext'] + data_JV2['Jext'])/10)
    idx = idx+1

plt.show()

