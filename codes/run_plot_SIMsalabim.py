###################################################
########## Make plots from Var_file ###############
###################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_settings_screen
from scipy import stats
import subprocess,shutil,os
from VLC_useful_func import sci_notation,run_SIMsalabim,SIMsalabim_nrj_diag,SIMsalabim_JVs_plot,make_df_JV
import warnings
from scipy import constants
# Don't show warnings
warnings.filterwarnings("ignore")
color_nTL = '#c7e6a3'
color_pTL = '#8cb5f0'
color_pero = '#ba0000'
color_electrode ='#999999'
ext_save_pic = '.jpg'
Vt = 0.0256 # at 25 C
kb = constants.value(u'Boltzmann constant')
k = constants.value(u'Boltzmann constant in eV/K')
h = constants.value(u'Planck constant')
me = constants.value(u'electron mass')

## Inputs
sytem = 'Windows'
size_fig = (16, 12)
num_fig = 0
path2SIMsalabim = os.getcwd()


## Prepare strings to run
Ls = [350e-9,450e-9]
str_lst,labels,JVexp_lst,JV_files,Var_files = [],[],[],[],[]

for i in Ls:
    str_lst.append('-L {:.2e} -JV_file JV{:.2e}nm.dat -Var_file Var{:.2e}nm.dat'.format(i,i,i))
    labels.append('{:.0f} nm'.format(i*1e9))
    JV_files.append(r'\JV{:.2e}nm.dat'.format(i))
    Var_files.append(r'\Var{:.2e}nm.dat'.format(i))
    JVexp_lst.append('')


colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1))
th_TL_left = 30 # nm
th_TL_right = 10 # nm

run_simu = False
# JV_file plots
plot_JVs,plot_exp = True,False

# Var_file plots
plot_nrj_diag = False

idx = 0

for Simu_str,exp_name,lbl,JV_file_name,Var_file_name in zip(str_lst,JVexp_lst,labels,JV_files,Var_files):
    
    # Run simulation
    if run_simu:
        run_SIMsalabim(Simu_str,sytem)

    # Import files
    path = os.getcwd()+''
    JV_names = ['Vext','Vint','Jext','Jint','P','recLan','recSRH','Jbimo','JSRH_bulk','JSRH_LI','JSRH_RI','Jph','Jn_l','Jp_l','Jn_r','Jp_r']
    
    data_JV = make_df_JV(path+JV_file_name)
    if plot_exp:
        data_JVexp = pd.read_csv(exp_name,delim_whitespace=True)
    else:
        data_JVexp = pd.DataFrame()
    
    ########################################################
    ################## JVs_file ############################
    ########################################################
    if plot_JVs:
        num_fig = 1
        num_fig_JVs = num_fig
        f_JVs = plt.figure(num_fig_JVs,figsize=size_fig)
        SIMsalabim_JVs_plot(num_fig_JVs,data_JV,[],[],0,lbl,colors[idx],plot_exp,data_JVexp,True,'JVs_plot'+ext_save_pic) 


    ########################################################
    ################## Var_file ############################
    ########################################################
    ## Energy diagram plot
    if plot_nrj_diag:
        Var_names = ['x','V','n','p','Evac','Ec','Ev','phin','phip','ntrap','ptrap','nid','pid','nion','pion','mun','mup','rec','dp','Gm','Jn','Jp']
        Var_file_name = Var_file_name
        data_var = pd.read_csv(path+Var_file_name,delim_whitespace=True)
        data_var['x'] =data_var['x'] * 1e9
   
        num_fig = 2
        num_fig_nrj_diag = num_fig
        f_nrj_diag = plt.figure(num_fig_nrj_diag,figsize=size_fig)
        SIMsalabim_nrj_diag(num_fig_nrj_diag,data_var,th_TL_left,th_TL_right,'Energy_diagram'+ext_save_pic,True)
        

    idx = idx+1

plt.show()