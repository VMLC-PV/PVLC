###################################################
########## Make plots for SIMsalabim ###############
###################################################
# by Vincent M. Le Corre
# Package import
import subprocess,shutil,os,tqdm,parmap,multiprocessing,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# package by VLC
from VLC_useful_func import *
import plot_settings_screen


def run_plot_SIMsalabim():

    ## General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
            max_jobs = 1
            slash = '/'
            try:
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass
    else:
        slash = '/'

    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/DDSuite_v422/SIMsalabim'
    run_simu = True #Rerun simulation

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (10, 8)
    num_fig = 0
    # JV_file plots
    plot_JVs = True # Make JV plot
    plot_exp = False # Add experimental data to JV plot
    if plot_JVs:
        num_fig = num_fig + 1
        num_JV_plot = num_fig
        f_JVs = plt.figure(num_JV_plot,figsize=size_fig)
    # Var_file plots
    plot_nrj_diag = False # Make energy diagram plot
    plot_densities = False # Make density plot
    if plot_nrj_diag:
        num_fig = num_fig + 1
        num_nrj_diag_plot = num_fig
        f_nrj_diag = plt.figure(num_nrj_diag_plot,figsize=size_fig)

    if plot_densities:
            num_fig = num_fig + 1
            num_dens_plot = num_fig
            f_nrj_diag = plt.figure(num_dens_plot,figsize=size_fig)


    ## Prepare strings to run
    # Fixed string
    # fixed_str = '-Nc 3.946E+27 -mun_0 2.295E-7 -mup_0 2.476E-7 -W_L 4.110 -W_R 5.538 -Bulk_tr 3.431E+18 -Etrap 4.187 -kdirect 1.351E-17 -Gehp 1.210E+28 -Rseries 6.288E-5 -Rshunt 9.367' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    # fixed_str = '-Nc 1.212E+27 -mun_0 2.888E-8 -mup_0 3.054E-8 -W_L 4.3 -W_R 5.405 -Bulk_tr 0E+19 -Etrap 4.453 -kdirect 8.491E-18 -Rseries 1.759E-5 -Rshunt 7.622E-1 -Gehp 1.269E+28 -UseExpData 0 -Vmin -0.5 -CB 4.04 -VB 5.51 -VB_RTL 5.51 -Cn 1e-10 -Cp 1e-10 -Nc_RTL 1.212E+27 -Nc_LTL 1.212E+27 -eps_r_LTL 3.5 -eps_r_RTL 3.5 -nu_int_LTL 1e3 -nu_int_RTL 1e3 -until_Voc 1 -eps_r 25' # -Nc 1.212E+26 
    fixed_str = '-kdirect 10e-18 -mun_0 2e-8 -mup_0 16e-8'


    # Parameters to vary
    parameter1 = {'name':'L','values':[100e-9]}
    parameter2 = {'name':'L_LTL','values':[0e-9]}
    parameter3 = {'name':'L_RTL','values':[0e-9]}
    # parameter4 = {'name':'St_R','values':[0,1e10,1e12,1e13,1e14]}
    # parameter4 = {'name':'CB_LTL','values':[3.9,4.04,4.1,4.2,4.3]}

    L_LTL = parameter2['values'][0] # needed for nrj_diag plot
    L_RTL = parameter3['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter2,parameter3]#,parameter4] 
    # parameters = [parameter2 ,parameter5,parameter6] 

    
    str_lst,labels,JVexp_lst,JV_files,Var_files,sys_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]
    # JVexp_lst= ['PM6_3918_dark.txt','PM6_3918_int3.txt','PM6_3918_int10.txt','PM6_3918_int33.txt','PM6_3918_int100.txt','PM6_3918_int330.txt','PM6_3918_am15_long.txt','PM6_3918_int800.txt']

    for param in parameters:
        val.append(param['values'])
        nam.append(param['name'])

        if param['name'] == 'L_LTL':
            L_LTL = param['values'][0] # needed for nrj_diag plot

        if param['name'] == 'L_RTL':
            L_RTL = param['values'][0]  # needed for nrj_diag plot
    
    idx = 0
    for i in list(itertools.product(*val)):
        str_line = ''
        lab = ''
        JV_name = 'JV'
        Var_name = 'Var'
        for j,name in zip(i,nam):
            str_line = str_line +'-'+name+' {:.2e} '.format(j)
            lab = lab+name+' {:.2e} '.format(j)
            JV_name = JV_name +'_'+name +'_{:.2e}'.format(j)
            Var_name = Var_name +'_'+ name +'_{:.2e}'.format(j)
        str_lst.append(fixed_str+ ' ' +str_line+ '-JV_file '+JV_name+ '.dat -Var_file '+Var_name+'.dat')# -ExpJV '+JVexp_lst[idx])
        JV_files.append(Path(path2SIMsalabim) / str(JV_name+ '.dat'))
        Var_files.append(Path(path2SIMsalabim) / str(Var_name+ '.dat'))
        labels.append(lab)
        JVexp_lst.append('')
        sys_lst.append(system)
        path_lst.append(path2SIMsalabim)
        idx = idx + 1
    # print(str_lst)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots
    colors = ['r']
    
    # Run simulation
    if run_simu:
        # Run SIMsalabim
        run_multiprocess_simu(run_SIMsalabim,max_jobs,str_lst,sys_lst,path_lst)


    idx = 0
    perf_exp,perf_simu = [],[]
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name in zip(str_lst,JVexp_lst,labels,JV_files,Var_files):

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(path2SIMsalabim  / exp_name,delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()


            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=plot_exp,data_JVexp=data_JVexp,xlimits=[-0.5,1.1],ylimits=[-25,5],save_yes=True)
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=2,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp)

        ## Plot Var_file
        if plot_nrj_diag or plot_densities:
            data_var = pd.read_csv(Var_file_name,delim_whitespace=True) # Load Var_file   

        # Energy diagram plot
        if plot_nrj_diag:
            SIMsalabim_nrj_diag(num_nrj_diag_plot ,data_var,L_LTL,L_RTL,legend=False,Background_color=False,no_axis=False,pic_save_name='Energy_diagram'+ext_save_pic)

        # Carrier density plot
        if plot_densities:
            # What do we plot?
            dens2plot = ['n','p']
            line_type = ['-','--','-.']    
            # control colorbar
            colorbar = 'lin'
            if idx == 0 and colorbar != 'None':
                colorbar_display = True
            else:
                colorbar_display = False 
            # add legend for the density type
            plt.figure(num_dens_plot)
            leg_line = []
            for d,l in zip(dens2plot,line_type):
                leg_line.append(mlines.Line2D([], [], color='k', marker='None',markersize=15, label=d,linestyle=l)) # Create a legend for the first line.
            first_legend = plt.legend(handles=leg_line, loc='upper right',frameon=False)
            plt.gca().add_artist(first_legend) # Add the legend manually to the current Axes.
            #list(np.linspace(np.min(data_var['Vext']),np.max(data_var['Vext']),5))
            # SIMsalabim_dens_plot(num_dens_plot,data_var,Vext=list(np.linspace(np.min(data_var['Vext']),np.max(data_var['Vext']),5)),y=dens2plot,line_type=line_type,plot_type=2,colors=colors[idx],labels=labels[idx],legend=False,colorbar_type = colorbar,colorbar_display=colorbar_display)
            SIMsalabim_dens_plot(num_dens_plot,data_var,Vext=[np.max(data_var['Vext'])],y=dens2plot,line_type=line_type,plot_type=2,colors=colors[idx],labels=labels[idx],legend=False,colorbar_type = colorbar,colorbar_display=False)
            
            
        idx = idx+1
  


if __name__ == '__main__':

    run_plot_SIMsalabim()
    plt.show()