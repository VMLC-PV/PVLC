###################################################
########## Make plots from Var_file ###############
###################################################
# by Vincent M. Le Corre
# Package import
import subprocess,shutil,os,tqdm,parmap,multiprocessing,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# Package by VLC
from VLC_useful_func import *
from tVG_gen import *
import plot_settings_screen

def run_plot_zimt():
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

    curr_dir = os.getcwd()
    path2ZimT = Path(os.getcwd()) /'Simulation_program/DDSuite_v405_peroPL/ZimT'
    run_simu = True #Rerun simulation

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (16, 12)
    num_fig = 0
    # tj_file plots
    plot_tjs = True # Make JV plot
    if plot_tjs:
        num_fig = num_fig + 1
        num_tj_plot = num_fig
        f_tjs = plt.figure(num_tj_plot,figsize=size_fig)
    # Var_file plots
    plot_nrj_diag = False # Make energy diagram plot
    plot_densities = True # Make density plot
    if plot_nrj_diag:
        num_fig = num_fig + 1
        num_nrj_diag_plot = num_fig
        f_nrj_diag = plt.figure(num_nrj_diag_plot,figsize=size_fig)

    if plot_densities:
            num_fig = num_fig + 1
            num_dens_plot = num_fig
            f_nrj_diag = plt.figure(num_dens_plot,figsize=size_fig)


    ## Prepare tVG file (Can be any of the experiment described in tVG_gen.py)
    name_tVG = 'tVG_PL.txt'
    # zimt_TPC(1e-8,1e-6,1e31,0,1e-8,tpulse=5e-8,time_exp =True,tVG_name=path2ZimT/name_tVG) # Simulation input for zimt_TPC (see tVG_gen.py) 
    # zimt_voltage_step(1e-4,3,1.15,0,3.35e27,steps=100,trf = 10e-9,time_exp =True,tVG_name=path2ZimT/name_tVG)
    zimt_TrPL(1e-9,1e-6,1e20,1e15,0,1e-8,tpulse=6e-9,time_exp =True,tVG_name=path2ZimT/name_tVG)
    ## Prepare strings to run
    # Fixed string
    fixed_str = '-tVG_file '+name_tVG
    # Parameters to vary
    parameter1 = {'name':'L','values':[400e-9]}
    # parameter2 = {'name':'Lang_pre','values':[0.1]}
    parameter3 = {'name':'L_LTL','values':[0e-9]}
    parameter4 = {'name':'L_RTL','values':[0e-9]}
    th_TL_left = parameter3['values'][0] # needed for nrj_diag plot
    th_TL_right = parameter4['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter3,parameter4] 

    str_lst,labels,JVexp_lst,tj_files,Var_files,sys_lst,path_lst = [],[],[],[],[],[],[]
    val,nam = [],[]
    for param in parameters:
        val.append(param['values'])
        nam.append(param['name'])

    for i in list(itertools.product(*val)):
        str_line = fixed_str+ ' '
        lab = ''
        tj_name = 'tj'
        Var_name = 'Var'
        for j,name in zip(i,nam):
            str_line = str_line +'-'+name+' {:.2e} '.format(j)
            lab = lab+name+' {:.2e} '.format(j)
            tj_name = tj_name +'_'+name +'_{:.2e}'.format(j)
            Var_name = Var_name +'_'+ name +'_{:.2e}'.format(j)
        str_lst.append(str_line + '-tj_file '+tj_name+ '.dat -Var_file '+Var_name+'.dat')
        tj_files.append(Path(path2ZimT) / str(tj_name+ '.dat'))
        Var_files.append(Path(path2ZimT) / str(Var_name+ '.dat'))
        labels.append(lab)
        JVexp_lst.append('')
        sys_lst.append(system)
        path_lst.append(path2ZimT)
    print(str_lst)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots
    

    ## Run simulation
    if run_simu:
        # Run SIMsalabim
        run_multiprocess_simu(run_zimt,max_jobs,str_lst,sys_lst,path_lst)


    ## Plots
    idx = 0
    for Simu_str,exp_name,lbl,tj_file_name,Var_file_name in zip(str_lst,JVexp_lst,labels,tj_files,Var_files):

        ## Plot tjs
        if plot_tjs:
            data_tj = pd.read_csv(path2ZimT/tj_file_name,delim_whitespace=True)
            zimt_tj_plot(num_tj_plot,data_tj,x='t',y=['Jdir'],xlimits=[],ylimits=[],plot_type=2,labels=labels[idx],colors=colors[idx],line_type = ['-'],mark='',legend=False,save_yes=False,pic_save_name='transient.jpg')
            # zimt_tj_plot(num_tj_plot,data_tj,x='t',y=['Jdir'],xlimits=[],ylimits=[],plot_type=2,labels=labels[idx],colors=colors[idx],line_type = ['-'],mark='',legend=False,save_yes=False,pic_save_name='transient.jpg')  
  


        ## Plot Var_file
        if plot_nrj_diag or plot_densities:
            data_var = pd.read_csv(Var_file_name,delim_whitespace=True) # Load Var_file
            
        # Energy diagram plot
        # if plot_nrj_diag:
        #     SIMsalabim_nrj_diag(num_nrj_diag_plot ,data_var,th_TL_left,th_TL_right,legend=False,Background_color=False,no_axis=False,pic_save_name='Energy_diagram'+ext_save_pic)

        # Carrier density plot
        if plot_densities:
             # What do we plot?
            dens2plot = ['n','p']
            line_type = ['-','--']    
            # control colorbar
            colorbar = 'log'
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

            zimt_dens_plot(num_dens_plot,data_var,time=list(np.geomspace(1e-8,1e-6,10)),y=['n','p'],plot_type=2,colors=colors[idx],labels=labels[idx],legend=False,colorbar_type = colorbar,colorbar_display=colorbar_display)
            
        idx = idx+1



if __name__ == '__main__':

    run_plot_zimt()
    plt.show()