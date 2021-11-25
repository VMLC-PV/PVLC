###################################################
########## Make plots for ZimT ###############
###################################################
# by Vincent M. Le Corre
# Package import
import os,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# package by VLC
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *

def run_plot_zimt():

    ## General Inputs
    warnings.filterwarnings("ignore")           # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                 # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    do_multiprocessing = False                      # Use multiprocessing
    if system == 'Windows':                     # cannot easily do multiprocessing in Windows
            max_jobs = 1
            do_multiprocessing = False
            try:                                # kill all running jobs to avoid conflicts
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass

    curr_dir = os.getcwd()
    path2ZimT = Path(os.getcwd()) /'Simulation_program/SIMsalabim_v425/ZimT'
    run_simu = True #Rerun simulation

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (10, 8)
    num_fig = 0
    # tj_file plots
    plot_tjs = True # Make tj plot
    if plot_tjs:
        num_fig += 1
        num_tj_plot = num_fig
        f_tjs = plt.figure(num_tj_plot,figsize=size_fig)
    # Var_file plots
    plot_nrj_diag = False # Make energy diagram plot
    plot_densities = True # Make density plot
    if plot_nrj_diag:
        num_fig += 1
        num_nrj_diag_plot = num_fig
        f_nrj_diag = plt.figure(num_nrj_diag_plot,figsize=size_fig)

    if plot_densities:
            num_fig += 1
            num_dens_plot = num_fig
            f_nrj_diag = plt.figure(num_dens_plot,figsize=size_fig)


    ## Prepare tVG file (Can be any of the experiment described in tVG_gen.py)
    name_tVG = 'tVG.txt'
    # zimt_TPC(1e-8,1e-6,1e28,0,1e-8,5e-8,time_exp =True,tVG_name=path2ZimT/name_tVG) # Simulation input for zimt_TPC (see tVG_gen.py) 
    zimt_CELIV(1e-8,10e-6,0,-1/1e-6,1e22,5e-8,1e-7,0,width_pulse = 6e-9,time_exp=True,steps=100,tVG_name=path2ZimT/name_tVG) # Simulation input for zimt_CELIV (see tVG_gen.py)
    # zimt_voltage_step(1e-4,3,1.15,0,3.35e27,steps=100,trf = 10e-9,time_exp =True,tVG_name=path2ZimT/name_tVG)
    # zimt_TrPL(1e-9,1e-6,1e20,1e15,0,1e-8,tpulse=6e-9,time_exp =True,tVG_name=path2ZimT/name_tVG)

    ## Prepare strings to run
    # Fixed string
    fixed_str = '-tVG_file '+name_tVG
    # Parameters to vary
    parameters = []
    parameters.append({'name':'L','values':[140e-9]})
    # parameters.append({'name':'L_LTL','values':[30e-9]})
    # parameters.append({'name':'L_RTL','values':[10e-9]})
    # parameters.append({'name':'W_L','values':[4.1]})

    ParFileDic = ReadParameterFile(f"{path2ZimT}/device_parameters.txt") # Read device parameters


    str_lst,labels,JVexp_lst,tj_files,Var_files,code_name_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]

    for param in parameters: # initalize lists of values and names
        val.append(param['values'])
        nam.append(param['name'])

    
    idx = 0
    for i in list(itertools.product(*val)): # iterate over all combinations of parameters
        str_line = ''
        lab = ''
        tj_name = 'tj'
        Var_name = 'Var'
        for j,name in zip(i,nam):
            str_line = str_line +'-'+name+' {:.2e} '.format(j)
            lab = lab+name+' {:.2e} '.format(j)
            tj_name = tj_name +'_'+name +'_{:.2e}'.format(j)
            Var_name = Var_name +'_'+ name +'_{:.2e}'.format(j)
        str_lst.append(fixed_str+ ' ' +str_line+ '-tj_file '+tj_name+ '.dat -Var_file '+Var_name+'.dat')
        tj_files.append(Path(path2ZimT) / str(tj_name+ '.dat'))
        Var_files.append(Path(path2ZimT) / str(Var_name+ '.dat'))
        labels.append(lab)
        code_name_lst.append('zimt')
        path_lst.append(path2ZimT)
        idx += 1
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots


    ## Run simulation
    if run_simu: 
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)


    ## Plots
    idx = 0
    for Simu_str,lbl,tj_file_name,Var_file_name in zip(str_lst,labels,tj_files,Var_files):

        ## Plot tjs
        if plot_tjs:
            data_tj = pd.read_csv(path2ZimT/tj_file_name,delim_whitespace=True)

            zimt_tj_plot(num_tj_plot,data_tj,x='t',y=['Jext'],xlimits=[],ylimits=[],plot_type=1,labels=labels[idx],colors=colors[idx],line_type = ['-'],mark='',legend=False,save_yes=False,pic_save_name='transient.jpg')
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
            
        idx += 1



if __name__ == '__main__':

    run_plot_zimt()
    plt.show()