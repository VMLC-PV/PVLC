###################################################
########## Make plots for SIMsalabim ###############
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
from VLC_units.plots.SimSS_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *



def run_plot_SIMsalabim():

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


    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/SIMsalabim_v425/SimSS'
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
    plot_nrj_diag = True # Make energy diagram plot
    plot_densities = True # Make density plot
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
    fixed_str = ''


    # Parameters to vary
    parameters = []
    parameters.append({'name':'L','values':[140e-9]})
    parameters.append({'name':'L_LTL','values':[30e-9]})
    parameters.append({'name':'L_RTL','values':[10e-9]})
    parameters.append({'name':'W_L','values':[4.1,4.2]})

    ParFileDic = ReadParameterFile(f"{path2SIMsalabim}/device_parameters.txt") # Read device parameters

    
    str_lst,labels,JVexp_lst,JV_files,Var_files,code_name_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]
    # JVexp_lst= ['PM6_3918_dark.txt','PM6_3918_int3.txt','PM6_3918_int10.txt','PM6_3918_int33.txt','PM6_3918_int100.txt','PM6_3918_int330.txt','PM6_3918_am15_long.txt','PM6_3918_int800.txt']

    for param in parameters: # initalize lists of values and names
        val.append(param['values'])
        nam.append(param['name'])

    
    idx = 0
    for i in list(itertools.product(*val)): # iterate over all combinations of parameters
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
        code_name_lst.append('SimSS')
        path_lst.append(path2SIMsalabim)
        idx += 1
    # print(str_lst)

    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots

    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)



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


            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=plot_exp,data_JVexp=data_JVexp,xlimits=[-0.5,1.1],ylimits=[-25,5])
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=2,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp)

        ## Plot Var_file
        if plot_nrj_diag or plot_densities:
            data_var = pd.read_csv(Var_file_name,delim_whitespace=True) # Load Var_file   

        # Energy diagram plot
        if plot_nrj_diag:
            L_LTL = ChosePar('L_LTL',GetParFromStr(Simu_str),ParFileDic) # Needed for nrj_diag plot
            L_RTL = ChosePar('L_RTL',GetParFromStr(Simu_str),ParFileDic) # Needed for nrj_diag plot
            SIMsalabim_nrj_diag(num_nrj_diag_plot ,data_var,L_LTL,L_RTL,legend=False,Background_color=True,no_axis=False,pic_save_name='Energy_diagram'+ext_save_pic)

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
            
            
        idx += 1
  


if __name__ == '__main__':
    
    run_plot_SIMsalabim()
    plt.show()