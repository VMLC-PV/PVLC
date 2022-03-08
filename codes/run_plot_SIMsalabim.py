#########################################################
########## Run and plot SIMsalabim output ###############
#########################################################
# by Vincent M. Le Corre
# Package import
import os,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# Import homemade package by VLC
from VLC_units.plots.SimSS_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
import VLC_units.plots.plot_settings_screen



def run_plot_SIMsalabim(fixed_str = None, parameters = None, path2SIMsalabim = None, run_simu = True, plot_JVs = True, plot_nrj_diag = True, plot_densities = True, plot_exp = False, JVexp_lst = None,  verbose = True):
    """Run and plot SIMsalabim outputs.

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    parameters : dict, optional
        Dictionary with the input for the SIMsalabim function, we will simulate all permutations, by default None.
    path2SIMsalabim : str, optional
        Path to Simss, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_JVs : bool, optional
        Make JV plot?, by default True
    plot_nrj_diag : bool, optional
        Make energy diagram plot, by default True
    plot_densities : bool, optional
        Make density plot, by default True
    plot_exp : bool, optional
        Add experimental data to JV plot, by default False
    JVexp_lst : list, optional
        List of experimental JV data filenames in the path2SIMsalabim directory, by default None
    verbose : bool, optional
        Verbose?, by default True
    
    Returns
    -------
    JV_files,Var_files,scPars_files : list,list
        Lists of JV and variable and performance filenames
    """
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


    curr_dir = os.getcwd()              # Current working directory
    if path2SIMsalabim  is None:
        path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_guide/SimSS')                    # Path to SimSS
                                                                                            # Rerun simulation

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ''  # add any fixed string to the simulation command

    # Parameters to vary
    if parameters is None:
        parameters = []
        parameters.append({'name':'Gfrac','values':[1]})

    ParFileDic = ReadParameterFile(f"{path2SIMsalabim}/device_parameters.txt") # Read device parameters

    # Initialize     
    str_lst,labels,JV_files,Var_files,scPars_files,code_name_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]
    if JVexp_lst is None:
        JVexp_lst = []
        add_exp = False
        # JVexp_lst= ['PM6_3918_dark.txt','PM6_3918_int3.txt','PM6_3918_int10.txt','PM6_3918_int33.txt','PM6_3918_int100.txt','PM6_3918_int330.txt','PM6_3918_am15_long.txt','PM6_3918_int800.txt']
    else:
        add_exp = True

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (12,7.5)
    num_fig = 0
    # JV_file plots
    if plot_JVs:
        num_fig += 1
        num_JV_plot = num_fig
        f_JVs = plt.figure(num_JV_plot,figsize=size_fig)
    # Var_file plots
    if plot_nrj_diag:
        num_fig += 1
        num_nrj_diag_plot = num_fig
        f_nrj_diag = plt.figure(num_nrj_diag_plot,figsize=size_fig)

    if plot_densities:
            num_fig += 1
            num_dens_plot = num_fig
            f_nrj_diag = plt.figure(num_dens_plot,figsize=size_fig)
    
    for param in parameters: # initalize lists of values and names
        val.append(param['values'])
        nam.append(param['name'])

    
    idx = 0
    for i in list(itertools.product(*val)): # iterate over all combinations of parameters
        str_line = ''
        lab = ''
        JV_name = 'JV'
        Var_name = 'Var'
        scPars_name = 'scPars'
        for j,name in zip(i,nam):
            str_line = str_line +'-'+name+' {:.2e} '.format(j)
            lab = lab+name+' {:.2e} '.format(j)
            JV_name = JV_name +'_'+name +'_{:.2e}'.format(j)
            Var_name = Var_name +'_'+ name +'_{:.2e}'.format(j)
            scPars_name = scPars_name +'_'+ name +'_{:.2e}'.format(j)
        str_lst.append(fixed_str+ ' ' +str_line+ '-JV_file '+JV_name+ '.dat -Var_file '+Var_name+'.dat -scPars_file '+scPars_name+'.dat')# -ExpJV '+JVexp_lst[idx])
        JV_files.append(os.path.join(path2SIMsalabim , str(JV_name+ '.dat')))
        Var_files.append(os.path.join(path2SIMsalabim , str(Var_name+ '.dat')))
        scPars_files.append(os.path.join(path2SIMsalabim , str(scPars_name+ '.dat')))
        labels.append(lab)
        if not add_exp:
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
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name,scPars_file_name in zip(str_lst,JVexp_lst,labels,JV_files,Var_files,scPars_files):

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(os.path.join(path2SIMsalabim, exp_name),delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()


            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=plot_exp,data_JVexp=data_JVexp,xlimits=[-0.5,1.1],ylimits=[-25,10],save_yes=True,pic_save_name=os.path.join(path2SIMsalabim , str('JV'+ext_save_pic)))
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=2,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=False,data_JVexp=data_JVexp,xlimits=[-0.5,2],ylimits=[1e-4,1e4],save_yes=True,pic_save_name=os.path.join(path2SIMsalabim , str('JV'+ext_save_pic)))

        ## Plot Var_file
        if plot_nrj_diag or plot_densities:
            data_var = pd.read_csv(Var_file_name,delim_whitespace=True) # Load Var_file   

        # Energy diagram plot
        if plot_nrj_diag:
            L_LTL = float(ChosePar('L_LTL',GetParFromStr(Simu_str),ParFileDic)) # Needed for nrj_diag plot
            L_RTL =  float(ChosePar('L_RTL',GetParFromStr(Simu_str),ParFileDic)) # Needed for nrj_diag plot
            SIMsalabim_nrj_diag(num_nrj_diag_plot ,data_var,L_LTL,L_RTL,legend=False,Background_color=True,no_axis=False,pic_save_name=os.path.join(path2SIMsalabim,'Energy_diagram'+ext_save_pic))

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
    
    return JV_files,Var_files,scPars_files
  


if __name__ == '__main__':
    
    parameters = []
    parameters.append({'name':'Gfrac','values':list(np.geomspace(1e-3,5,15))})
    # parameters.append({'name':'Rshunt','values':[-7.995E-1,7.995E-1,7.995E-2,7.995E-3]})
    # parameters.append({'name':'Rseries','values':[0,3.063E-4,3.063E-3,3.063E-2]})
    # fixed_strs = ['-Rshunt -7.995E-1','-Rshunt 7.995E-1','-Rshunt 7.995E-2','-Rshunt 7.995E-3']
    fixed_strs = ['-Rseries 0','-Rseries 3.063E-4','-Rseries 3.063E-3','-Rseries 3.063E-2']
    colors = plt.cm.viridis(np.linspace(0,1,max(len(fixed_strs),4)+1))
    count = 0
    for fix in fixed_strs:
        JV_files,Var_files,scPars_files=run_plot_SIMsalabim(fixed_str=fix,parameters=parameters,plot_densities=False, plot_nrj_diag=False)

        Voc = []

        for i in scPars_files:
                perf = pd.read_csv(i,delim_whitespace=True)
                Voc.append(perf['Voc'])
            
        plt.figure(12,figsize=(12,7.5))
        plt.semilogx(list(np.geomspace(1e-3,5,15)),Voc,color=colors[count])
        count += 1
    plt.xlabel('Intensity [suns]')
    plt.ylabel('V$_{OC}$ [V]')

    # plt.ylabel('FF')
    plt.xlim([5e-4,5])
    # plt.ylim([0.5,0.9])
    plt.ylim([0.35,1.05])
    plt.tight_layout()

    plt.savefig(os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_guide/SimSS','Voc_Gfrac.png'))

    plt.show()