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

    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/AutoFit125/DDSuite_v407_Xiaoyan/SIMsalabim'
    run_simu = True #Rerun simulation

    ## Figures control
    ext_save_pic = '.jpg'
    size_fig = (10, 8)
    num_fig = 0
    # JV_file plots
    plot_JVs = True # Make JV plot
    plot_exp = True # Add experimental data to JV plot
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
    
    ##3918
    ID = 3918
    fixed_str = '-Nc 1.342E+27 -mun_0 3.685E-7 -mup_0 2.312E-7 -W_L 4.110 -W_R 5.476 -Bulk_tr 5.638E+19 -Etrap 4.403 -kdirect 6.620E-18 -Rseries 4.297E-5 -Rshunt 5.873E-1 -Gehp 1.220E+28 -UseExpData 0 -Vmin -0.5' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

 
    ## 3936
    # ID = 3936
    # fixed_str = '-Nc 2.986E+27 -mun_0 1.157E-7 -mup_0 2.052E-7 -W_L 4.226 -W_R 5.511 -Bulk_tr 3.175E+18 -Etrap 4.435 -kdirect 4.500E-17 -Rseries 4.251E-5 -Rshunt 8.677E-1 -Gehp 1.281E+28 -UseExpData 0 -Vmin -0.5' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    
    
    # ID = 3937
    # fixed_str = '-Nc 1.123E+27 -mun_0 3.624E-7 -mup_0 2.244E-7 -W_L 4.110 -W_R 5.524 -Bulk_tr 8.243E+20 -Etrap 4.339 -kdirect 2.491E-18 -Rseries 2.345E-5 -Rshunt 1.038 -Gehp 1.340E+28 -UseExpData 0 -Vmin -0.5' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    # Parameters to vary
    parameter1 = {'name':'L','values':[140e-9]}
    parameter2 = {'name':'Gfrac','values':[0,0.00877, 0.0295, 0.0967, 0.2903605440421194, 0.943, 1, 2.258]}
    parameter3 = {'name':'L_LTL','values':[30e-9]}
    parameter4 = {'name':'L_RTL','values':[10e-9]}
    L_LTL = parameter3['values'][0] # needed for nrj_diag plot
    L_RTL = parameter4['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter2,parameter3,parameter4] 

    
    str_lst,labels,JVexp_lst,JV_files,Var_files,sys_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]
    # JVexp_lst= ['PM6_3918_dark.txt','PM6_3918_int3.txt','PM6_3918_int10.txt','PM6_3918_int33.txt','PM6_3918_int100.txt','PM6_3918_int330.txt','PM6_3918_am15_long.txt','PM6_3918_int800.txt']
    JVexp_lst= ['PM6_'+str(ID)+'_dark.txt','PM6_'+str(ID)+'_int3.txt','PM6_'+str(ID)+'_int10.txt','PM6_'+str(ID)+'_int33.txt','PM6_'+str(ID)+'_int100.txt','PM6_'+str(ID)+'_int330.txt','PM6_'+str(ID)+'_am15.txt','PM6_'+str(ID)+'_int800.txt']

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
        str_lst.append(fixed_str+ ' ' +str_line+ '-JV_file '+JV_name+ '.dat -Var_file '+Var_name+'.dat -ExpJV '+JVexp_lst[idx])
        JV_files.append(Path(path2SIMsalabim) / str(JV_name+ '.dat'))
        Var_files.append(Path(path2SIMsalabim) / str(Var_name+ '.dat'))
        labels.append(lab)
        # JVexp_lst.append('')
        sys_lst.append(system)
        path_lst.append(path2SIMsalabim)
        idx = idx + 1
    # print(str_lst)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots

    
    # Run simulation
    if run_simu:
        # Run SIMsalabim
        run_multiprocess_simu(run_SIMsalabim,max_jobs,str_lst,sys_lst,path_lst)


    idx = 0
    perf_exp,perf_simu = [],[]
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name,Gfrac in zip(str_lst,JVexp_lst,labels,JV_files,Var_files,parameter2['values']):

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(path2SIMsalabim  / exp_name,delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()

            if Gfrac > 0:
                perf_simu.append([get_Voc(data_JV['Vext'],data_JV['Jext']),get_Jsc(data_JV['Vext'],data_JV['Jext']),get_FF(data_JV['Vext'],data_JV['Jext']),get_PCE(data_JV['Vext'],data_JV['Jext'],suns=Gfrac)/10])
                perf_exp.append([get_Voc(data_JVexp['V'],data_JVexp['J']),get_Jsc(data_JVexp['V'],data_JVexp['J']),get_FF(data_JVexp['V'],data_JVexp['J']),get_PCE(data_JVexp['V'],data_JVexp['J'],suns=Gfrac)/10])
            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp,xlimits=[-0.5,1.1],ylimits=[-50,5],save_yes=True,pic_save_name=path2SIMsalabim/'JVfit.jpg')
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
            line_type = ['-','--']    
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

            SIMsalabim_dens_plot(num_dens_plot,data_var,Vext=list(np.linspace(np.min(data_var['Vext']),np.max(data_var['Vext']),5)),y=dens2plot,line_type=line_type,plot_type=2,colors=colors[idx],labels=labels[idx],legend=False,colorbar_type = colorbar,colorbar_display=colorbar_display)
            
            
        idx = idx+1
    perf_exp = np.asarray(perf_exp)
    perf_simu = np.asarray(perf_simu)

    lights = parameter2['values'].copy()
    lights.pop(0)

    plt.figure(2,figsize=size_fig)
    plt.subplot(221)
    plt.semilogx(lights,perf_exp[:,3],linestyle='None',marker='o',markersize=10,markerfacecolor='w',label='Exp')
    plt.semilogx(lights,perf_simu[:,3],linestyle='-',marker='None',markersize=10,markerfacecolor='w',label='Simu')
    plt.legend(loc='best',frameon=False)
    plt.xlabel('Suns')
    plt.ylabel('PCE [%]')
    # plt.ylim([8,15])
    plt.grid(b=True,which='both')
    plt.subplot(222)
    plt.semilogx(lights,perf_exp[:,2],linestyle='None',marker='o',markersize=10,markerfacecolor='w')
    plt.semilogx(lights,perf_simu[:,2],linestyle='-',marker='None',markersize=10,markerfacecolor='w')
    plt.xlabel('Suns')
    plt.ylabel('FF')
    # plt.ylim([0.5,1])
    plt.grid(b=True,which='both')
    plt.subplot(223)
    plt.semilogx(lights,perf_exp[:,0],linestyle='None',marker='o',markersize=10,markerfacecolor='w')
    plt.semilogx(lights,perf_simu[:,0],linestyle='-',marker='None',markersize=10,markerfacecolor='w')
    plt.xlabel('Suns')
    plt.ylabel('V$_{OC}$ [V]')
    # plt.ylim([0.6,0.9])
    plt.grid(b=True,which='both')
    plt.subplot(224)
    plt.loglog(lights,-perf_exp[:,1]/10,linestyle='None',marker='o',markersize=10,markerfacecolor='w')
    plt.loglog(lights,-perf_simu[:,1]/10,linestyle='-',marker='None',markersize=10,markerfacecolor='w')
    plt.xlabel('Suns')
    plt.ylabel('J$_{SC}$ mA cm$^{-2}$')
    # plt.ylim([0,26])
    plt.tight_layout()
    plt.grid(b=True,which='both')
    plt.savefig(path2SIMsalabim/'perf_fit.jpg')


if __name__ == '__main__':

    run_plot_SIMsalabim()
    plt.show()