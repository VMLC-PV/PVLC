############################################################
########## Degradation simulation with SIMsalabim ##########
############################################################
# by Vincent M. Le Corre
# Package import
import subprocess,shutil,os,tqdm,parmap,multiprocessing,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.optimize
from pathlib import Path
# package by VLC
from VLC_useful_func import *
from tVG_gen import *
import plot_settings_screen


def run_plot_SIMsalabim(a,b):

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
    path2ZimT = Path(os.getcwd()) /'Simulation_program/DDSuite_v405/ZimT'
    run_simu_JV = True #Rerun simulation
    run_simu_CELIV = False #Rerun simulation

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
    
    # CELIV plots 
    plot_CELIVs = False
    if plot_JVs:
        num_fig = num_fig + 1
        num_CELIV_plot = num_fig
        f_CELIVs = plt.figure(num_CELIV_plot,figsize=size_fig)


    ## Prepare strings to run
    # Fixed string
    fixed_str = '-Nc 2.986E+27 -mup_0 2.052E-7 -W_L 4.226 -W_R 5.511 -Bulk_tr 3.175E+18 -Etrap 4.435 -kdirect 4.500E-17 -Rseries 4.251E-5 -Rshunt 8.677E-1 -Gehp 1.281E+28 ' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    # Parameters to vary
    # time = np.arange(0,100,10)
    time = np.geomspace(1e-1,1000,13)
    time = np.insert(time,0,0)
    # a = MonoExpDecay(time, 1, 1e-7, 1e-10)
    # b = StretchedExp(time, 1, 3/12, 1e-8, 1e-9)
    mob = a + b
    lang = MonoExpInc(time, 10, 0.001, 0.002,)
    parameter1 = {'name':'L','values':[140e-9]}
    parameter2 = {'name':'mun_0','values':list(mob)}
    # parameter2 = {'name':'Lang_pre','values':list(lang)}
    parameter3 = {'name':'L_LTL','values':[20e-9]}
    parameter4 = {'name':'L_RTL','values':[30e-9]}
    L_LTL = parameter3['values'][0] # needed for nrj_diag plot
    L_RTL = parameter4['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter2,parameter3,parameter4] 

    str_lst,labels,JVexp_lst,JV_files,Var_files,scPars_files,sys_lst,path_lst,val,nam,tVG_files_CELIV,tj_files_CELIV,path_lst_CELIV,str_lst_CELIV = [],[],[],[],[],[],[],[],[],[],[],[],[],[]

    for param in parameters:
        val.append(param['values'])
        nam.append(param['name'])

        if param['name'] == 'L_LTL':
            L_LTL = param['values'][0] # needed for nrj_diag plot

        if param['name'] == 'L_RTL':
            L_RTL = param['values'][0]  # needed for nrj_diag plot
    
    ## CELIV parameters
    tmin = 1e-8
    tmax = 50e-6
    Voffset = 1
    slope = -2/60e-6
    Gen = 1e21                                              # Max generation rate for the gaussian laser pulse
    tpulse = 5e-8
    tstep = 1e-7
    tdelay = 0
    
    
    for i in list(itertools.product(*val)):
        str_line = fixed_str
        lab = ''
        JV_name = 'JV'
        Var_name = 'Var'
        tj_name = 'tj_CELIV'
        tVG_name = 'tVG_CELIV'
        scPars_name = 'scPars'
        for j,name in zip(i,nam):
            str_line = str_line +'-'+name+' {:.2e} '.format(j)
            lab = lab+name+' {:.2e} '.format(j)
            JV_name = JV_name +'_'+name +'_{:.2e}'.format(j)
            Var_name = Var_name +'_'+ name +'_{:.2e}'.format(j)
            scPars_name = scPars_name +'_'+ name +'_{:.2e}'.format(j)
            tj_name = tj_name +'_'+ name +'_{:.2e}'.format(j)
            tVG_name = tVG_name +'_'+ name +'_{:.2e}'.format(j)

        str_lst.append(str_line + '-JV_file '+JV_name+ '.dat -Var_file '+Var_name+'.dat -scPars_file '+scPars_name+'.dat')
        JV_files.append(Path(path2SIMsalabim) / str(JV_name+ '.dat'))
        Var_files.append(Path(path2SIMsalabim) / str(Var_name+ '.dat'))
        scPars_files.append(Path(path2SIMsalabim) / str(scPars_name+ '.dat'))


        labels.append(lab)
        JVexp_lst.append('')
        sys_lst.append(system)
        path_lst.append(path2SIMsalabim)


        # prepare for CELIV
        str_lst_CELIV.append(str_line + '-tj_file '+tj_name+ '.dat -tVG_file '+tVG_name+'.txt')
        tVG_files_CELIV.append(Path(path2ZimT) / str(tVG_name+ '.txt'))
        tj_files_CELIV.append(Path(path2ZimT) / str(tj_name+ '.dat'))
        path_lst_CELIV.append(path2ZimT)



    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots


    # Run simulation
    if run_simu_JV:
        # Run SIMsalabim
        print('Running SIMsalabim')
        run_multiprocess_simu(run_SIMsalabim,max_jobs,str_lst,sys_lst,path_lst)

        
    print('\n')
    print('Start plotting JVs:')
    idx = 0
    Voc,Jsc,FF,PCE = [],[],[],[]
    Jmax_CELIV,t_Jmax_CELIV,tau_CELIV = [],[],[]
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name,scPars_file_name,tj_file_name,tVG_file_CELIV_name in zip(str_lst,JVexp_lst,labels,JV_files,Var_files,scPars_files,tj_files_CELIV,tVG_files_CELIV):
        
        data_scpars = pd.read_csv(scPars_file_name,delim_whitespace=True) # Load scPars_file
        Jsc.append(float(data_scpars['Jsc']))
        Voc.append(float(data_scpars['Voc']))
        FF.append(float(data_scpars['FF']))
        PCE.append(float((abs(data_scpars['Jsc'])/10)*data_scpars['Voc']*data_scpars['FF']))


        ## Plot JVs
        if plot_JVs:
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(exp_name,delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()

            SIMsalabim_JVs_plot(num_JV_plot,data_JV,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False)
        
        if run_simu_CELIV: # prepare tVG for CELIV simulation that starts at Voc
            # zimt_CELIV(tmin,tmax,float(data_scpars['Voc']),slope,Gen,tpulse,tstep,tdelay,width_pulse = 6e-9,time_exp=True,steps=100,tVG_name=tVG_file_CELIV_name)
            zimt_CELIV(tmin,tmax,Voffset,slope,Gen,tpulse,tstep,tdelay,width_pulse = 6e-9,time_exp=True,steps=100,tVG_name=tVG_file_CELIV_name)  
        
        idx = idx+1
    print('\n')
    print('Done plotting JVs.')
    Jsc,Voc,FF,PCE = np.asarray(Jsc),np.asarray(Voc),np.asarray(FF),np.asarray(PCE)
    ## Run transients
    #    
    if run_simu_CELIV:
        print('\n')
        print('Running CELIV')
        run_multiprocess_simu(run_zimt,max_jobs,str_lst_CELIV,sys_lst,path_lst_CELIV)
    
    print('\n')
    print('Start plotting CELIVs:')
    idx = 0
    for tj_file_CELIV_name,tVG_file_CELIV_name in zip(tj_files_CELIV,tVG_files_CELIV):
        ## Plot CELIVs
        if plot_CELIVs:
            
            data_tj = pd.read_csv(tj_file_CELIV_name,delim_whitespace=True)
            Jmax_CELIV.append(data_tj['Jext'][data_tj['Jext'].idxmin()]/10)
            t_pic = data_tj['t'][data_tj['Jext'].idxmin()]
            t_Jmax_CELIV.append(t_pic)
            data2fit = data_tj.copy()
            data2fit['t'] = data2fit['t']-t_pic 
            data2fit = data2fit[data2fit.t >= t_pic]

            # perform the fit
            p0 = (1e-6, min(data2fit['Jext']), data2fit['Jext'].iloc[-1]) # start with values near those we expect
            params, cv = scipy.optimize.curve_fit(MonoExpDecay,data2fit['t'], data2fit['Jext'], p0)
            k, A, B = params
            tau_CELIV.append(1/k)
            plt.figure(num_CELIV_plot)
            plt.plot(data2fit['t'], MonoExpDecay(data2fit['t']-t_pic, k,A,B)/10, '--',color=colors[idx])
            zimt_tj_plot(num_CELIV_plot,data_tj,colors=colors[idx],plot_type=0,legend=False)
        idx = idx+1
    print('\n')
    print('Done plotting JVs.')
    
   

    # JV_file plots
    plot_effs = True # Make JV plot
    if plot_effs:
        num_fig = num_fig + 1
        num_eff_plot = num_fig
        f_effs = plt.figure(num_eff_plot,figsize=size_fig)
    
        plt.figure(num_eff_plot)
        # x_var =  parameter2['values']
        # x_label  = parameter2['name']
        x_var = time
        x_label = 'Time [s]'
    
        plt.subplot(221)
        plt.semilogx(x_var,PCE/max(PCE),linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('Norm. PCE')
        plt.ylim([0.2,1.2])
        plt.grid(b=True,which='both')
        plt.subplot(222)
        plt.semilogx(x_var,FF,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('FF')
        # plt.ylim([0,1])
        plt.grid(b=True,which='both')
        plt.subplot(223)
        plt.semilogx(x_var,Voc,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('V$_{OC}$ [V]')
        # plt.ylim([0.7,1.3])
        plt.grid(b=True,which='both')
        plt.subplot(224)
        plt.semilogx(x_var,abs(np.asarray(Jsc)),linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('J$_{SC}$ mA cm$^{-2}$')
        # plt.ylim([10,30])
        # plt.legend(loc='best', title="Scan direction:",fontsize=35)
        plt.tight_layout()
        plt.grid(b=True,which='both')


    # JV_file plots
    plot_CELIV_val = False # Make JV plot
    if plot_CELIV_val :
        num_fig = num_fig + 1
        num_CELIV_val = num_fig
        f_effs = plt.figure(num_CELIV_val,figsize=size_fig)
    
        plt.figure(num_CELIV_val)
        # x_var =  parameter2['values']
        # x_label  = parameter2['name']
        x_var = time
        x_label = 'Time [s]'
    
        plt.subplot(131)
        plt.plot(x_var,Jmax_CELIV,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('J$_{pic}$ [mA cm$^{-2}$]')
        # plt.ylim([2,25])
        plt.grid(b=True,which='both')
        plt.subplot(132)
        plt.plot(x_var,t_Jmax_CELIV,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('Pic position [s]')
        # plt.ylim([0,1])
        plt.grid(b=True,which='both')
        plt.subplot(133)
        plt.plot(x_var,tau_CELIV,linestyle='-',marker='o',markersize=10,markerfacecolor='w')
        plt.xlabel(x_label)
        plt.ylabel('Lifetime [s]')
        # plt.ylim([0,1])
        plt.grid(b=True,which='both')

        plt.tight_layout()
        plt.grid(b=True,which='both')
    




if __name__ == '__main__':

    

    
    time = np.geomspace(1e-1,300,13)
    time = np.insert(time,0,0)
    a = MonoExpDecay(time, 1, 5e-8, 1e-10)
    b = StretchedExp(time, 1, 3/12, 4.5e-8, 3e-8)
    run_plot_SIMsalabim(a,b)
    plt.figure(110)
    plt.plot(time,a,label='mono')
    plt.plot(time,b,label = 'stretch')
    plt.loglog(time,a+b,label='both')
    plt.legend()
    plt.show()