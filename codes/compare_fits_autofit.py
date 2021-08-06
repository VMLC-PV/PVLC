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
    
    ## DT-Y6
    ID = 3918
    age = 'fresh'
    pixel = 3
    Gfracs = [0,0.008774872893385296, 0.0295868067205203, 0.09674039280460424, 0.2903605440421194, 0.9438975920716443,1, 2.2588639120448035]
    # fixed_str = '-Nc 1.212E+27 -mun_0 2.888E-7 -mup_0 3.054E-7 -W_L 4.040 -W_R 5.405 -Bulk_tr 6.001E+19 -Etrap 4.453 -kdirect 8.491E-18 -Rseries 1.759E-5 -Rshunt 7.622E-1 -Gehp 1.269E+28 -UseExpData 0 -Vmin -0.5 -CB 4.04 -VB 5.51 -CB_LTL 4.04 -VB_RTL 5.51' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    # fixed_str = '-Nc 8.000E+27 -mun_0 1.000E-7 -mup_0 1.502E-7 -W_L 4.040 -W_R 5.510 -Bulk_tr 2.310E+18 -Etrap 4.551 -kdirect 3.713E-18 -Rseries 1.000E-5 -Rshunt 7.687E-1 -Gehp 1.193E+28 -UseExpData 0 -Vmin -0.5 -CB 4.04 -VB 5.51 -CB_LTL 4.04 -VB_RTL 5.51'

    fixed_str = '-Nc 1.000E+26 -mun_0 1.000E-7 -mup_0 1.502E-7 -W_L 4.220 -W_R 5.510 -Bulk_tr 1.310E+19 -Etrap 4.551 -kdirect 3.713E-18 -Rseries 1.000E-5 -Rshunt 7.687E-1 -Gehp 1.193E+28 -UseExpData 0 -Vmin -0.5 -CB 4.22 -VB 5.51 -CB_LTL 4.22 -VB_RTL 5.51'

    # age = 'aged'
    # pixel = 3
    # Gfracs = [0,0.008125240951214394, 0.02704509630728539, 0.08851876613458132, 0.26551237818238604, 0.861085219869508,1, 2.051557459849598]
    # fixed_str = '-Nc 8.212E+27 -mun_0 2.888E-7 -mup_0 3.054E-7 -W_L 4.040 -W_R 5.405 -Bulk_tr 9.001E+19 -Etrap 4.453 -kdirect 8.491E-18 -Rseries 1.759E-5 -Rshunt 4.622E-1 -Gehp 1.329E+28 -UseExpData 1 -Vmin -0.5 -CB 4.04 -VB 5.51 -CB_LTL 4.04 -VB_RTL 5.51' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    # fixed_str = '-Nc 7.100E+27 -mun_0 2.084E-7 -mup_0 3.996E-7 -W_L 4.048 -W_R 5.489 -Bulk_tr 1.010E+20 -Etrap 4.579 -kdirect 1.655E-18 -Rseries 1.118E-5 -Rshunt 6.600E-1 -Gehp 1.360E+28 -UseExpData 0 -Vmin -0.5 -CB 4.04 -VB 5.51 -CB_LTL 4.04 -VB_RTL 5.51'

    # pixel = 1
    # Gfracs = [0,0.00866783401714054, 0.029238534889572412, 0.09543346872544842, 0.2862747313125504, 0.9306521928901167,1, 2.226981200864729]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    ## BTP-4F-12
    # ID = 3937
    # age = 'fresh'
    # pixel = 2
    # Gfracs = [0,0.0090015933301563, 0.029796718028699196, 0.09717869847550556, 0.29087180672086876, 0.943525832647383, 1, 2.2428264040319434]

    # fixed_str = '-Nc 9.656E+26 -mun_0 2.151E-7 -mup_0 4.719E-7 -W_L 4.061 -W_R 5.491 -Bulk_tr 1.224E+20 -Etrap 4.367 -kdirect 3.269E-18 -Rseries 4.736E-5 -Rshunt 1.105 -Gehp 1.303E+28 -UseExpData 0 -Vmin -0.5 -CB 4.06 -VB 5.51 -CB_LTL 4.06 -VB_RTL 5.51' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)
    # Parameters to vary
    
    # pixel = 3
    # Gfracs = [0,0.009016735212175199, 0.029660936641206906, 0.0968623846255899, 0.29014650417906485, 0.9410713947273658, 1, 2.241490248848625]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    ## Y6-BO-4F
    # ID = 3916
    # age = 'fresh'
    # pixel = 4
    # Gfracs = [0,0.00899285762350159, 0.02958847222355735, 0.09647014717333, 0.2897613119670855, 0.9399050246570602, 1, 2.2396445153662032]
    # fixed_str = '-Nc 1.036E+27 -mun_0 4.051E-7 -mup_0 1.907E-7 -W_L 4.060 -W_R 5.50 -Bulk_tr 9.024E+20 -Etrap 4.281 -kdirect 1.744E-18 -Rseries 2.889E-5 -Rshunt 1.000 -Gehp 1.301E+28 -UseExpData 0 -Vmin -0.5 -CB 4.05 -VB 5.51 -CB_LTL 4.05 -VB_RTL 5.51'

    # pixel = 2
    # Gfracs = [0,0.009009535782328836, 0.029625731507062908, 0.09649315148591558, 0.2897571956906255, 0.9400858561495319,1, 2.238969531299488]
    # fixed_str = '-Nc 1.036E+27 -mun_0 4.051E-7 -mup_0 1.907E-7 -W_L 4.060 -W_R 5.50 -Bulk_tr 9.024E+20 -Etrap 4.281 -kdirect 1.744E-18 -Rseries 2.889E-5 -Rshunt 1.000 -Gehp 1.301E+28 -UseExpData 0 -Vmin -0.5 -CB 4.05 -VB 5.51 -CB_LTL 4.05 -VB_RTL 5.51'


    ## o-IDFBR
    # ID = 3939
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.007703621989336275, 0.030391616105901822, 0.10494576208861922, 0.31542562971134397, 0.994778451921309, 1, 2.2770546056260343]
    # fixed_str = '-Nc 2.560E+26 -mun_0 2.708E-7 -mup_0 5.000E-7 -W_L 3.700 -W_R 5.510 -Bulk_tr 1.523E+20 -Etrap 4.413 -kdirect 2.281E-18 -Rseries 4.263E-3 -Rshunt 0.1964 -Gehp 3.480E+27 -UseExpData 0 -Vmin -0.5 -CB 3.7 -VB 5.51 -CB_LTL 3.7 -VB_RTL 5.51'
    # fixed_str = '-St_L 1e15 -St_R 1e15 -Nc 5.000E+26 -mun_0 1.095E-7 -mup_0 1.145E-8 -W_L 3.8 -W_R 5.45 -kdirect 1.220E-14 -Rseries 2.949E-3 -Rshunt 9.170E-1 -Gehp 3.555E+27 -Bulk_tr 0e18 -accDens 0.1 -UseExpData 0 -Vmin -0.5 -CB 3.8 -VB 5.51 -CB_LTL 3.8 -VB_RTL 5.51 -Cn 1e-9 -Cp 1e-9'
    # fixed_str = '-Nc 1.000E+26 -mun_0 1.095E-7 -mup_0 1.145E-8 -W_L 4.132 -W_R 5.292 -kdirect 1.220E-18 -Rseries 2.949E-3 -Rshunt 9.170E-1 -Gehp 3.555E+27 -UseExpData 0 -Vmin -0.5 -CB 3.7 -VB 5.51 -accDens 0.5 -L 100e-9 -L_LTL 0 -L_RTL 0'
    # fixed_str = '-rms_mode log -rms_threshold 0.7  -Nc 1E+26 -mun_0 8.002E-8 -mup_0 2.626E-8 -W_L 4.073 -W_R 5.237 -kdirect 1.682E-18 -Rseries 2.416E-3 -Rshunt 1.804E-1 -Gehp 2.699E+27 -accDens 0.5 -L 100e-9 -L_LTL 0 -L_RTL 0 -Bulk_tr 0E+20'
    
    parameter1 = {'name':'L','values':[140e-9]}
    parameter2 = {'name':'Gfrac','values':Gfracs}
    parameter3 = {'name':'L_LTL','values':[30e-9]}
    parameter4 = {'name':'L_RTL','values':[10e-9]}
    L_LTL = parameter3['values'][0] # needed for nrj_diag plot
    L_RTL = parameter4['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter2,parameter3,parameter4] 

    
    str_lst,labels,JVexp_lst,JV_files,Var_files,sys_lst,path_lst,val,nam = [],[],[],[],[],[],[],[],[]
    # JVexp_lst= ['PM6_3918_dark.txt','PM6_3918_int3.txt','PM6_3918_int10.txt','PM6_3918_int33.txt','PM6_3918_int100.txt','PM6_3918_int330.txt','PM6_3918_am15_long.txt','PM6_3918_int800.txt']
    JVexp_lst= ['PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_dark.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int3.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int10.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int33.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int100.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int330.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_am15.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int800.txt']

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
    print(str_lst)
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
            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp,xlimits=[-0.5,1.4],ylimits=[-45,10],save_yes=True,pic_save_name=path2SIMsalabim/'JVfit.jpg')
            # data_JV['Jext'] =abs(data_JV['Jext'])
            # data_JVexp['J'] =abs(data_JVexp['J'])
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=3,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp)

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