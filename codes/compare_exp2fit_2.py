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
# from VLC_useful_func import *
# Import homemade package by VLC
from VLC_units.plots.SimSS_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.SCLC.SCLC_func import *
import VLC_units.plots.plot_settings_screen
from VLC_units.useful_functions.aux_func import *
# import plot_settings_screen

def GetParFromStr(str2run,L = 100e-9,Cn=1e-13,Cp=1e-13):

    str2run = ' '.join(str2run.split()) #remove extra white space

    new = str2run.split()

    Names= new[::2]
    Values = new[1::2]
    for i,j in enumerate(Names):
        Names[i] = Names[i].replace('-', '')

    for i,j in enumerate(Values):
        Values[i] = float(Values[i])

    return pd.DataFrame({'Parameters':Names,'Values':Values})



def run_plot_SIMsalabim2():

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

    # path2SIMsalabim = os.path.join(os.getcwd() ,'Simulation_program/DDSuite_v409_crack/SIMsalabim')
    # path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/DDSuite_v409_crack/SIMsalabim')
    path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_Xioyan/SimSS')
    run_simu = 1 #Rerun simulation 

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
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.008774872893385296, 0.0295868067205203, 0.09674039280460424, 0.2903605440421194, 0.9438975920716443,1, 2.2588639120448035]
    
    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 3.063E-5 -Rshunt 7.995E-1 -Gehp 1.217E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.353 -Nc 5.645E+27 -W_L 4.060 -St_R 6.496E+9 -St_L 9.067E+8 -CB 4.04' #Best 1%

      
    # age = 'aged'
    # pixel = 3
    # Gfracs = [0,0.008125240951214394, 0.02704509630728539, 0.08851876613458132, 0.26551237818238604, 0.861085219869508,1, 2.051557459849598]
    
    # fixed_str = '-mun_0 1.107E-7 -mup_0 1.269E-7 -kdirect 3.483E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.212E+28 -Etrap 4.664 -Bulk_tr 7.477E+18 -W_R 5.204 -Nc 5.519E+27 -W_L 4.319 -CB 4.04'

    # fixed_str = '-Nc_LTL 5e24 -Nc_RTL 5e24 -CB_LTL 4.34 -W_L 4.34 -St_L 5e14 -Gehp 1.412E+28 -tolJ 1e-2'
    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.412E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.204 -Nc 5.645E+27 -W_L 4.219 -St_R 6.496E+15 -St_L 9.067E+15 -CB 4.04' 


    # pixel = 1
    # Gfracs = [0,0.00866783401714054, 0.029238534889572412, 0.09543346872544842, 0.2862747313125504, 0.9306521928901167,1, 2.226981200864729]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    ## BTP-4F-12
    ID = 3937
    # age = 'fresh'
    # pixel = 2
    # Gfracs = [0,0.0090015933301563, 0.029796718028699196, 0.09717869847550556, 0.29087180672086876, 0.943525832647383, 1, 2.2428264040319434]


    # fixed_str = '-mun_0 2.312E-8 -mup_0 2.247E-7 -kdirect 8.000E-18 -Rseries 5.468E-5 -Rshunt 2.521 -Gehp 1.341E+28 -Etrap 4.587 -Bulk_tr 4.894E+18 -W_R 5.367 -Nc 5.224E+27 -W_L 4.110 -St_R 3.152E+6 -St_L 4.313E+8 -CB 4.06' #keep'Parameters to vary
    
    # pixel = 3
    # Gfracs = [0,0.009016735212175199, 0.029660936641206906, 0.0968623846255899, 0.29014650417906485, 0.9410713947273658, 1, 2.241490248848625]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    age = 'aged'
    pixel = 2
    Gfracs = [0,0.007814523992331497, 0.026764391395073882, 0.08785804792342802, 0.26391800166747525, 0.8572188961690461, 1,2.042784401098461]
 
    # fixed_str = '-mun_0 5.471E-8 -mup_0 1.734E-8 -kdirect 7.703E-18 -Rseries 1.092E-5 -Rshunt 1.241 -Gehp 1.384E+28 -Etrap 4.626 -Bulk_tr 6.462E+18 -W_R 5.442 -Nc 6.311E+27 -W_L 4.156 -St_R 2.733E+5 -St_L 4.750E+10'

    fixed_str = '-kdirect 2.826e-17 -Bulk_tr 2.011e+19 -mun_0 5.682e-08 -mup_0 8.056e-08 -Rseries 5.332e-05 -Rshunt 1.241 -W_L 4.190e+00 -W_R 5.476e+00 -Nc 6.055e+27 -CB 4.06 -Gehp 1.394E+28'


    ## Y6-BO-4F
    # ID = 3916
    # age = 'fresh'
    # pixel = 4
    # Gfracs = [0,0.00899285762350159, 0.02958847222355735, 0.09647014717333, 0.2897613119670855, 0.9399050246570602, 1, 2.2396445153662032]

    # fixed_str = '-mun_0 2.442E-8 -mup_0 1.829E-7 -kdirect 7.999E-18 -Rseries 2.652E-5 -Rshunt 1.370 -Gehp 1.326E+28 -Etrap 4.558 -Bulk_tr 6.827E+18 -W_R 5.485 -Nc 5.799E+27 -W_L 4.074 -St_R 1.336E+8 -St_L 9.698E+9 -CB 4.05'

    # pixel = 2
    # Gfracs = [0,0.009009535782328836, 0.029625731507062908, 0.09649315148591558, 0.2897571956906255, 0.9400858561495319,1, 2.238969531299488]
    # fixed_str = '-mun_0 7.002E-8 -mup_0 2.733E-8 -kdirect 8.000E-18 -Rseries 1.182E-5 -Rshunt 1.560 -Gehp 1.320E+28 -Etrap 4.683 -Bulk_tr 6.336E+18 -W_R 5.424 -Nc 5.739E+27 -W_L 4.085 -St_R 2.074E+5 -St_L 6.013E+4 -CB 4.05'


    # o-IDFBR
    # ID = 3939
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.007703621989336275, 0.030391616105901822, 0.10494576208861922, 0.31542562971134397, 0.994778451921309, 1, 2.2770546056260343]
    # fixed_str = '-mun_0 1.678E-7 -mup_0 8.802E-7 -kdirect 2.816E-18 -Rseries 1.591E-6 -Rshunt 1.210 -Gehp 3.657E+27 -Etrap 4.319 -Bulk_tr 3.476E+20 -W_R 5.152 -Nc 1.120E+26 -W_L 3.826 -St_R 8.787E+5 -St_L 8.865E+4 -CB 3.7'
    
    parameter1 = {'name':'L','values':[140e-9]}
    parameter2 = {'name':'Gfrac','values':Gfracs}
    parameter3 = {'name':'L_LTL','values':[30e-9]}
    parameter4 = {'name':'L_RTL','values':[10e-9]}
    L_LTL = parameter3['values'][0] # needed for nrj_diag plot
    L_RTL = parameter4['values'][0] # needed for nrj_diag plot
    parameters = [parameter1,parameter2,parameter3,parameter4] 

    
    str_lst,labels,JVexp_lst,JV_files,Var_files,sys_lst,path_lst,val,nam,code_name_lst = [],[],[],[],[],[],[],[],[],[]
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
        code_name_lst.append('SimSS')
        # code_name_lst.append('SIMsalabim')
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
        run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        # run_multiprocess_simu(run_SIMsalabim,max_jobs,str_lst,sys_lst,path_lst)


    idx = 0
    perf_exp,perf_simu = [],[]
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name,Gfrac in zip(str_lst,JVexp_lst,labels,JV_files,Var_files,parameter2['values']):

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(os.path.join(path2SIMsalabim , exp_name),delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()

            if Gfrac > 0:
                perf_simu.append([get_Voc(data_JV['Vext'],data_JV['Jext']),get_Jsc(data_JV['Vext'],data_JV['Jext']),get_FF(data_JV['Vext'],data_JV['Jext']),get_PCE(data_JV['Vext'],data_JV['Jext'],suns=Gfrac)/10])
                perf_exp.append([get_Voc(data_JVexp['V'],data_JVexp['J']),get_Jsc(data_JVexp['V'],data_JVexp['J']),get_FF(data_JVexp['V'],data_JVexp['J']),get_PCE(data_JVexp['V'],data_JVexp['J'],suns=Gfrac)/10])
            SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp,xlimits=[-0.5,1.1],ylimits=[-50,10],save_yes=True,pic_save_name=os.path.join(path2SIMsalabim ,'JVfit.jpg'))
            # SIMsalabim_JVs_plot(num_JV_plot,data_JV,plot_type=0,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=True,data_JVexp=data_JVexp,xlimits=[-0.5,1.35],ylimits=[-15,10],save_yes=True,pic_save_name=path2SIMsalabim/'JVfit.jpg')
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
    plt.savefig(os.path.join(path2SIMsalabim ,'perf_fit.jpg'))

    # FOMS = CalcVLCFOM(fixed_str,os.path.join(path2SIMsalabim ,'Device_parameters.txt'))
    # print('Figures of Merits:')
    # print('delta_b = {:.2e}'.format(FOMS[0]))
    # print('delta_tr = {:.2e}'.format(FOMS[2]))
    # print('Delta_b = {:.2e}'.format(FOMS[1]))
    # print('Delta_tr = {:.2e}'.format(FOMS[3]))

    # FOMsave = pd.DataFrame({'delta_b':[FOMS[0]],'delta_tr':[FOMS[2]],'Delta_b':[FOMS[1]],'Delta_tr':[FOMS[3]],'Vbi':[FOMS[4]],'Rseries':[FOMS[5]],'Rshunt':[FOMS[6]]})
    # FOMsave.to_csv(Path(path2SIMsalabim / 'FOMS.txt'),index=False)


if __name__ == '__main__':

    run_plot_SIMsalabim2()
    plt.show()