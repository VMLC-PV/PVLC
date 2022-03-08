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
from VLC_units.SCLC.SCLC_func import *
import VLC_units.plots.plot_settings_screen
from VLC_units.useful_functions.aux_func import *

from scipy.interpolate import make_interp_spline


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
        path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_SCLC/SimSS')                    # Path to SimSS
                                                                                            # Rerun simulation

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ['']  # add any fixed string to the simulation command

    # Parameters to vary
    if parameters is None:
        parameters = []
        parameters.append({'name':'Gfrac','values':[0]})

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
    idx2 = 0
    for fix in fixed_str: # loop over fixed strings

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
            str_lst.append(fix+ ' ' +str_line+ '-JV_file '+JV_name+'_'+str(idx2)+'.dat -Var_file '+Var_name+'.dat -scPars_file '+scPars_name+'_'+str(idx2)+'.dat')# -ExpJV '+JVexp_lst[idx])
            JV_files.append(os.path.join(path2SIMsalabim , str(JV_name+'_'+str(idx2)+'.dat')))
            Var_files.append(os.path.join(path2SIMsalabim , str(Var_name+ '.dat')))
            scPars_files.append(os.path.join(path2SIMsalabim , str(scPars_name+'_'+str(idx2)+ '.dat')))
            labels.append(lab)
            if not add_exp:
                JVexp_lst.append('')
            code_name_lst.append('SimSS')
            path_lst.append(path2SIMsalabim)
            idx += 1
        idx2 += 1

    # print(str_lst)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(str_lst),4)+1)) # prepare color for plots

    # Run simulation
    if run_simu:
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True)


    # Read parameters file (useful latter)
    ParFileDic = ReadParameterFile(os.path.join(path2SIMsalabim , 'device_parameters.txt'))
    idx = 0
    V1f,J1f,V2f,J2f,Vinf,Jinf,n1,n2,ninf,mu = [],[],[],[],[],[],[],[],[],[]
    for Simu_str,exp_name,lbl,JV_file_name,Var_file_name,scPars_file_name in zip(str_lst,JVexp_lst,labels,JV_files,Var_files,scPars_files):

        ## Plot JVs
        if plot_JVs:
            # print(JV_file_name)
            data_JV = make_df_JV(JV_file_name) # Load JV_file
            if plot_exp: # Load Exp JV
                data_JVexp = pd.read_csv(os.path.join(path2SIMsalabim, exp_name),delim_whitespace=True)
            else:
                data_JVexp = pd.DataFrame()


            SCLC_res1=Make_SCLC_plot(num_JV_plot,data_JV,plot_type=3,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=plot_exp,data_JVexp=data_JVexp,pic_save_name=os.path.join(path2SIMsalabim , str('JV'+ext_save_pic)),save_yes=True,show_tangent=[4],ylimits=[1e-6,1e7],xlimits=[1e-2,30])
            

            plt.figure(num_JV_plot+1,figsize=size_fig)
            Spline = make_interp_spline(SCLC_res1[0],SCLC_res1[2])
            V_ = np.geomspace(SCLC_res1[0].min(), SCLC_res1[0].max(), 500)
            slope_ = Spline(V_)
            # plt.semilogx(SCLC_res1[0],SCLC_res1[2],color=colors[idx])
            plt.semilogx(V_,slope_,color=colors[idx])
            # plt.axvline(SCLC_res1[13],linestyle=':',color=colors[idx])

            plt.figure(num_JV_plot)
            
            ## Get input values
            ParStrDic = GetParFromStr(Simu_str)
            Bulk_tr = float(ChosePar('Bulk_tr',ParStrDic ,ParFileDic))
            CNI = float(ChosePar('CNI',ParStrDic ,ParFileDic))
            CPI = float(ChosePar('CPI',ParStrDic ,ParFileDic))
            eps_r = float(ChosePar('eps_r',ParStrDic ,ParFileDic))
            mun_0 = float(ChosePar('mun_0',ParStrDic ,ParFileDic))
            L = float(ChosePar('L',ParStrDic ,ParFileDic))
            L_LTL = float(ChosePar('L_LTL',ParStrDic ,ParFileDic))
            L_RTL = float(ChosePar('L_RTL',ParStrDic ,ParFileDic))
            T = float(ChosePar('T',ParStrDic ,ParFileDic))
            Nc = float(ChosePar('Nc',ParStrDic ,ParFileDic))
            CB = float(ChosePar('CB',ParStrDic ,ParFileDic))
            W_L = float(ChosePar('W_L',ParStrDic ,ParFileDic))
            W_R = float(ChosePar('W_R',ParStrDic ,ParFileDic))
            Vbi = abs(W_L-W_R)
            Vmax = float(ChosePar('Vmax',ParStrDic ,ParFileDic))
            phi = abs(CB-W_L)
            Rseries = float(ChosePar('Rseries',ParStrDic ,ParFileDic))
            ## Calc Voltages
            Vnet = calc_Vnet_with_ions(CNI,Bulk_tr,L-L_LTL-L_RTL,eps_r)
            Vtfl = calc_Vtfl(Bulk_tr,L-L_LTL-L_RTL,eps_r)
            Vsat = calc_Vsat(L-L_LTL-L_RTL,Nc,phi,eps_r,T)
            Ntmin = calc_nt_min(L,eps_r,T)
            Vmin = calc_Vtfl(Ntmin,L-L_LTL-L_RTL,eps_r)
            # plt.axvline(Vmin,linestyle='--',color='k')
            # plt.axvline(Vnet,linestyle='--',color=colors[idx])
            # plt.axvline(Vtfl,linestyle='-.',color=colors[idx])
            # if Vmax > Vsat:
            #     plt.axvline(Vsat,linestyle=':',color=colors[idx])

            # Calculate densities
            n1.append(calc_net_charge(SCLC_res1[9],L,eps_r))
            n2.append(calc_net_charge(SCLC_res1[11],L,eps_r))
            ninf.append(calc_net_charge(SCLC_res1[13],L,eps_r))


            # fit Mott-Gurney to get mobility
            # data_JV_mott = data_JV[data_JV['Vext']>=Vtfl]
            data_JV_mott = data_JV
            if Vsat > Vtfl:
                data_JV_mott = data_JV_mott[data_JV_mott['Vext']<=Vsat]
                # print(len(data_JV_mott))
                data_JV_mott['Vext'] = data_JV_mott['Vext']# - Rseries*data_JV_mott['Jext'] - Vbi

                # SCLC_res2=Make_SCLC_plot(num_JV_plot,data_JV_mott,plot_type=3,x='Vext',y=['Jext'],colors=colors[idx],labels=labels[idx],legend=False,plot_jvexp=plot_exp,data_JVexp=data_JVexp,pic_save_name=os.path.join(path2SIMsalabim , str('JV'+ext_save_pic)),save_yes=True,show_tangent=[4],ylimits=[1e-6,1e7],xlimits=[1e-2,30],line_type = ['None'],mark='x')
                Mott_Gurney_fit = fit_MottGurney(data_JV_mott['Vext'],data_JV_mott['Jext'],mun_0,eps_r,Vbi,L,var2fit=['mu'])
                mu.append(Mott_Gurney_fit[0])
                print('Mott_Gurney_fit:',Mott_Gurney_fit)
            else:
                # mu.append(np.nan)
                mu.append(0)

            print('Vnet = {:.2e}'.format(Vnet))
            print('Vtfl = {:.2e}'.format(Vtfl))
            print('Vsat = {:.2e}'.format(Vsat))

            

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
        plt.figure(num_JV_plot+1)
        plt.xlabel('Applied Voltage [V]')
        plt.ylabel('Slope')
        plt.tight_layout()
        plt.grid(b=True,which='both')
        plt.savefig( os.path.join(path2SIMsalabim,'slope.jpg'),dpi=100,transparent=True)
    return JV_files,Var_files,scPars_files,n1,n2,ninf,Ntmin,mu,mun_0
  


if __name__ == '__main__':
    
    path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_SCLC/SimSS')

    # fixed_strs,net = [],[]
    # xaxis = 'Trap Density [cm$^{-3}$]'
    # loop = np.geomspace(5e21,1e22,6)
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr '+str(i)+' -CNI 0 -CPI 0')
    #     net.append(i)
    # loop = loop/1e6 # convert to cm^-3
    # net = np.asarray(net)/1e6 # convert to cm^-3
    
    # fixed_strs,net = [],[]
    # xaxis = 'Net Charge Density [cm$^{-3}$]'
    # loop = np.geomspace(5e20,5e21,10)
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr '+str(1e22)+' -CNI '+str(i)+' -CPI '+str(i))
    #     net.append(1e22-i)
    # loop = loop/1e6 # convert to cm^-3
    # net = np.asarray(net)/1e6 # convert to cm^-3

    # fixed_strs,net = [],[]
    # xaxis = 'Series Resistance [$\Omega$ m$^{2}$]'
    # loop = np.geomspace(1e-10,1e-4,6)
    # for i in np.geomspace(1e-10,1e-4,6):
    #     fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -Rseries '+str(i))
    #     net.append(1e22-3e21)
    # net = np.asarray(net)/1e6 # convert to cm^-3

    # fixed_strs,net = [],[]
    # loop = np.linspace(3.9,4.2,7, endpoint=True)
    # xaxis = 'Injection Barrier [eV]'
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -W_L '+str(i)+' -W_R '+str(i))
    #     net.append(1e22-3e21)
    # loop = loop-3.9

    # fixed_strs,net = [],[]
    # loop = np.linspace(3.9,4.2,7, endpoint=True)
    # xaxis = 'V$_{BI}$ [V]'
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -W_L 3.9 -W_R '+str(i))
    #     net.append(1e22-3e21)
    # loop = loop - 3.9

    # fixed_strs,net = [],[]
    # loop = np.geomspace(1e-10,1e-4,7, endpoint=True)
    # xaxis = '$\mu_{TL}$ [cm$^{2}$V$^{-1}$s$^{-1}$]'
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -tolJ 1e-2 -L 820e-9 -L_LTL 10e-9 -L_RTL 10e-9 -mob_LTL '+str(i)+' -mob_RTL '+str(i))
    #     net.append(1e22-3e21)

    
    # fixed_strs,net = [],[]
    # loop = np.asarray([3,10,25])
    # xaxis = '$\epsilon_{TL}$'
    # for i in loop:
    #     fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -tolJ 1e-2 -L 820e-9 -L_LTL 10e-9 -L_RTL 10e-9 -mob_LTL 1e-6 -mob_RTL 1e-6 -eps_r_LTL '+str(i)+' -eps_r_RTL '+str(i))
    #     net.append(1e22-3e21)

    fixed_strs,net = [],[]
    loop = np.linspace(3.9,4.0,3)
    xaxis = 'Offset [eV]'
    for i in loop:
        fixed_strs.append('-Bulk_tr 1e22 -CNI 3e21 -CPI 3e21 -L 820e-9 -L_LTL 10e-9 -L_RTL 10e-9 -mob_LTL 1e-6 -mob_RTL 1e-6 -tolJ 1e-2 -NP 5000 -W_L '+str(i)+' -W_R '+str(i)+' -CB_LTL '+str(i)+' -CB_RTL '+str(i))
        net.append(1e22-3e21)
    print(fixed_strs)
    JV_files,Var_files,scPars_files,n1,n2,ninf,Ntmin,mu,mun_0 = run_plot_SIMsalabim(fixed_str = fixed_strs,plot_nrj_diag=False,plot_densities=False,run_simu=1)
    colors = plt.cm.viridis(np.linspace(0,1,max(len(net),3)+1))
    
    # print(n1)
    # print(n2)
    # print(ninf)
    ind = -1
    plt.figure(10,figsize=(12,7.5))
 
    plt.scatter(loop,np.asarray(n1)/1e6, c=colors[:ind], s=100, marker='s', label='V$_1$')
    plt.scatter(loop,np.asarray(n2)/1e6, c=colors[:ind], s=100, marker='^', label='V$_2$')
    plt.scatter(loop,np.asarray(ninf)/1e6, c=colors[:ind], s=100, marker='o', label='V$_{inf}$')
    # plt.scatter(net,np.asarray(n1)/1e6, c=colors[:ind], s=100, marker='s', label='V$_1$')
    # plt.scatter(net,np.asarray(n2)/1e6, c=colors[:ind], s=100, marker='^', label='V$_2$')
    # plt.scatter(net,np.asarray(ninf)/1e6, c=colors[:ind], s=100, marker='o', label='V$_{inf}$')
    plt.axhline(y=1e22/1e6,color='k',linestyle='-',label='N$_t$')
    plt.axhline(y=(1e22-3e21)/1e6,color='k',linestyle='--',label='N$_{net}$')

    plt.grid(b=True,which='both',zorder=-1)
    # plt.axvspan(6e13, Ntmin/1e6, color='lightgray', alpha=0.9,zorder=-2)
    # plt.axvspan(Ntmin/1e6, 3*Ntmin/1e6, color='lightgray', alpha=0.4,zorder=-3)
    # plt.loglog([6e13 ,3e16],[6e13 ,3e16],'k--',label='N$_{net}$')
    # plt.xlim([6e14 ,3e16])
    # plt.ylim([6e14 ,3e16])
    # plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(xaxis)
    plt.ylabel('Calc. Density [cm$^{-3}$]')
    plt.legend(fontsize=28,bbox_to_anchor=(0,1.02,1,0.2), loc="lower left", mode="expand", borderaxespad=0, ncol=5,frameon=False)  #loc='best',ncol=3
    plt.tight_layout()
    plt.savefig( os.path.join(path2SIMsalabim,'net_charge.jpg'),dpi=100,transparent=True)

    plt.figure(11,figsize=(12,7.5))
    ticks = []
    for i in loop:
        ticks.append('{:.2f}'.format(i))
        # ticks.append(sci_notation(i*1e4, sig_fig=-1))
    mu = np.asarray(mu)*1e4
    plt.bar(np.arange(len(mu)),mu,color=colors[:ind])
    # plt.bar(loop,mu*1e4,color=colors[:ind],width=loop[1]-loop[0])
    plt.axhline(mun_0*1e4,linestyle='-',color='k')
    plt.xticks(np.arange(len(mu)), ticks,fontsize=32)
    # plt.tick_params(
    # axis='x',          # changes apply to the x-axis
    # which='both',      # both major and minor ticks are affected
    # bottom=False,      # ticks along the bottom edge are off
    # top=False,         # ticks along the top edge are off
    # labelbottom=False) # labels along the bottom edge are off
    # plt.xscale('log')
    plt.yscale('log')
    plt.grid(b=True,which='both',axis='y',zorder=-1)
    plt.ylim([1e-4 ,max(mu)*10])
    plt.ylabel('Mobility [cm$^{2}$V$^{-1}$s$^{-1}$]')
    plt.xlabel(xaxis)
    plt.tight_layout()
    plt.savefig( os.path.join(path2SIMsalabim,'mobility.jpg'),dpi=100,transparent=True)
    
    


        


    plt.show()