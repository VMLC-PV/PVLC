###################################################
########## Simulate TDCF using ZimT ###############
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from time import time
from scipy import constants,interpolate
from scipy.integrate import simps
from scipy.optimize import curve_fit
from matplotlib.lines import Line2D
# Import homemade package by VLC
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *
from VLC_units.cleanup_folder.clean_folder import *
from VLC_units.useful_functions.aux_func import *
from VLC_useful_func import TDCF_fit

# Main Program
def TDCF(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = False, plot_tjs = True, move_ouput_2_folder = True, Store_folder = 'TDCF',clean_output = False,verbose = True):
    """Run TDCF simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_TDCF function (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'TDCF'
    clean_output : bool, optional
        Clean up output?, by default False
    verbose : bool, optional
        Verbose?, by default True
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
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass

    
    curr_dir = os.getcwd()              # Current working directory
    if path2ZimT is None:
        path2ZimT = os.path.join(os.getcwd(),'Simulation_program/SIMsalabim_v425/ZimT')                  # Path to ZimT in curr_dir

    make_fit = False                                            # Make fit dn/dt

    ## Physics constants
    q = constants.value(u'elementary charge')
    eps_0 = constants.value(u'electric constant')

    ## TDCF Inputs
    # see zimt_tdcf in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No TDCF input dictionary given, using default values')
          
        tmin = 1e-10                                                # First time step after 0
        tcol = 30e-6                                                 # tmax = tcol + tdelay  = final time step
        Vpres = [0.99]#[-0.5,0,0.2,0.4,0.6,0.7,0.8,0.85,0.9,1]#,0.85]      # Initial applied voltage (steady-state)
        Vcol = -2.5
        Gens = [1e21]                                               # Max generation rate for the gaussian laser pulse
        tpulse = 3e-9                                               # Time at which the pulse occurs, needs to be changed if width_pulse is modified
        tstep = 1e-9
        delays = np.geomspace(5e-9,1e-6,num=10) #[6e-9]                                             # Delays for TDCF simulation
        width_pulse = 6e-9
        time_exp = True
        steps = 200
    else:
        if 'tmin' in input_dic.keys():
            tmin = input_dic['tmin']
        if 'tcol' in input_dic.keys():
            tcol = input_dic['tcol']
        if 'Vpres' in input_dic.keys():
            Vpres = input_dic['Vpres']
        if 'Vcol' in input_dic.keys():
            Vcol = input_dic['Vcol']
        if 'Gens' in input_dic.keys():
            Gens = input_dic['Gens']
        if 'tpulse' in input_dic.keys():
            tpulse = input_dic['tpulse']
        if 'tstep' in input_dic.keys():
            tstep = input_dic['tstep']
        if 'delays' in input_dic.keys():
            delays = input_dic['delays']
        if 'width_pulse' in input_dic.keys():
            width_pulse = input_dic['width_pulse']
        if 'time_exp' in input_dic.keys():
            time_exp = input_dic['time_exp']
        if 'steps' in input_dic.keys():
            steps = input_dic['steps']

    ## Prepare strings to run
    # Fixed string
    if fixed_str is None:
        if verbose:
            print('No fixed string given, using default value')
        fixed_str = ''  # add any fixed string to the simulation command

    # Get thicknesses (needed later)
    ParFileDic = ReadParameterFile(f"{path2ZimT}/device_parameters.txt") # Read device parameters
    L = float(ChosePar('L',GetParFromStr(fixed_str),ParFileDic))
    L_LTL = float(ChosePar('L_LTL',GetParFromStr(fixed_str),ParFileDic))
    L_RTL = float(ChosePar('L_RTL',GetParFromStr(fixed_str),ParFileDic))

    ## Figures control
    num_fig = 0
    size_fig = (16, 10)
    colors = cm.viridis((np.linspace(0,1,max(len(Vpres),4)+1)) ) # Color range for plotting
    save_fig = True

    # plot input tVG_file (light) 
    plot_input_tVG = True                                      
    if plot_input_tVG and run_simu:
        num_fig_input_tVG = num_fig
        num_fig = num_fig + 1
        f_input_tVG = plt.figure(num_fig_input_tVG,figsize=size_fig)
    
    # plot photocurrent (Jtdcf)
    plot_PhotoCurr_TDCF = True                                      
    if plot_PhotoCurr_TDCF:
        num_fig_PhotoCurr = num_fig
        num_fig = num_fig + 1
        f_PhotoCurr = plt.figure(num_fig_PhotoCurr,figsize=size_fig)

    # plot extracted charges (calculated by integrating photocurrent)
    plot_extract_charges = True                                      
    if plot_extract_charges:
        num_fig_extract_charges = num_fig
        num_fig = num_fig + 1
        f_extract_charges = plt.figure(num_fig_extract_charges,figsize=size_fig)
    
    # plot extracted charges vs Vpre
    plot_extract_charges_vs_V = True
    if plot_extract_charges and plot_extract_charges_vs_V:
        num_fig_extract_charges_vs_V  = num_fig
        num_fig = num_fig + 1
        f_extract_charges_vs_V  = plt.figure(num_fig_extract_charges_vs_V ,figsize=size_fig)

    # plot normalized densities (ntot, npre and ncol) 
    plot_Norm_Dens = True                                      
    if plot_Norm_Dens:
        num_fig_Norm_Dens = num_fig
        num_fig = num_fig + 1
        f_Norm_Dens = plt.figure(num_fig_Norm_Dens,figsize=size_fig)
    
    # plot differential density (dntot/dt) 
    plot_Diff_Dens = True                                      
    if plot_Diff_Dens:
        num_fig_Diff_Dens = num_fig
        num_fig = num_fig + 1
        f_Diff_Dens = plt.figure(num_fig_Diff_Dens,figsize=size_fig)
                                 
   
    # Initialize 
    npre,ncol,ntot,str_lst,code_name_lst,path_lst,tj_lst,tVG_lst,tVG_Gen_lst,Qext_Qgen = [],[],[],[],[],[],[],[],[],[]
    idx = 0
    start = time()
    
    if run_simu:
        # Run simulation in the dark
        # Generate tVG files for Dark simulation
        for Vpre in Vpres:
            for tdelay in delays:
                zimt_tdcf(tmin,tdelay+tcol,Vpre,Vcol,0,tpulse,tstep,tdelay,width_pulse = width_pulse,tVp = 10e-9,time_exp=time_exp,steps=steps,tVG_name=os.path.join(path2ZimT,'tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt'.format(tdelay,Vpre))) 
                str_lst.append(fixed_str+' -tVG_file tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt -tj_file tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre,tdelay,Vpre))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt'.format(tdelay,Vpre))
                tj_lst.append('tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre))
            
            # Run simulation with light pulse
            for Gen in Gens:
                for tdelay in delays:
                    str_lst.append(fixed_str+' -tVG_file tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt -tj_file tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre,Gen,tdelay,Vpre))
                    zimt_tdcf(tmin,tdelay+tcol,Vpre,Vcol,Gen,tpulse,tstep,tdelay,width_pulse = width_pulse,tVp = 10e-9,time_exp=time_exp,steps=steps,tVG_name=os.path.join(path2ZimT,'tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre)))
                    code_name_lst.append('zimt')
                    path_lst.append(path2ZimT)
                    tVG_lst.append('tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre))
                    tVG_Gen_lst.append('tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre))
                    tj_lst.append('tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre))
        
        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose=verbose)

            
        ## Move output folder to new folder
        if move_ouput_2_folder: # Move outputs to Store_folder
            Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
            Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

        
    ########################################################
    ################## Plotting ############################
    ########################################################
    if run_simu:
        if plot_input_tVG:
            plt.figure(num_fig_input_tVG)
            ax1 = plt.axes()
            ax2 = ax1.twinx()
            for i in tVG_Gen_lst:
                VG_file = pd.read_csv(os.path.join(path2ZimT,Store_folder,i),delim_whitespace=True)
                ax1.semilogx(VG_file['t'],VG_file['Gehp'])
                ax2.semilogx(VG_file['t'],VG_file['Vext'],'--')
    for Vpre in Vpres:
        for tdelay in delays:
            data_tjdark = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre)),delim_whitespace=True)
            tswith = tdelay + tpulse
            npre_dumb,ncol_dumb,ntot_dumb = [],[],[]
            for Gen in Gens:
                data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre)),delim_whitespace=True)

                # Interpolate Jdark in case some time step did not converge
                tck = interpolate.splrep(data_tjdark['t'], data_tjdark['Jext'], s=0)
                data_tj['Jdark'] =interpolate.splev(data_tj['t'], tck, der=0)

                # Calculate the photocurrent Jtdcf by substracting light vs dark simulation    
                # data_tj['Jtdcf'] = abs(data_tj['Jext']-data_tjdark['Jext'])
                data_tj['Jtdcf'] = abs(data_tj['Jext']-data_tj['Jdark'])

                # Calculate Densities
                pre = data_tj[data_tj.t < tswith]
                col = data_tj[data_tj.t >= tswith]
                npre_dumb.append(abs(simps(pre['Jtdcf']/(q*(L-L_LTL-L_RTL)),x=pre['t'])))
                ncol_dumb.append(abs(simps(col['Jtdcf']/(q*(L-L_LTL-L_RTL)),x=col['t'])))
                ntot_dumb.append(abs(simps(data_tj['Jtdcf']/(q*(L-L_LTL-L_RTL)),x=data_tj['t'])))
                
                # Make plot
                if plot_PhotoCurr_TDCF:
                    label = 'Gen {:.1e} delay {:.1e}'.format(Gen,tdelay)
                    zimt_tj_plot(num_fig_PhotoCurr,data_tj,y=['Jtdcf'],xlimits= [0,max(data_tj['t'])],colors=colors[idx],plot_type=1,save_yes=save_fig,legend=False,pic_save_name = os.path.join(path2ZimT,Store_folder,'transient.jpg'))

                if plot_extract_charges:
                    plt.figure(num_fig_extract_charges)
                    extrat_charge = integrate.cumtrapz(data_tj['Jtdcf'], data_tj['t'], initial=0) #calc number of extracted charges
                    data_tj['Qext']= extrat_charge/(q*(L-L_LTL-L_RTL))
                    max_charge = max(data_tj['Qext'])
                    # data_tj['t'] = data_tj['t']-min(data_tj['t'])
                    data_tj['t1'] = data_tj['t'] - tdelay - tpulse
                    tend_charge = max(data_tj['t1'])
                    plt.semilogx(data_tj['t1'],data_tj['Qext'],color=colors[idx],label='V$_p$$_r$$_e$ = {:.2f}'.format(Vpre))

                if plot_extract_charges and plot_extract_charges_vs_V:
                    gen_charge = integrate.cumtrapz(data_tj['Jphoto'], data_tj['t'], initial=0)/(q*(L-L_LTL-L_RTL)) #calc number of photogenerated charges
                    plt.figure(num_fig_extract_charges_vs_V)
                    Qext_Qgen.append(max(data_tj['Qext'])/max(gen_charge))
                    plt.plot(Vpre,max(data_tj['Qext'])/max(gen_charge),color=colors[idx],linestyle='none',marker='o',markeredgecolor=colors[idx],markersize=10,markerfacecolor=colors[idx],markeredgewidth = 3)
                                        

            pre = data_tj[data_tj.t < tswith]
            col = data_tj[data_tj.t >= tswith]
            # Store Densities
            if npre == []:
                npre = np.asarray(npre_dumb)
                ncol = np.asarray(ncol_dumb)
                ntot = np.asarray(ntot_dumb)
            else:
                npre = np.vstack([npre,np.asarray(npre_dumb)])
                ncol = np.vstack([ncol,np.asarray(ncol_dumb)])
                ntot = np.vstack([ntot,np.asarray(ntot_dumb)])
        
        idx = idx+1


    if plot_Norm_Dens or plot_Diff_Dens:
        if ntot.ndim > 1 and len(delays) > 1: # Check if there is more than one delay
            
            # Calculate dn/dt and norm densities
            dndt, ntot_norm, npre_norm, ncol_norm = [],[],[],[]
            for i in range(len(Gens)):
                dn = np.diff(ntot[:,i])
                dt = np.diff(delays)
                dndt.append(dn/dt)
                npre_norm.append(npre[:,i]/max(ntot[:,i])) 
                ncol_norm.append(ncol[:,i]/max(ntot[:,i])) 
                ntot_norm.append(ntot[:,i]/max(ntot[:,i])) 

            dndt = np.transpose(np.asarray(dndt))
            npre_norm = np.transpose(np.asarray(npre_norm))
            ncol_norm = np.transpose(np.asarray(ncol_norm))
            ntot_norm = np.transpose(np.asarray(ntot_norm))

            # Make plots
            if plot_Norm_Dens:
                plt.figure(num_fig_Norm_Dens)
                patchs = [Line2D([0], [0], linestyle='-', color='k' ,label='n$_{tot}$'),Line2D([0], [0], linestyle='-.', color='k' ,label='n$_{pre}$'),Line2D([0], [0], linestyle='--', color='k' ,label='n$_{col}$')]
                colors2 = cm.viridis((np.linspace(0,1,max(len(Gens),4)+1)) ) 
                idx = 0
                for i in range(len(Gens)):
                    plt.semilogx(delays,ntot_norm[:,i],color=colors2[idx],linestyle='-')
                    plt.semilogx(delays,npre_norm[:,i],color=colors2[idx],linestyle='-.')
                    plt.semilogx(delays,ncol_norm[:,i],color=colors2[idx],linestyle='--')
                    patchs.append(Line2D([0], [0], linestyle='-', color=colors2[idx] ,label='G = '+sci_notation(Gens[i],sig_fig=1)+' m$^{-3}$ s$^{-1}$'))
                    idx = idx + 1
                

            if plot_Diff_Dens:    
                plt.figure(num_fig_Diff_Dens)
                plt.loglog(ncol[1::],abs(dndt),'o',label='')
                plt.loglog(np.geomspace(1e18,1e23,num=10), abs(TDCF_fit(np.geomspace(1e18,1e23,num=10), 1e-17,2e22)), 'g--',label='fit: k$_2$='+sci_notation(1e-17,sig_fig=1)+' m$^3$ s$^{-1}$, n$_{BG}$='+sci_notation(6e21,sig_fig=1)+'m$^{-3}$')
                # print(abs(TDCF_fit(np.geomspace(ncol.min(),ncol.max(),num=10), 1e-17,2e22)))
                # Make fit
                if make_fit ==True:
                    popt, pcov = curve_fit(TDCF_fit, np.asarray(ncol[1:]), np.asarray(abs(dndt)),p0=[1e-18,1e22])
                    plt.loglog(ncol[1:], abs(TDCF_fit(np.asarray(ncol[1:]), *popt)), 'b--')#,label='fit: a=%5.3f, b=%5.3f' % tuple(popt))

                

        else:
            print('\n /!\ Only one delay was simulated so dn/dt cannot be calculated \n') 
            plot_Diff_Dens = False
            plot_Norm_Dens = False


    ## Figure axis and Falselabel control
    # plot input tVG
    if plot_input_tVG and run_simu:
        plt.figure(num_fig_input_tVG)
        ax1.set_xlabel("Time [s]")
        ax1.set_ylabel('Generation rate [m$^{-3}$ s$^{-1}]$')
        ax2.set_ylabel('Voltage [V]')

    # plot extracted charges (calculated by integrating photocurrent)                    
    if plot_extract_charges:
        plt.figure(num_fig_extract_charges)
        # plt.vlines(L**2/(2*mumax*(2.5-0.86)),0,2*max_charge,label='t$_{max}$ = L$^2$/(2 $\mu_{max}$ (V$_{ext}$ - V$_{OC}$))',colors='r')
        # plt.vlines(L**2/(2*mumin*(2.5-0.86)),0,2*max_charge,label='t$_{min}$ = L$^2$/(2 $\mu_{min}$ (V$_{ext}$ - V$_{OC}$))',colors='b')
        # plt.vlines(50*5e-7*3.5*eps_0/L,0,2*max_charge,label='t$_{RC}$ = R*C$_{geo}$',colors='k') 
        plt.legend(loc='best',frameon=False,fontsize = 30,ncol=2)
        plt.grid(b=True,which='both')
        plt.xlim(5e-9,tend_charge)
        plt.ylim(0,2*max_charge)
        plt.xlabel('t [s]')
        plt.ylabel(r'n$_{ext}$ [m$^{-3}$]')
        plt.tight_layout()
        if save_fig:
            plt.savefig(os.path.join(path2ZimT,Store_folder,'chargevstime.jpg'),dpi=100,transparent=True)  

    # plot extracted over generated charges ratio                    
    if plot_extract_charges and plot_extract_charges_vs_V:
        plt.figure(num_fig_extract_charges_vs_V)
        plt.grid(b=True,which='both')
        plt.xlim(min(Vpres)-0.1,max(Vpres)+0.1)
        plt.ylim(0,1.1)
        plt.xlabel('V$_{pre}$ [V]')
        plt.ylabel('Q$_{ext}$/Q$_{gen}$')
        plt.tight_layout()
        if save_fig:
            plt.savefig(os.path.join(path2ZimT,Store_folder,'chargevsvoltage.jpg'),dpi=100,transparent=True) 


    
    # plot normalized densities (ntot, npre and ncol)                                     
    if plot_Norm_Dens:
        plt.figure(num_fig_Norm_Dens)
        plt.legend(handles=patchs,loc='best',frameon=False,fontsize = 30)
        plt.grid(b=True,which='both')
        plt.xlabel('Time [s]')
        plt.ylabel('Normalized Carrier Density')
        plt.tight_layout()
        if save_fig:
            plt.savefig(os.path.join(path2ZimT,Store_folder,'norm_charges.jpg'),dpi=310,transparent=True)
    
    # plot differential density (dntot/dt)                                     
    if plot_Diff_Dens:
        plt.figure(num_fig_Diff_Dens)
        plt.legend(loc='best',frameon=False,fontsize = 30)
        plt.grid(b=True,which='both')
        plt.xlabel('n$_{col}$ [m$^{-3}$]')
        plt.ylabel(r'$\| \frac{dn}{dt}\|$ [m$^{-3}$ s$^{-1}$]')
        plt.tight_layout()
        if save_fig:
            plt.savefig(os.path.join(path2ZimT,Store_folder,'dndt.jpg'),dpi=100,transparent=True)


    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds    
    return Vpres,Qext_Qgen,num_fig
    


if __name__ == '__main__':
    
    TDCF()

    plt.show()
    