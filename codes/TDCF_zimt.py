###################################################
########## Simulate TDCF using ZimT ###############
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,platform,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.lines import Line2D
from scipy.integrate import simps
from scipy import integrate, interpolate
from scipy.optimize import curve_fit
from scipy import constants
from time import time
from itertools import repeat,product
# Don't show warnings
warnings.filterwarnings("ignore")
# Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,TDCF_fit,Store_output_in_folder,clean_up_output
from tVG_gen import zimt_tdcf

# Main Program
def TDCF(L,L_LTL=0,L_RTL=0,num_fig=0,str2run='',path2ZimT='',Store_folder=''):
    """Run time-delayed collection field (TDCF) simulations

    Parameters
    ----------
    L : float
        Device total thickness (unit: m)

    L_LTL : float
        Left transport layer thickness (unit: m), by default 0
    
    L_RTL : float
        Right transport layer thickness (unit: m), by default 0
    
    num_fig : int, optional
        Starting figure number for the first plot, for the other figure the figure number will be incremented as num_fig + 1, by default 0

    str2run : str, optional
        Input string to run with ZimT,  by default ''

    path2ZimT : str, optional
        Path to folder containing ./ZimT or zimt.exe in current directory, by default = ''

    Store_folder : str, optional
        Path to folder where the output of the simulation is stored in path2ZimT directory, by default = ''

    Returns
    -------
    

    """
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

    curr_dir = os.getcwd()                                      # Current working directory
    path2ZimT = path2ZimT+slash                                 # Path to ZimT in curr_dir
    Store_folder = Store_folder+slash                           # Path to folder containing output data in curr_dir
    run_simu = False                                            # Rerun simu?  
    move_ouput_2_folder = True                                  # Move (True) output of simulation to Store_folder
    clean_output = False                                        # Clean output after simulation (delete files)
    make_fit = False                                            # Make fit dn/dt

    ## Physics constants
    q = constants.value(u'elementary charge')
    eps_0 = constants.value(u'electric constant')

    ## Simulation input for zimt_tdcf (see tVG_gen.py)   
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

    ## Figures control
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
    npre,ncol,ntot,str_lst,sys_lst,path_lst,tj_lst,tVG_lst,tVG_Gen_lst,Qext_Qgen = [],[],[],[],[],[],[],[],[],[]
    idx = 0
    start = time()
    
    if run_simu:
        # Run simulation in the dark
        # Generate tVG files for Dark simulation
        for Vpre in Vpres:
            for tdelay in delays:
                zimt_tdcf(tmin,tdelay+tcol,Vpre,Vcol,0,tpulse,tstep,tdelay,width_pulse = width_pulse,tVp = 10e-9,time_exp=time_exp,steps=steps,tVG_name=curr_dir+slash+path2ZimT+'tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt'.format(tdelay,Vpre)) 
                str_lst.append('-FailureMode 2 -tolJ 1e-3 -L '+str(L)+str2run+' -tVG_file tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt -tj_file tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre,tdelay,Vpre))
                sys_lst.append(system)
                path_lst.append(curr_dir+slash+path2ZimT)
                tVG_lst.append('tVG_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.txt'.format(tdelay,Vpre))
                tj_lst.append('tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre))
            
            # Run simulation with light pulse
            for Gen in Gens:
                for tdelay in delays:
                    str_lst.append('-FailureMode 2 -tolJ 1e-3 -L '+str(L)+str2run+' -tVG_file tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt -tj_file tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre,Gen,tdelay,Vpre))
                    zimt_tdcf(tmin,tdelay+tcol,Vpre,Vcol,Gen,tpulse,tstep,tdelay,width_pulse = width_pulse,tVp = 10e-9,time_exp=time_exp,steps=steps,tVG_name=curr_dir+slash+path2ZimT+'tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre))
                    sys_lst.append(system)
                    path_lst.append(curr_dir+slash+path2ZimT)
                    tVG_lst.append('tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre))
                    tVG_Gen_lst.append('tVG_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.txt'.format(Gen,tdelay,Vpre))
                    tj_lst.append('tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre))
        # print(str_lst)            
        # Run ZimT
        # str_lst = str_lst[::-1] # reverse list order to start with longest delays
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_zimt,list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)
        p.close()
        p.join()
            
        ## Move output folder to new folder
        if move_ouput_2_folder: # Move outputs to Store_folder
            Store_output_in_folder(tVG_lst,Store_folder,curr_dir+slash+path2ZimT)
            Store_output_in_folder(tj_lst,Store_folder,curr_dir+slash+path2ZimT)

        
    ########################################################
    ################## Plotting ############################
    ########################################################
    if run_simu:
        if plot_input_tVG:
            plt.figure(num_fig_input_tVG)
            ax1 = plt.axes()
            ax2 = ax1.twinx()
            for i in tVG_Gen_lst:
                VG_file = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+i,delim_whitespace=True)
                ax1.semilogx(VG_file['t'],VG_file['Gehp'])
                ax2.semilogx(VG_file['t'],VG_file['Vext'],'--')
    for Vpre in Vpres:
        for tdelay in delays:
            data_tjdark = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_TDCF_dark_delay_{:.1e}_Vpre_{:.2f}.dat'.format(tdelay,Vpre),delim_whitespace=True)
            tswith = tdelay + tpulse
            npre_dumb,ncol_dumb,ntot_dumb = [],[],[]
            for Gen in Gens:
                data_tj = pd.read_csv(curr_dir+slash+path2ZimT+Store_folder+'tj_TDCF_Gen_{:.1e}_delay{:.1e}_Vpre_{:.2f}.dat'.format(Gen,tdelay,Vpre),delim_whitespace=True)

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
                    zimt_tj_plot(num_fig_PhotoCurr,data_tj,y=['Jtdcf'],xlimits= [0,max(data_tj['t'])],colors=colors[idx],plot_type=1,save_yes=save_fig,legend=False,pic_save_name = curr_dir+slash+path2ZimT+Store_folder+'transient.jpg')

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
                plt.loglog(np.geomspace(1e18,1e23,num=10), abs(TDCF_fit(np.geomspace(1e18,1e23,num=10), 1.5e-16,6e21)), 'g--',label='fit: k$_2$='+sci_notation(1.5e-16,sig_fig=1)+' m$^3$ s$^{-1}$, n$_{BG}$='+sci_notation(6e21,sig_fig=1)+'m$^{-3}$')
                # print(abs(TDCF_fit(np.geomspace(ncol.min(),ncol.max(),num=10), 1e-17,2e22)))
                # Make fit
                if make_fit ==True:
                    popt, pcov = curve_fit(TDCF_fit, np.asarray(ncol[1:]), np.asarray(abs(dndt)),p0=[3e-18,2e22])
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
            plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'chargevstime.jpg',dpi=300,transparent=True)  

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
            plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'chargevsvoltage.jpg',dpi=300,transparent=True) 


    
    # plot normalized densities (ntot, npre and ncol)                                     
    if plot_Norm_Dens:
        plt.figure(num_fig_Norm_Dens)
        plt.legend(handles=patchs,loc='best',frameon=False,fontsize = 30)
        plt.grid(b=True,which='both')
        plt.xlabel('Time [s]')
        plt.ylabel('Normalized Carrier Density')
        plt.tight_layout()
        if save_fig:
            plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'norm_charges.jpg',dpi=300,transparent=True)
    
    # plot differential density (dntot/dt)                                     
    if plot_Diff_Dens:
        plt.figure(num_fig_Diff_Dens)
        plt.legend(loc='best',frameon=False,fontsize = 30)
        plt.grid(b=True,which='both')
        plt.xlabel('n$_{col}$ [m$^{-3}$]')
        plt.ylabel(r'$\| \frac{dn}{dt}\|$ [m$^{-3}$ s$^{-1}$]')
        plt.tight_layout()
        if save_fig:
            plt.savefig(curr_dir+slash+path2ZimT+Store_folder+'dndt.jpg',dpi=300,transparent=True)


    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',curr_dir+slash+path2ZimT+Store_folder)
        clean_up_output('tVG',curr_dir+slash+path2ZimT+Store_folder)
        print('Ouput data was deleted from '+curr_dir+slash+path2ZimT+Store_folder)
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds    
    return Vpres,Qext_Qgen,num_fig
    


if __name__ == '__main__':
    
    TDCF(L=140e-9,L_LTL=20e-9,L_RTL=20e-9,num_fig=0,path2ZimT = 'Simulation_program/DDSuite_v403_OPV/ZimT',Store_folder='TDCF')

    plt.show()
    