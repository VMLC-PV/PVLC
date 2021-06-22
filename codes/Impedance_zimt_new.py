###################################################
########## Simulate Impedance using ZimT ##########
###################################################
## by Vincent M. Le Corre
## Package import
import os,sys,platform,tqdm,parmap,multiprocessing,warnings,cmath
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mplcolors
from matplotlib.ticker import LogFormatter
from matplotlib.lines import Line2D
from scipy.integrate import simps
from scipy.optimize import curve_fit, leastsq
from scipy import interpolate
from scipy import constants
from time import time
from itertools import repeat
from impedance.models.circuits import circuits as impcir
from impedance import preprocessing as imppre
from impedance import visualization as impvis
## Don't show warnings
warnings.filterwarnings("ignore")
## Homemade package import
import plot_settings_screen
from VLC_useful_func import sci_notation,run_zimt,zimt_tj_plot,Store_output_in_folder,clean_up_output,preprocess_Impedance_data,get_complex_impedance,fit_sin_func,sin_func
from tVG_gen import zimt_impedance

## Main Program
def Impedance(): 
    ## General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()          # Operating system
    print(system)
    max_jobs = os.cpu_count()-2                       # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
            max_jobs = 1
            slash = '\\'
            try:
                os.system('taskkill.exe /F /IM zimt.exe')
                print('taskkill excecuted')
            except:
                print('no taskkill excecuted')
    else:
        slash = '/'

    curr_dir = os.getcwd()                           # Current working directory
    path2ZimT = 'Simulation_program/DDSuite_v400/ZimT'+slash  # Path to ZimT

    ## Physics constants
    q = constants.value(u'elementary charge')
    eps_0 = constants.epsilon_0

    ## Simulation input
    run_simu = True                                           # Rerun simu? (bool)
    plot_tjs = 0                                           # make tJ and tV plots? (bool)
    plot_output = 1                                        # make Nyquist and Bode plots? (bool)
    move_ouput_2_folder = True                             # (bool)
    Store_folder = 'Impedance'+slash
    clean_output = False                                   # Make output plots? (bool)
    make_fit = False                                       # make fit dn/dt (bool)
    calc_capacitance = 1                                   # calculate capacitance? (bool)
    nFcm2 = 1                                              # Capacitance unit (bool)
    L = 100e-9                                             # Device thickness (m)
    L_LTL = 20e-9 # m
    L_RTL = 20e-9 # m
    eps = 4                                              # active layer permittivity, same as in device_parameters.txt
    C_geo = (eps*eps_0/(L-L_LTL-L_RTL))                             # geometric capacitance (Ohm)
    freqs1 = np.geomspace(1e2,1e4,num=15,endpoint=False)
    freqs2 = np.geomspace(1e4,1e9,num=20,endpoint=True)
    freqs = np.append(freqs1, freqs2)                      # frequencies to simulate (Hz)
    freqs_interp = np.geomspace(1e2,1e7,num=2000)          # same but more dense range of frequencies (for interpolation)
    # print(freqs)
    Vapps = [0]#[-.5, .5]  #[-1, 0, 0.5]  # [-2, -1, 0, 0.2, 0.4, 0.6, 0.8]         # Applied voltage (V)
    Vamp = 0.01                       # Amplitude voltage perturbation
    Gen = 0e27                        # Average generation rate 
    sun = 1                      # generation rate at 1 sun

    ## Figure control
    plottitle = 'ZimT: Organic Solar Cell  {:.0f}nm  {:.1f}sun'.format(L*1e9, Gen/sun)  # full title of plots
    savetitle = 'OSC_{:.0f}nm_{:.1f}sun'.format(L*1e9, Gen/sun)    # start of filename of plots
    size_fig = (13, 9)
    size_fig_nyquist = (20, 9)
    colors = cm.viridis((np.linspace(0, 1, max(len(freqs), 4)+1)))  # Freq colors (virid-ish)
    colors1 = cm.winter(np.linspace(0, 1, max(len(Vapps), 4)+1))    # Vapp colors (blue-ish)
    colors2 = cm.autumn(np.linspace(0, 1, max(len(Vapps), 4)+1))    # Vapp colors (red-ish) for second axis
    colors1[1] = colors1[0]  # Color adjustment for dark colors: skip second color
    colors2[1] = colors2[0]  # same
    colors1 = colors1[1:]    # same
    colors2 = colors2[1:]    # same
    color1 = colors1[0]      # Color of left axis
    color2 = colors2[0]      # Color of right axis

    ## Initialize 
    str_lst,sys_lst,path_lst,tj_lst,tVG_lst = [],[],[],[],[]
    start = time()
    print('Vapps:', str(Vapps))

    ###########################################################################
    ### Run ZimT ##############################################################
    ###########################################################################

    if run_simu:
        print('sim loop...')
        # Generate tVG files 
        for freq in freqs:
            for Vapp in Vapps:
                
                CB = 4 # eV
                VB = 6 # eV
                W_L = 5 # eV
                W_R = 5 # eV
                mu = 1e-7 # m^2/Vs, zero field mobility
                zimt_impedance(Vapp,Vamp,freq,Gen,steps=200,tVG_name=path2ZimT
                               +'tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq))
                str_lst.append('-eps_r '+str(eps)+' -L '+str(L)+' -L_LTL '+str(L_LTL)+' -L_RTL '+str(L_RTL)+' -tVG_file '
                               +'tVG_{:.2f}V_f_{:.1e}Hz.txt '.format(Vapp,freq)
                               +'-tj_file '
                               +'tj_{:.2f}V_f_{:.1e}Hz.dat '.format(Vapp,freq))
                sys_lst.append(system)
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_{:.2f}V_f_{:.1e}Hz.txt'.format(Vapp,freq))
                tj_lst.append('tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapp,freq))
        
        # Run ZimT
        # str_lst = str_lst[::-1]  # reverse list order to start with longest delays
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_zimt,list(zip(str_lst,sys_lst,path_lst)),
                                 pm_pool=p, pm_processes=max_jobs,pm_pbar=True)
        p.close()
        p.join()
    
    print('Calculation time {:.2f} s'.format(time() - start)) # Time in ,seconds

    
    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    ###########################################################################
    ### Calculate Complex Impedance ###########################################
    ###########################################################################

    Zs,ReZ,ImZ,Zmag = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    Cap,R,phase = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    fZ_file = [[] for f in Vapps]
    Jm, tm, Jmo, tmo = [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps], [[] for f in Vapps]
    print('Z loop...')
    for idx in range(len(Vapps)):
        for freq in freqs:
            with open(path2ZimT+Store_folder  # open in read mode
                      +'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapps[idx],freq), 'r') as file:
                data_tj2 = pd.read_csv(file, delim_whitespace=True)
            Jmo[idx].append(data_tj2['Jext'])  # save original data for comparison
            tmo[idx].append(data_tj2['t'])
            data_tj2 = preprocess_Impedance_data(data_tj2,freq)
            Jm[idx].append(data_tj2['Jext'])  # save preprocessed data for comparison
            tm[idx].append(data_tj2['t'])
            comp = get_complex_impedance(data_tj2,freq)
            ## save
            Zs[idx].append(comp)
            ReZ[idx].append(comp.real)
            ImZ[idx].append(comp.imag)
            Zmag[idx].append(abs(comp))
            phase[idx].append(cmath.phase(comp))
            Cap[idx].append((1/comp).imag/(2*np.pi*freq))
        ## save f and Z in file for later use of impedance.py
        fZ = np.transpose([freqs, np.asarray(ReZ[idx]), np.asarray(ImZ[idx])])
        fZ = fZ[np.asarray(ImZ[idx])<0]  # only save realistic data: Im(Z)<0
        fZ_file[idx] = path2ZimT+Store_folder+'fZ_{:.2f}V.txt'.format(Vapps[idx])
        np.savetxt(fZ_file[idx], fZ, delimiter=',')  # ',' needed for impedance.preprocessing

    ## convert lists to arrays to enable calculations
    Zs, ReZ, ImZ, Cap = np.asarray(Zs), np.asarray(ReZ), np.asarray(ImZ), np.asarray(Cap)

    ###########################################################################
    ### Fits ##################################################################
    ###########################################################################

    if calc_capacitance:
        ## interpolate ImZ and ReZ over freq to gain Cap_max
        f_max, Cap_max, interpI, interpR, Rp, Rs = [], [], [], [], [], []
        for idx in range(len(Vapps)):
            intI = interpolate.interp1d(freqs, -ImZ[idx], kind='cubic')
            intR = interpolate.interp1d(freqs, ReZ[idx], kind='cubic')
            fmax = freqs_interp[np.argmax(intI(freqs_interp))]  # interpolated freq of max(-ImZ)
            rs = np.min(ReZ[idx])   # assumes Rs=min(ReZ)
            rp = 2*(intR(fmax)-rs)  # assumes ReZ=Rs+Rp/2 @ max(-ImZ)
            capmax = 1/(2*np.pi*fmax*rp)
            ## save
            interpI.append(intI)
            interpR.append(intR)
            f_max.append(fmax)
            Rp.append(rp)
            Rs.append(rs)
            Cap_max.append(capmax)

        ## fit to equivalent circuit with impedance.py
        circuit_model = 'R0-p(R1,C1)'      # important: no blanks ' ' in string
        initial_guess = [1e3, 1e5, C_geo]  # guess for circuit components (SI)
        Z_fit, circ_fit = [[] for f in Vapps], []
        for idx in range(len(Vapps)):
            ## import from file OR directly:
            # frequencies, Z = freqs, [ReZ[idx], ImZ[idx]]
            frequencies, Z = imppre.readFile(fZ_file[idx])

            ## define and calculate circuit components
            circuit = impcir.CustomCircuit(circuit=circuit_model, initial_guess=initial_guess)
            circuit.fit(frequencies, Z)
            print(Vapps[idx], 'V')
            print(circuit)  # print initial guess and fit parameters
            print('geometric capacitance:', C_geo, 'F', '\n')  # print C_geo for comparison

            ## save
            Z_fit[idx] = circuit.predict(freqs_interp)   # predict Z for more dense freq array
            circ_fit.append(circuit.parameters_)         # fit parameters for circuit components

            ## Nyquist plot for every Vapp
            fig, ax = plt.subplots()
            impvis.plot_nyquist(ax, Z, fmt='o')                       # pot zimt data
            impvis.plot_nyquist(ax, np.asarray(Z_fit[idx]), fmt='-')  # plot circuit fit
            plt.legend(['Zimt Data', 'Circuit Fit'])
            plt.title('Nyquist at {:.1f}V'.format(Vapps[idx]))
            plt.axis('equal')     # 1:1 axis aspect
            ax.axhline(0, c='k')  # plot x-axis for better overview

        ## convert lists to arrays to enable calculations
        Z_fit, circ_fit, Cap_max = np.asarray(Z_fit), np.asarray(circ_fit), np.asarray(Cap_max)
      
    ###########################################################################
    ### J/t, V/t graphs #######################################################
    ###########################################################################

    if plot_tjs:
        print('tJ loop...')
        lines = ['-', '--', ':', '-.', 'x-', 'x--', 'x:', 'x-.', 'o-', 'o--', 'o:',
                'o-.', '+-', '+--', '+-.', '+:']  # linestyles for up to 16 Vapps
        for idx in range(len(Vapps)):
            for i in range(len(freqs)):
                with open(path2ZimT+Store_folder  # open in read mode
                        +'tj_{:.2f}V_f_{:.1e}Hz.dat'.format(Vapps[idx],freqs[i]), 'r') as file:
                    data_tj = pd.read_csv(file, delim_whitespace=True)
                ## norm data to interval (-1, 1)
                data_tj['Jext_norm'] = data_tj['Jext']/max(data_tj['Jext'])*10
                data_tj['Vext_norm'] = data_tj['Vext']/max(data_tj['Vext'])*10
                data_tj['t_norm'] = data_tj['t']/max(data_tj['t'])*10
                if i%int((len(freqs)-1)/2) == 0:    # plot 3 lines per Vapp
                    ## plot normed J/t
                    zimt_tj_plot(100, data_tj, x='t_norm', y=['Jext_norm'],
                                ylimits= [-1.1,1.1],
                                labels=sci_notation(freqs[i],sig_fig=0)+' Hz, '+str(Vapps[idx])+' V',
                                colors=colors[i], line_type=[lines[idx]],
                                plot_type=0, save_yes=True, legend=True, figsize=size_fig,
                                title=plottitle+'  Normed J')
                    ## plot normed V/t
                    zimt_tj_plot(101, data_tj, x='t_norm', y=['Vext_norm'],
                                ylimits= [-1.1,1.1],
                                labels=sci_notation(freqs[i],sig_fig=0)+' Hz '+str(Vapps[idx])+' V',
                                colors=colors[i], line_type=[lines[idx]],
                                plot_type=0, save_yes=True, legend=True, figsize=size_fig,
                                title=plottitle+' Normed V')

                ## check whether the sin fit works
                if i%int((len(freqs)-1)/2) == 0:  # test sin fit of 3 curves per Vapp
                    try:
                        Jf_amp, Jf_freq, Jf_phi, Jf_off = fit_sin_func(np.asarray(tm[idx][i]), np.asarray(Jm[idx][i]), freqs[i])
                        Jf = Jf_amp*np.sin(Jf_freq*2*np.pi*tm[idx][i] + Jf_phi) + Jf_off  # fitted J
                        plt.figure()
                        plt.plot(tmo[idx][i], Jmo[idx][i], 'o', label='original')
                        plt.plot(tm[idx][i], Jm[idx][i], 'o', label='preprocessed')
                        plt.plot(tm[idx][i], Jf, label='fitted')
                        plt.xlabel('t')
                        plt.ylabel('J')
                        plt.legend()
                        plt.title('freq {:.1e}Hz Vapp {:.1f}V'.format(freqs[i], Vapps[idx]))
                    except:
                        print('tJ fit plot failed.')

    ###########################################################################
    ### Plots #################################################################
    ###########################################################################

    if plot_output:
        print('plotting...')
        ## Nyquist plot Im(Z)/Re(Z)
        try: 
            ax = [[] for _ in range(len(Vapps))]
            for idx in range(len(Vapps)):
                fig = plt.figure(figsize=size_fig_nyquist, frameon=False)
                ax[idx] = fig.add_subplot(aspect='equal')
                ## circuit fitted and interpolated graphs and lines at Rs+Rp
                intplot = ax[idx].plot(interpR[idx](freqs_interp), interpI[idx](freqs_interp),
                                   c=colors2[idx], label=str(Vapps[idx])+' V, interpolated fit', zorder=-idx)
                intRsRp = ax[idx].axvline(Rs[idx]+Rp[idx], c=colors2[idx], zorder=-idx)
                intRs   = ax[idx].axvline(Rs[idx], c=colors2[idx], zorder=-idx)
                fitplot = ax[idx].plot(Z_fit[idx].real, -Z_fit[idx].imag, '--',
                          c=colors2[idx], zorder=-idx, label='equivalent circuit fit')
                fitRsRp = ax[idx].axvline(circ_fit[idx][0]+circ_fit[idx][1], ls='--', c=colors2[idx], zorder=-idx)
                fitRs   = ax[idx].axvline(circ_fit[idx][0], ls='--', c=colors2[idx], zorder=-idx)
                ## zimt data
                cmap = cm.viridis
                cmap.set_bad('black',0.)
                # freqs_plot = [freqs for _ in Vapps]
                scatter = ax[idx].scatter(ReZ[idx], -ImZ[idx], c=freqs, cmap=cmap, norm=mplcolors.LogNorm(), 
                                    zorder = 100, label='ZimT')
                ## color bar and layout
                formatter = LogFormatter(10, labelOnlyBase=False) 
                cbar = plt.colorbar(scatter, format=formatter, shrink=1)
                cbar.ax.set_ylabel('Frequency  f  [Hz]')
                ax[idx].legend(fontsize = 'small')
                ax[idx].set_ylim(bottom=0, top=.5*plt.xlim()[1])
                plt.xlabel('Resistance  Z\'  '+r'[$\Omega$]')
                plt.ylabel(r'-$\,$Reactance  -$\,$Z'+'\'\'  '+r'[$\Omega$]')
                plt.title(plottitle)
                # plt.tight_layout()
                plt.savefig(path2ZimT+Store_folder+savetitle+'_Impedance_Nyquist_'+str(Vapps[idx])+'V.png',
                            dpi=300, bbox_inches='tight', transparent=True)
        except:
            print('Nyquist Plot failed.')

        ## Bode plot Z/f and phi/f
        try:
            fig, ax1 = plt.subplots(figsize=size_fig)
            ax1.set_xlabel('Frequency  f  [Hz]')
            ax1.set_ylabel(r'Impedance  |Z|  [$\Omega$]', color=color1)
            ax1.set_yscale('log')
            ax2 = ax1.twinx()  # instantiate second axes that shares the same x-axis
            ax2.set_ylabel(r'Phase  $\phi$  [rad]', color=color2)
            for idx in range(len(Vapps)):
                ax1.semilogx(freqs, Zmag[idx],'o-',color=colors1[idx], label='{:.1f} V'.format(Vapps[idx]))
                ax2.semilogx(freqs, phase[idx],'o-',color=colors2[idx], label='{:.1f} V'.format(Vapps[idx]))
            ax1.legend(title='|Z|', fontsize='small', loc='upper center')  # set location to prevent overlay of legends
            ax2.legend(title=r'$\phi$', fontsize='small')
            plt.title(plottitle)
            plt.tight_layout()
            plt.savefig(path2ZimT+Store_folder+savetitle+'_Impedance_Z_phi_f.png',
                        dpi=300,transparent=True, bbox_inches='tight')
        except:
            print('Bode plot Z,phi/f failed')
        
        ## Bode plot Re(Z)/f
        try:
            fig2, ax3 = plt.subplots(figsize=size_fig)
            ax3.set_xlabel('Frequency  f  [Hz]')
            ax3.set_ylabel('Resistance  Z\'  '+r'[$\Omega$]')
            ax3.set_yscale('log')
            for i in range(len(Vapps)):
                ax3.semilogx(freqs, ReZ[i],'o-',color=colors1[i], label='{:.1f} V'.format(Vapps[i]))
            plt.legend(fontsize='small')
            plt.title(plottitle)
            plt.savefig(path2ZimT+Store_folder+savetitle+'_Impedance_ReZ_f.png',
                        dpi=300,transparent=True, bbox_inches='tight')
        except:
            print('Bode plot Re(Z)/f failed')

        ## Bode plot Im(Z)/f
        try:
            fig3, ax4 = plt.subplots(figsize=size_fig)
            ax4.set_xlabel('Frequency  f  [Hz]')
            plt.ylabel(r'-$\,$Reactance  -$\,$Z'+'\'\'  '+r'[$\Omega$]')
            ax4.set_yscale('log')
            for i in range(len(Vapps)):
                ax4.semilogx(freqs, -ImZ[i],'o-',color=colors2[i], label='{:.1f} V'.format(Vapps[i]))
            plt.legend(fontsize='small')
            plt.title(plottitle)
            plt.savefig(path2ZimT+Store_folder+savetitle+'_Impedance_ImZ_f.png',
                        dpi=300,transparent=True, bbox_inches='tight')
        except:
            print('Bode plot Im(Z)/f failed')

        ## Bode plot C/f
        if nFcm2:
            Cap_max = Cap_max/(1e-9/1e-4)  # to get nF/cm^2
            Cap = Cap/(1e-9/1e-4)  # to get nF/cm^2
            C_geo = C_geo/(1e-9/1e-4)  # to get nF/cm^2
            for idx in range(len(Vapps)):
                circ_fit[idx][2] = circ_fit[idx][2]/(1e-9/1e-4)
                    
        try:
            fig4, ax5 = plt.subplots(figsize=size_fig)
            ax5.set_xlabel('Frequency  f  [Hz]')
            ax5.set_ylabel('Capacitance  C  [F]')
            if nFcm2:
                ax5.set_ylabel(r'Capacitance  C  [nF/cm$^2$]')
                ax5.set_ylim(-5, 83)
                ax5.set_xlim(1e2, 1e9)
            else:
                ax5.set_yscale('log')
            for idx in range(len(Vapps)):
                ax5.axhline(circ_fit[idx][2], ls='--', c=colors1[idx])
                ## Cap_max will be totally off when semi-circle hasn't reached its maximum
                if np.isclose(circ_fit[idx][2], Cap_max[idx], rtol=.1):
                    ax5.axhline(Cap_max[idx], ls=':', c=colors1[idx])
                ax5.semilogx(freqs, Cap[idx], 'o-', color=colors1[idx], label='{:.1f} V, ZimT'.format(Vapps[idx]))

            ## black lines for legend
            ax5.plot([], [], ls='--', c='k', label='equivalent circuit fit')
            ax5.plot([], [], ls=':', c='k', label='Nyquist interpolation')
            ax5.axhline(C_geo, c='k', label='geometric capacitance')
            plt.legend(fontsize='small')
            plt.title(plottitle)
            plt.savefig(path2ZimT+Store_folder+savetitle+'_Impedance_C_f.png',
                        dpi=300,transparent=True, bbox_inches='tight')
        except:
            print('Bode plot C/f failed')

    plt.show()


    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',path2ZimT+Store_folder)
        clean_up_output('tVG',path2ZimT+Store_folder)
        print('Ouput data was deleted from '+path2ZimT+Store_folder)

    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds


if __name__ == '__main__':
    Impedance()
