###################################################
########## Simulate BACE using ZimT ###############
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
# Import homemade package by VLC
from VLC_units.plots.ZimT_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from VLC_units.make_tVG.tVG_gen import *
from VLC_units.cleanup_folder.clean_folder import *
from VLC_units.useful_functions.aux_func import *

# Main Program
def BACE(fixed_str = None, input_dic = None, path2ZimT = None, run_simu = False, plot_tjs = True, move_ouput_2_folder = True, Store_folder = 'BACE',clean_output = False,verbose = True):
    """Run BACE simulation using ZimT

    Parameters
    ----------
    fixed_str : str, optional
        Add any fixed string to the simulation command for zimt, by default None.
    input_dic : dict, optional
        Dictionary with the input for the zimt_BACE function (see tVG_gen.py), by default None.
    path2ZimT : str, optional
        Path to ZimT, by default None
    run_simu : bool, optional
        Rerun simu?, by default True
    plot_tjs : bool, optional
        make plot ?, by default True
    move_ouput_2_folder : bool, optional
        Move output to folder?, by default True
    Store_folder : str, optional
        Folder name for storing output, by default 'BACE'
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
    
    # Physics constants
    q = constants.value(u'elementary charge')

    ## BACE Inputs
    # see zimt_BACE in tVG_gen.py
    if input_dic is None:
        if verbose:
            print('No BACE input dictionary given, using default values')
        tmin = 1e-8                     # First time step after 0 [s]
        tmax = 3e-6                     # Final time step [s]
        Gens = [1.44e28]                # Average generation rate at t = 0 [m^3/s]
        Vpres = [0.77]                  # Pre-bias potential [V]
        Vextrs =[0,-1,-2,-3,-4]         # Extraction potential [V]
        steps = 100                     # Number of steps
    else:
        if 'tmin' in input_dic.keys():
            tmin = input_dic['tmin']
        if 'tmax' in input_dic.keys():
            tmax = input_dic['tmax']
        if 'Gens' in input_dic.keys():
            Gens = input_dic['Gens']
        if 'Vpres' in input_dic.keys():
            Vpres = input_dic['Vpres']
        if 'Vextrs' in input_dic.keys():
            Vextrs = input_dic['Vextrs']
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

    # Initialize 
    code_name_lst,str_lst,nextr,path_lst,tj_lst,tVG_lst = [],[],[],[],[],[]
    idx = 0
    start = time()
    p = multiprocessing.Pool(max_jobs)

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    colors = cm.viridis((np.linspace(0,1,max(len(Vextrs),4)+1)) ) # Color range for plotting
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

   
    if run_simu:    
        # Generate tVG files and str_lst for Dark simulation
        for Vpre in Vpres:
            for Vextr in Vextrs:
                zimt_BACE(tmin,tmax,0,Vpre,Vextr,time_exp=True,steps=steps,tVG_name=os.path.join(path2ZimT,'tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Vpre,Vextr)))
                str_lst.append(fixed_str+' -tVG_file tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt -tj_file tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr,Vpre,Vextr))
                code_name_lst.append('zimt')
                path_lst.append(path2ZimT)
                tVG_lst.append('tVG_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Vpre,Vextr))
                tj_lst.append('tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr))
            

        # Generate tVG files and str_lst for light pulse simulation
        for Gen in Gens:
            for Vpre in Vpres:
                for Vextr in Vextrs:
                    zimt_BACE(tmin,tmax,Gen,Vpre,Vextr,time_exp=True,steps=steps,tVG_name=os.path.join(path2ZimT,'tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Gen,Vpre,Vextr)))
                    str_lst.append('-FailureMode 1 -L '+str(L)+' -tVG_file tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt -tj_file tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr,Gen,Vpre,Vextr))
                    code_name_lst.append('zimt')
                    path_lst.append(path2ZimT)
                    tVG_lst.append('tVG_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.txt'.format(Gen,Vpre,Vextr))
                    tj_lst.append('tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr))
            

        # Run ZimT
        if do_multiprocessing:
            run_multiprocess_simu(run_code,code_name_lst,path_lst,str_lst,max_jobs)
        else:
            for i in range(len(str_lst)):
                run_code(code_name_lst[i],path_lst[i],str_lst[i],show_term_output=True,verbose=verbose)
            
        print('Calculation time {:.2f} s'.format(time() - start)) # Time in seconds

    ## Move output folder to new folder
    if move_ouput_2_folder: # Move outputs to Store_folder
        Store_output_in_folder(tVG_lst,Store_folder,path2ZimT)
        Store_output_in_folder(tj_lst,Store_folder,path2ZimT)

    ########################################################
    ################## Plotting ############################
    ########################################################
    if plot_tjs:
        for Gen in Gens:
            for Vpre in Vpres:
                nextr_dumb = []
                idx = 0
                for Vextr in Vextrs: 
                    data_tjdark = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_BACE_dark_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Vpre,Vextr)),delim_whitespace=True)
                    data_tj = pd.read_csv(os.path.join(path2ZimT,Store_folder,'tj_BACE_G_{:.1e}_Vpre_{:.2f}_Vext{:.2f}.dat'.format(Gen,Vpre,Vextr)),delim_whitespace=True)
                    data_tj['JBACE'] = abs(data_tj['Jext']-data_tjdark['Jext'])
                    zimt_tj_plot(num_fig_tjs,data_tj,y=['JBACE'],labels='G '+sci_notation(Gen, sig_fig=1)+' Vpre {:.2f} Vext {:.2f}'.format(Vpre,Vextr),colors=colors[idx],plot_type=0,save_yes=True,legend=False,pic_save_name =os.path.join(path2ZimT,Store_folder,'transient.jpg')) 
                    nextr.append(abs(simps(data_tj['JBACE']/(q*L),x=data_tj['t'])))
                    idx = idx + 1
    # print(nextr)         
            # Make plots
            f_Q = plt.figure(10,figsize=size_fig)
        
            idx = 0
            for i in range(len(Gens)):
                for j in range(len(Vpres)):
                    plt.plot(Vextrs,nextr,'o',markersize=10)
                    idx = idx + 1

            plt.grid(b=True,which='both')
            plt.xlabel('Extraction voltage [V]')
            plt.ylabel('Carrier Density [m$^{-3}$]')
            plt.tight_layout()
            plt.savefig(os.path.join(path2ZimT,Store_folder,'extrac_charges.jpg'),dpi=100,transparent=True)
    
    plt.show()
    ## Clean-up outputs from folder
    if clean_output: # delete all outputs
        clean_up_output('tj',os.path.join(path2ZimT,Store_folder))
        clean_up_output('tVG',os.path.join(path2ZimT,Store_folder))
        print('Ouput data was deleted from '+os.path.join(path2ZimT,Store_folder))
    
    print('Elapsed time {:.2f} s'.format(time() - start)) # Time in seconds
if __name__ == '__main__':
    BACE()
    
    
    

