######################################################################
############### Function to run simulation code ######################
######################################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.legend import Legend
from scipy import stats,optimize,constants
import subprocess,shutil,os,glob,tqdm,parmap,multiprocessing,random,sys,platform
from itertools import repeat
import warnings
import sys
# package by VLC
from VLC_units.simu.compile_prog import *
# Don't show warnings
warnings.filterwarnings("ignore")


def run_code(name_prog,path2prog,str2run='',show_term_output=False,verbose=False):
    """Run program 'name_prog' in the folder specified by 'path2prog'.

    Parameters
    ----------
    name_prog : str
        name of the program to run.
    
    path2prog : str
        path to the folder containing zimt 
        ('./zimt' in Linux and 'zimt.exe' in Windows ).

    st2run : str
        String to run for the name_prog.
    
    show_term_output : bool, optional
        If True, show the terminal output of the program, by default False.
    
    verbose : bool, optional
        Verbose?, by default False.
    
    Returns
    -------
    
    """
    
    System = platform.system()                  # Operating system
    is_windows = (System == 'Windows')          # Check if we are on Windows
    path2prog = str(path2prog)                  # Convert to string
    curr_dir = os.getcwd()                      # Get current directory
    os.chdir(path2prog)                         # Change directory to working directory

    if show_term_output == True:
        output_direct = None
    else:
        output_direct = subprocess.DEVNULL
    
    if is_windows:
        cmd_list = name_prog.lower()+'.exe ' + str2run
        if not os.path.isfile(path2prog+'\\'+name_prog.lower()+'.exe'):
            fpc_prog(name_prog,path2prog,show_term_output=False,force_fpc=False,verbose=verbose)
    else : 
        cmd_list = './'+name_prog.lower()+' ' + str2run
        if not os.path.isfile('./'+path2prog+'\\'+name_prog.lower()):
            fpc_prog(name_prog,path2prog,show_term_output=False,force_fpc=False,verbose=verbose)
    try:
        
        subprocess.check_call(cmd_list.split(), encoding='utf8', stdout=output_direct, cwd=path2prog, shell=is_windows)
    except subprocess.CalledProcessError:
        print(path2prog)
        raise ChildProcessError
    os.chdir(curr_dir)                          # Change directory back to original directory

def run_multiprocess_simu(prog2run,code_name_lst,path_lst,str_lst,max_jobs=max(1,os.cpu_count()-1)):
    """run_multiprocess_simu runs simulations in parrallel (if possible) on max_jobs number of threads

    Parameters
    ----------
    prog2run : function
        name of the function that runs the simulations

    code_name_lst : list of str
        list of names of the codes to run
    
    str_lst : list of str
        List containing the strings to run for the simulations

    path_lst : list of str
        List containing the path to the folder containing the simulation program
    
    max_jobs : int, optional
        Number of threads used to run the simulations, by default os.cpu_count()-1
    """
    p = multiprocessing.Pool(max_jobs)
    results = parmap.starmap(prog2run,list(zip(code_name_lst,path_lst,str_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)
    p.close()
    p.join()

if __name__ == '__main__':
    
    System = platform.system()                  # Operating system
    is_windows = (System == 'Windows')          # Check if we are on Windows
    
    if is_windows:
        run_code('SimSS','c:/Users/lecor/Desktop/GitHub/PVLC/codes/Simulation_program/SIMsalabim_v425/SimSS',str2run='-L 110e-9')
    else:
        run_code('SimSS','/mnt/c/Users/lecor/Desktop/GitHub/PVLC/codes/Simulation_program/SIMsalabim_v425/SimSS',str2run='-L 110e-9')