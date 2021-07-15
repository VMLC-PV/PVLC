from pathlib import Path
import warnings, platform, os, multiprocessing, parmap
from numpy.core.function_base import linspace
from pathlib import Path
from tVG_gen import zimt_TAS
from VLC_useful_func import run_zimt, store_output_in_folder, calcZC
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy import fftpack


def TAS():
    # General Inputs
    warnings.filterwarnings("ignore")   # Don't show warnings
    system = platform.system()                  # Operating system
    max_jobs = os.cpu_count()-2                        # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
            max_jobs = 1
            try:
                os.system('taskkill.exe /F /IM zimt.exe')
            except:
                pass

    path2ZimT = Path('/mnt/c/Users/schmidt/Documents/Scripts/PVLC/codes/Simulation_program/ZimT/')     # Path to ZimT in curr_dir
    move_output = True
    store_folder = Path('TAS')

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

    Vdc = 1
    Vamp = 10e-3
    num_probe =  3
    num_periods = 4
    Gen = 0
    tVG_path = path2ZimT
    tVG_name='tVG.txt'

    fprobe_span = np.logspace(0, 6, num_probe) 

    tZ, tC = [], []

    for fprobe in fprobe_span: 

        print("Running simulation for frequency {}".format(fprobe))

        idx_probe_start = zimt_TAS(Vdc, Vamp, fprobe, num_periods, Gen, tVG_path, tVG_name='tVG.txt')
        str_lst, sys_lst, path_lst, tVG_lst, tj_lst = [], [], [], [], []

        params = { "tVG_file": "tVG.txt"           
                    }
        
        params_str = ''
        for key, value in params.items():
            params_str += "-" + key + " " + str(value) + " "

        str_lst.append(params_str)
        sys_lst.append(system)
        path_lst.append(path2ZimT)
        tVG_lst.append('tVG.txt')
        tj_lst.append('tj.dat')
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_zimt, list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)

        p.close()
        p.join()

        if move_output: 
            store_output_in_folder(tVG_lst,store_folder,path2ZimT)
            store_output_in_folder(tj_lst,store_folder,path2ZimT)
        
    for idx , _ in enumerate(tj_lst):
        data_tj = pd.read_csv(Path(path2ZimT, store_folder, tj_lst[idx]), delim_whitespace=True)
        
        tV = data_tj['Vext'][idx_probe_start:].to_numpy()
        tJ = data_tj['Jext'][idx_probe_start:].to_numpy()
        t = data_tj['t'][idx_probe_start:].to_numpy()
        plt.plot(t, tV)
        plt.show()
        plt.plot(t, tJ)
        plt.show()
        Z, C = calcZC(tV, tJ, fprobe, equalCircuit='RCp')

        tZ.append(Z)
        tC.append(C)

    plt.plot(fprobe_span, tZ)
    plt.show()
if __name__ == "__main__":
    TAS()
