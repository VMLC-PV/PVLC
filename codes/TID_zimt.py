from pathlib import Path
import warnings, platform, os, multiprocessing, parmap
from pathlib import Path
from tVG_gen import zimt_TID
from VLC_useful_func import run_zimt, store_output_in_folder
import pandas as pd
import matplotlib.pyplot as plt
import math


def TID():
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
    move_output = False
    store_folder = Path('TID')

    # Figures control
    size_fig = (16, 12)
    num_fig_tjs= 0
    f_tjs = plt.figure(num_fig_tjs,figsize=size_fig)

    Vfill = 1
    tfill = 5
    fill_steps = 100
    Vdrift = 0
    tdrift = 2
    dt_probe = 50e-2 # time between probing points 
    Vamp = 10e-3
    freq = 10e3
    drift_steps = 21
    Gen = 0
    tVG_path = path2ZimT
    tVG_name='tVG.txt'

    if not dt_probe == 0:
        zimt_TID(Vfill, tfill, fill_steps, Vdrift, tdrift, dt_probe, Vamp, freq, drift_steps, Gen, tVG_path, tVG_name)
    else:
        zimt_TID(Vfill, tfill, fill_steps, Vdrift, tdrift, dt_probe, Vamp, freq, drift_steps, Gen, tVG_path, tVG_name)



    str_lst, sys_lst, path_lst, tVG_lst, tj_lst = [], [], [], [], []

    params = {"Rin": 0, 
                "L": 1.4e-7,
                "W_L": 5,
                "W_R": 5, 
                "tVG_file": "tVG.txt"           
                }
    
    params_str = ''
    for key, value in params.items():
        params_str += "-" + key + " " + str(value) + " "

    str_lst.append(params_str)
    sys_lst.append(system)
    path_lst.append(path2ZimT)
    tVG_lst.append('tVG.txt')
    tj_lst.append('tj.dat')
    # p = multiprocessing.Pool(max_jobs)
    
    # results = parmap.starmap(run_zimt, list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)

    # p.close()
    # p.join()

    if move_output: 
        store_output_in_folder(tVG_lst,store_folder,path2ZimT)
        store_output_in_folder(tj_lst,store_folder,path2ZimT)

    for idx , _ in enumerate(tj_lst):
        data_tj = pd.read_csv(Path(path2ZimT, store_folder, tj_lst[idx]), delim_whitespace=True)
        idx_drift = data_tj[data_tj["t"] == tfill].index.values[0] + 1
        tV_drift = data_tj['Vext'][idx_drift:].to_numpy()
        tJ_drift = data_tj['Jext'][idx_drift:].to_numpy()
        t_drift = data_tj['t'][idx_drift:].to_numpy()

        period = 1/freq
        num_periods = math.floor((t_drift[-1]-t_drift[0])/period)
        

        # for p in range(0, num_periods, 1000):
        #     idx_start = p*drift_steps
        #     idx_end = idx_start + 2*drift_steps
        #     print(tV_drift[idx_start])
    

if __name__ == "__main__":
    TID()