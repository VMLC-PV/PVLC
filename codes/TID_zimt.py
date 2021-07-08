from pathlib import Path
import warnings, platform, os, multiprocessing, parmap
from pathlib import Path
from tVG_gen import zimt_TID
from VLC_useful_func import run_zimt


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

    zimt_TID(Vfill = 1,
            tfill = 5,
            fill_steps = 100,
            Vdrift = 0,
            tdrift = 1,
            Vamp = 10e-3,
            freq = 1,
            drift_steps = 100,
            Gen = 0,
            tVG_name='tVG.txt')

    str_lst, sys_lst, path_lst, tVG_lst = [], [], [], []

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
    p = multiprocessing.Pool(max_jobs)
    
    results = parmap.starmap(run_zimt, list(zip(str_lst,sys_lst,path_lst)), pm_pool=p, pm_processes=max_jobs,pm_pbar=True)

    p.close()
    p.join()



if __name__ == "__main__":
    TID()
