############################################################
############# Big simulation with SIMsalabim ###############
############################################################
# by Vincent M. Le Corre
# Package import
import subprocess,shutil,os,tqdm,parmap,multiprocessing,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy.optimize
from pathlib import Path
from threading import Thread, Lock, current_thread
from queue import Queue
from multiprocessing import Process, Queue
import sys, select, os
# package by VLC
from VLC_useful_func import *
from tVG_gen import *
import plot_settings_screen


def make_str_cmd(parameters):
    str2run = ''
    for par in parameters:
        if par['scale'] == 'int':
            str_val = ' {:.0f}'.format(get_random_value(par['min'],par['max'],par['scale']))
        else:
            str_val =  ' {:.3e}'.format(get_random_value(par['min'],par['max'],par['scale']))
        str2run = str2run + ' -' + par['name'] + str_val
    return str2run

def worker(q, lock,parameters,System,path):
    global numIt
    while True:
        value = q.get()  # blocks until the item is available
        
        # do stuff...
        with lock:
            # prevent printing at the same time with this lock
            print(f"in {current_thread().name} got {value}")
            numIt = numIt + 1
        JVname = ' -JV_file JV_'+str(int(numIt))+'.dat'
        run_SIMsalabim(make_str_cmd(parameters)+JVname,System,path)
        # ...

        # For each get(), a subsequent call to task_done() tells the queue
        # that the processing on this item is complete.
        # If all tasks are done, q.join() can unblock
        q.task_done()


def big_simu():

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

    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/AutoFit125/DDSuite_v407_Xiaoyan/SIMsalabim'

    ## Define parameters
    parameters = [] #list with parameter input
    parameters.append({'name' : 'L', 'min' : 140e-9, 'max' : 150e-9 , 'scale' : 'lin'})
    parameters.append({'name' : 'L_RTL', 'min' : 10e-9, 'max' : 15e-9 , 'scale' : 'lin'})
    parameters.append({'name' : 'L_LTL', 'min' : 30e-9, 'max' : 35e-9 , 'scale' : 'lin'})
    parameters.append({'name' : 'kdirect', 'min' : 1e-18, 'max' : 1e-17 , 'scale' : 'log'})
    # parameters.append({'name' : 'Tr_type_B', 'min' : -1, 'max' : 1 , 'scale' : 'int'})

    itScan = 50
    
    q = Queue()
    lock = Lock()

    for i in range(max_jobs):
        t = Thread(name=f"Thread{i+1}", target=worker, args=(q, lock,parameters,system,path2SIMsalabim))
        t.daemon = True  # dies when the main thread dies
        t.start()

    
    # t2 = 
    idx = 0
    

    # i = 0
    # print( "I'm doing stuff. Press Enter to stop me!")
    # while idx < itScan :
    #     os.system('cls' if os.name == 'nt' else 'clear')
        
    #     # print(i)
    #     q.put(idx)

    #     idx = idx + 1
    #     print( "I'm doing stuff. Press Enter to stop me!")
    #     print(idx)
    #     if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
    #         # line = raw_input()
    #         q.join()
    #         break
    # i += 1
    
    while idx < itScan:
    
        q.put(idx)

        idx = idx + 1
    
    # q.empty()
    t.join()  # Blocks until all items in the queue have been gotten and processed.

    print(numIt)
    print('main done')


if __name__ == '__main__':
    numIt = 0
    
    big_simu()




