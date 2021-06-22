###################################################
############### Compile SIMsalabim ################
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,subprocess,warnings,platform
from pathlib import Path


def compile_zimt():
    # Don't show warnings
    warnings.filterwarnings("ignore")
    System = platform.system()                  # Operating system
    curr_dir =  os.getcwd()                     # Current working directory
    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/DDSuite_v403_OPV/ZimT' # Path to SIMsalabim in curr_dir


    os.chdir(path2SIMsalabim) # Go to SIMsalabim directory

    # Compile depending on the operating system
    if System == 'Windows':
        subprocess.call('fpc zimt.pas')
    elif System == 'Linux':
        subprocess.check_call(('fpc zimt.pas').split())
    else: print('Wrong system input')

    os.chdir(Path(curr_dir)) # Come back to current directory


if __name__ == '__main__':

    compile_zimt()