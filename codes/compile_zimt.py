###################################################
################# Compile zimt ####################
###################################################
# by Vincent M. Le Corre
# Package import
import os
import sys
import subprocess
import warnings
from pathlib import Path
# Don't show warnings
warnings.filterwarnings("ignore")

slash = '/'
curr_dir = os.getcwd()                      # Current working directory
path2zimt = 'Simulation_program/ZimT043_BETA'+slash    # Path to SIMsalabim in curr_dir


os.chdir(Path(curr_dir+slash+path2zimt+slash+'Units'))
subprocess.call('fpc InputOutputUtils.pp')
subprocess.call('fpc TypesAndConstants.pp')
subprocess.call('fpc NumericalUtils.pp')

os.chdir(Path(curr_dir+slash+path2zimt))
subprocess.call('fpc zimt.pas')
os.chdir(Path(curr_dir))
