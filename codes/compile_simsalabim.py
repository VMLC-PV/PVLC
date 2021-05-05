###################################################
######### Test SIMsalabim versus scpas ############
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
path2SIMsalabim = 'Simulation_program/SIMsalabimv385'+slash    # Path to SIMsalabim in curr_dir


os.chdir(Path(curr_dir+slash+path2SIMsalabim+slash+'Units'))
subprocess.call('fpc InputOutputUtils.pp')
subprocess.call('fpc TypesAndConstants.pp')
subprocess.call('fpc NumericalUtils.pp')

os.chdir(Path(curr_dir+slash+path2SIMsalabim))
subprocess.call('fpc SIMsalabim.pas')
os.chdir(Path(curr_dir))
