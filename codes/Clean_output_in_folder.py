#####################################################
#### Clean output in folder for SIMsalabim/ZimT #####
#####################################################
# by Vincent M. Le Corre
import os
from pathlib import Path
from VLC_useful_func import clean_up_output


curr_dir = os.getcwd()
# path2folder = 'Simulation_program/DDSuite_v403_OPV/SIMsalabim'
path2folder  = 'Simulation_program/DDSuite_v403_OPV/ZimT'

## Clean-up outputs from folder
clean_up_output('JV',Path(curr_dir) / path2folder)
clean_up_output('Var',Path(curr_dir) / path2folder)
clean_up_output('tj',Path(curr_dir) / path2folder)
clean_up_output('tVG',Path(curr_dir) / path2folder)
print('Ouput data was deleted from '+ str(Path(curr_dir) / path2folder))