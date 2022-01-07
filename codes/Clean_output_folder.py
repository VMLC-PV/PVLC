
######################################################################
####################### Cleaning output files ########################
######################################################################

# Import homemade package by VLC
from VLC_units.cleanup_folder.clean_folder import *

# folder to clean
folder2clean = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabim_v427/SimSS') 


clean_up_output('tj',folder2clean)
clean_up_output('tVG',folder2clean)
clean_up_output('JV',folder2clean)
clean_up_output('Var',folder2clean)
clean_up_output('scPars',folder2clean)


