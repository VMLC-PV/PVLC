###################################################
########## Make plots for SIMsalabim ###############
###################################################
# by Vincent M. Le Corre
# Package import
import os,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# Import homemade package by VLC
from VLC_units.plots.SimSS_plots import *
from VLC_units.simu.runsim import *
from VLC_units.simu.get_input_par import *
from run_plot_SIMsalabim import *


def compare_exp2fit():

    # path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/SIMsalabimv429_Xioyan/SimSS')
    path2SIMsalabim = os.path.join(os.getcwd() , 'Simulation_program/DDSuite_v409_crack/SIMsalabim')

    #### Define experiments to compare to simulation
    ## DT-Y6
    ID = 3918
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.008774872893385296, 0.0295868067205203, 0.09674039280460424, 0.2903605440421194, 0.9438975920716443,1, 2.2588639120448035]
    
    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 3.063E-5 -Rshunt 7.995E-1 -Gehp 1.217E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.353 -Nc 5.645E+27 -W_L 4.060 -St_R 6.496E+9 -St_L 9.067E+8 -CB 4.04' #Best 1%

      
    age = 'aged'
    pixel = 3
    Gfracs = [0,0.008125240951214394, 0.02704509630728539, 0.08851876613458132, 0.26551237818238604, 0.861085219869508,1, 2.051557459849598]
    
    fixed_str = '-NP 1000 -mun_0 1.107E-7 -mup_0 1.269E-7 -kdirect 3.483E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.412E+28 -Etrap 4.664 -Bulk_tr 4.910E+19 -W_R 5.304 -Nc 5.519E+27 -W_L 4.119 -St_R 6.994E+8 -St_L 1.344E+6 -CB 4.04 -Vmax 1'

    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.412E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.204 -Nc 5.645E+27 -W_L 4.219 -St_R 6.496E+15 -St_L 9.067E+15 -CB 4.04' 
    ## DT-Y6
    # ID = 3918
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.008774872893385296, 0.0295868067205203, 0.09674039280460424, 0.2903605440421194, 0.9438975920716443,1, 2.2588639120448035]
    
    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 3.063E-5 -Rshunt 7.995E-1 -Gehp 1.217E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.353 -Nc 5.645E+27 -W_L 4.060 -St_R 6.496E+9 -St_L 9.067E+8 -CB 4.04' #Best 1%

      
    # age = 'aged'
    # pixel = 3
    # Gfracs = [0,0.008125240951214394, 0.02704509630728539, 0.08851876613458132, 0.26551237818238604, 0.861085219869508,1, 2.051557459849598]
    
    # # fixed_str = '-mun_0 1.107E-7 -mup_0 1.269E-7 -kdirect 3.483E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.412E+28 -Etrap 4.664 -Bulk_tr 4.910E+19 -W_R 5.304 -Nc 5.519E+27 -W_L 4.119 -St_R 6.994E+8 -St_L 1.344E+6 -CB 4.04'

    # fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 6.487E-5 -Rshunt 6.844E-1 -Gehp 1.412E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.204 -Nc 5.645E+27 -W_L 4.219 -St_R 6.496E+15 -St_L 9.067E+15 -CB 4.04' 


    # pixel = 1
    # Gfracs = [0,0.00866783401714054, 0.029238534889572412, 0.09543346872544842, 0.2862747313125504, 0.9306521928901167,1, 2.226981200864729]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    ## BTP-4F-12
    # ID = 3937
    # age = 'fresh'
    # pixel = 2
    # Gfracs = [0,0.0090015933301563, 0.029796718028699196, 0.09717869847550556, 0.29087180672086876, 0.943525832647383, 1, 2.2428264040319434]


    # fixed_str = '-mun_0 2.312E-8 -mup_0 2.247E-7 -kdirect 8.000E-18 -Rseries 5.468E-5 -Rshunt 2.521 -Gehp 1.311E+28 -Etrap 4.587 -Bulk_tr 4.894E+18 -W_R 5.367 -Nc 5.224E+27 -W_L 4.110 -St_R 3.152E+6 -St_L 4.313E+8 -CB 4.06' #keep'
    # Parameters to vary
    
    # pixel = 3
    # Gfracs = [0,0.009016735212175199, 0.029660936641206906, 0.0968623846255899, 0.29014650417906485, 0.9410713947273658, 1, 2.241490248848625]
    # fixed_str = '' # can chose a custom string here (advice don't put the thicknesses here but in the parameters below)

    # age = 'aged'
    # pixel = 2
    # Gfracs = [0,0.007814523992331497, 0.026764391395073882, 0.08785804792342802, 0.26391800166747525, 0.8572188961690461, 1,2.042784401098461]
 
    # fixed_str = '-mun_0 5.471E-8 -mup_0 1.734E-8 -kdirect 7.703E-18 -Rseries 1.092E-5 -Rshunt 1.241 -Gehp 1.384E+28 -Etrap 4.626 -Bulk_tr 6.462E+18 -W_R 5.442 -Nc 6.311E+27 -W_L 4.156 -St_R 2.733E+5 -St_L 4.750E+10'


    ## Y6-BO-4F
    # ID = 3916
    # age = 'fresh'
    # pixel = 4
    # Gfracs = [0,0.00899285762350159, 0.02958847222355735, 0.09647014717333, 0.2897613119670855, 0.9399050246570602, 1, 2.2396445153662032]

    # fixed_str = '-mun_0 2.442E-8 -mup_0 1.829E-7 -kdirect 7.999E-18 -Rseries 2.652E-5 -Rshunt 1.370 -Gehp 1.326E+28 -Etrap 4.558 -Bulk_tr 6.827E+18 -W_R 5.485 -Nc 5.799E+27 -W_L 4.074 -St_R 1.336E+8 -St_L 9.698E+9 -CB 4.05'

    # pixel = 2
    # Gfracs = [0,0.009009535782328836, 0.029625731507062908, 0.09649315148591558, 0.2897571956906255, 0.9400858561495319,1, 2.238969531299488]
    # fixed_str = '-mun_0 7.002E-8 -mup_0 2.733E-8 -kdirect 8.000E-18 -Rseries 1.182E-5 -Rshunt 1.560 -Gehp 1.320E+28 -Etrap 4.683 -Bulk_tr 6.336E+18 -W_R 5.424 -Nc 5.739E+27 -W_L 4.085 -St_R 2.074E+5 -St_L 6.013E+4 -CB 4.05'


    # o-IDFBR
    # ID = 3939
    # age = 'fresh'
    # pixel = 3
    # Gfracs = [0,0.007703621989336275, 0.030391616105901822, 0.10494576208861922, 0.31542562971134397, 0.994778451921309, 1, 2.2770546056260343]
    # fixed_str = '-mun_0 1.678E-7 -mup_0 8.802E-7 -kdirect 2.816E-18 -Rseries 1.591E-6 -Rshunt 1.210 -Gehp 3.657E+27 -Etrap 4.319 -Bulk_tr 3.476E+20 -W_R 5.152 -Nc 1.120E+26 -W_L 3.826 -St_R 8.787E+5 -St_L 8.865E+4 -CB 3.7'
    
    parameters = []
    parameters.append({'name':'Gfrac','values':Gfracs})

    JVexp_lst= ['PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_dark.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int3.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int10.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int33.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int100.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int330.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_am15.txt','PM6_'+str(ID)+'_px_'+str(pixel)+'_'+age+'_int800.txt']

    run_plot_SIMsalabim(fixed_str = fixed_str, parameters = parameters,path2SIMsalabim = path2SIMsalabim , run_simu = True, plot_JVs = True, plot_nrj_diag = False, plot_densities = False, plot_exp = True, JVexp_lst = JVexp_lst,  verbose = True)
     
  


if __name__ == '__main__':
    
    compare_exp2fit()
    plt.show()