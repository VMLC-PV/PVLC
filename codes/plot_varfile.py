###################################################
########## Make plots from Var_file ###############
###################################################
# by Vincent M. Le Corre
# Package import
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plot_settings_screen
from scipy import stats
import subprocess,shutil,os
from VLC_useful_func import sci_notation,run_SIMsalabim
import warnings

# Don't show warnings
warnings.filterwarnings("ignore")
Vt = 0.0256 # at 25 C

# Inputs
sytem = 'Windows'
size_fig = (16, 12)
num_fig = 0
run_simu = True
Simu_str = '-Vmin 0 -Vmax 0 -UseExpData 0'
color_nTL = '#c7e6a3'
color_pTL = '#8cb5f0'
color_pero = '#ba0000'
color_electrode ='#999999'

# Run simulation

if run_simu:
    run_SIMsalabim(Simu_str,sytem)

# Import Var_file
path = ''
names = ['x','V','n','p','Evac','Ec','Ev','phin','phip','ntrap','ptrap','nid','pid','nion','pion','mun','mup','rec','dp','Gm','Jn','Jp']
Var_file_name = 'Var.dat'
data_var = pd.read_csv(path+Var_file_name,delim_whitespace=True)

data_var['x'] =data_var['x'] * 1e9
## Energy diagram plot
plot_nrj_diag =True
if plot_nrj_diag:
    line_thick = 3
    th_TL_left = 20
    th_TL_right = 20
    num_fig = num_fig+1
    f_nrj_diag = plt.figure(num_fig,figsize=size_fig)
    ax_nrj_diag = plt.axes()
    # ax_nrj_diag.plot('x','Evac',data=data_var,label = r'E$_{vac}$',linestyle='-',linewidth=2,color = 'k')
    ax_nrj_diag.plot('x','Ec',data=data_var,label = r'E$_{c}$',linestyle='-', linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','Ev',data=data_var,label = r'E$_{v}$',linestyle='-', linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','phin',data=data_var,label = r'E$_{fn}$',linestyle='--',linewidth=line_thick,color = 'k')
    ax_nrj_diag.plot('x','phip',data=data_var,label = r'E$_{fp}$',linestyle='--',linewidth=line_thick,color = 'k')

    TL_left = data_var[data_var['x']<th_TL_left]
    TL_right = data_var[data_var['x']>max(data_var['x'])-th_TL_right]
    AL = data_var[data_var['x']<max(data_var['x'])-th_TL_right]
    AL = AL[AL['x']>th_TL_left]
    ax_nrj_diag.fill_between(TL_left['x'],TL_left['Ec'],y2=0,color=color_nTL)
    ax_nrj_diag.fill_between(TL_left['x'],TL_left['Ev'],y2=-8,color=color_nTL)
    ax_nrj_diag.fill_between(TL_right['x'],TL_right['Ec'],y2=0,color=color_pTL)
    ax_nrj_diag.fill_between(TL_right['x'],TL_right['Ev'],y2=-8,color=color_pTL)
    ax_nrj_diag.fill_between(AL['x'],AL['Ec'],y2=0,color=color_pero)
    ax_nrj_diag.fill_between(AL['x'],AL['Ev'],y2=-8,color=color_pero)
    ax_nrj_diag.plot([-10,0],[min(data_var['phin']),min(data_var['phin'])],color='k')
    ax_nrj_diag.plot([max(data_var['x']),max(data_var['x'])+10],[max(data_var['phip']),max(data_var['phip'])],color='k')
    ax_nrj_diag.fill_between([-10,0],[min(data_var['phin']),min(data_var['phin'])],y2=-8,color=color_electrode)
    ax_nrj_diag.fill_between([max(data_var['x']),max(data_var['x'])+10],[max(data_var['phip']),max(data_var['phip'])],y2=-8,color=color_electrode)
    # plt.axhline(y=max(data_var['phip']), color='k', linestyle='-')
    # Hide axis and spines
    ax_nrj_diag.get_xaxis().set_visible(False)
    ax_nrj_diag.get_yaxis().set_visible(False)
    for sides in ['right','left','top','bottom']:
        ax_nrj_diag.spines[sides].set_visible(False)
    # legend
    plt.legend(loc='center',frameon=False,ncol = 2, bbox_to_anchor=(0.52,1.02),fontsize = 40)
    plt.tight_layout()
    plt.savefig('Energy_diagram.eps',dpi=600,transparent=True)

plt.show()

