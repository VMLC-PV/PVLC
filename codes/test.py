###################################################
########## Make plots for SIMsalabim ###############
###################################################
# by Vincent M. Le Corre
# Package import
import subprocess,shutil,os,tqdm,parmap,multiprocessing,platform,warnings,itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from pathlib import Path
# package by VLC
from VLC_useful_func import *
import plot_settings_screen

def GetParFromStr(str2run):
    """Get parameters from command string for SIMsalabim or ZimT

    Parameters
    ----------
    str2run : STR
        Command string for SIMsalabim or ZimT

    Returns
    -------
    dict
        Contains the parameters and values from the command string
    """    

    str2run = ' '.join(str2run.split()) #remove extra white space

    str2run = str2run.split()

    Names= str2run[::2]
    Values = str2run[1::2]
    ParStrDic = {}
    for i,j in enumerate(Names):
        Names[i] = Names[i].replace('-', '')
        Values[i] = float(Values[i])
        ParStrDic[Names[i]] = float(Values[i])
    # pd.DataFrame({'Parameters':Names,'Values':Values})

    return ParStrDic

def ReadParameterFile(path2file):
    """Get all the parameters from the 'Device_parameters.txt' file
    for SIMsalabim and ZimT
    Parameters
    ----------
    path2file : str
        Path to the 'Device_parameters.txt'

    Returns
    -------
    dict
        Contains the parameters and values from the 'Device_parameters.txt'
    """    
    
    lines = []
    ParFileDic = {}
    with open(path2file) as f:
        lines = f.readlines()

    count = 0
    for line in lines:
        line = line.replace(' ', '')
        if line[0] != '*' and (not line.isspace()):
            equal_idx = line.find('=')
            star_idx = line.find('*')
            # print(line[0:equal_idx] , line[equal_idx+1:star_idx])
            ParFileDic[line[0:equal_idx] ] = line[equal_idx+1:star_idx]
            count += 1
            # print(f'line {count}: {line}')   
    return ParFileDic

def ChosePar(parname,ParStrDic,ParFileDic):
    if parname in ParStrDic.keys():
        parval = ParStrDic[parname]
    else :
        parval = ParFileDic[parname]
    
    return parval

def CalcVLCFOM(str2run,path2DevFile):

    q = constants.value(u'elementary charge')
    eps_0 = constants.value(u'vacuum electric permittivity')

    ParStrDic = GetParFromStr(str2run)
    ParFileDic = ReadParameterFile(path2DevFile)

    Nc = float(ChosePar('Nc',ParStrDic,ParFileDic))
    L = float(ChosePar('L',ParStrDic,ParFileDic))
    eps_r = float(ChosePar('eps_r',ParStrDic,ParFileDic))
    Egap = abs(float(ChosePar('CB',ParStrDic,ParFileDic))-float(ChosePar('VB',ParStrDic,ParFileDic)))
    
    mun_0 = float(ChosePar('mun_0',ParStrDic,ParFileDic))
    mup_0 = float(ChosePar('mup_0',ParStrDic,ParFileDic))
    Gehp = float(ChosePar('Gehp',ParStrDic,ParFileDic))
    Bulk_tr = float(ChosePar('Bulk_tr',ParStrDic,ParFileDic))
    Cn = float(ChosePar('Cn',ParStrDic,ParFileDic))
    Cp = float(ChosePar('Cp',ParStrDic,ParFileDic))

    L_LTL = float(ChosePar('L_LTL',ParStrDic,ParFileDic))
    L_RTL = float(ChosePar('L_RTL',ParStrDic,ParFileDic))

    if int(ChosePar('UseLangevin',ParStrDic,ParFileDic)) == 1:
        Lang_pre = float(ChosePar('Lang_pre',ParStrDic,ParFileDic))
        gamma = Lang_pre*(q*(mun_0+mup_0))/(eps_r*eps_0)
    else:
        gamma = float(ChosePar('kdirect',ParStrDic,ParFileDic))


    Vint = Egap - 0.4
    Ceff = np.sqrt(Cn*Cp)
    mueff = np.sqrt(mun_0*mup_0)

    delta_b = (gamma*Nc**2)/Gehp
    Delta_b = (gamma*Gehp*(L-L_LTL-L_RTL)**4)/(mun_0*mup_0*Vint**2)
    delta_tr = (Cn * Bulk_tr * Nc**2)/(Gehp**2)    
    Delta_tr = (Ceff * Bulk_tr * (L-L_LTL-L_RTL)**2)/(mueff * Vint)

    return [delta_b,Delta_b,delta_tr,Delta_tr]





   



def test():

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

    path2SIMsalabim = Path(os.getcwd()) /'Simulation_program/DDSuite_v409_crack/SIMsalabim'

    fixed_str = '-mun_0 7.020E-8 -mup_0 4.315E-8 -kdirect 7.245E-18 -Rseries 3.063E-5 -Rshunt 7.995E-1 -Gehp 1.217E+28 -Etrap 4.774 -Bulk_tr 3.477E+18 -W_R 5.353 -Nc 5.645E+27 -W_L 4.060 -St_R 6.496E+9 -St_L 9.067E+8 -CB 4.04' #Best 1%

    # print(GetParFromStr(fixed_str))

    # print(ReadParameterFile(Path(path2SIMsalabim / 'Device_parameters.txt')))

    print(CalcVLCFOM(fixed_str,Path(path2SIMsalabim / 'Device_parameters.txt')))



if __name__ == '__main__':
    
    test()
    plt.show()