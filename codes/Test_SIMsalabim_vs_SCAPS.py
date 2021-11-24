###################################################
######### Test SIMsalabim versus scaps ############
###################################################
# by Vincent M. Le Corre
# Package import
import os,sys,platform,tqdm,parmap,multiprocessing,platform,subprocess,shutil,warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.lines as mlines
from time import time
from scipy import constants
from pathlib import Path
from scipy import stats
# Package by VLC
import plot_settings_screen
from VLC_useful_func import *

# Don't show warnings
warnings.filterwarnings("ignore")

# Main Program
def SIMsalabim_vs_scaps(path_scaps, str_lst, has_TL, labels, JV_files, Var_files,path2SIMsalabim,Store_Folder, JVscaps_lst, NPexp_lst, JV_plot_filename, np_plot_filename):
    """Run a comparison between pre-ran drift-diffusion simulation made using scaps (http://scaps.elis.ugent.be/)
    and the results from SIMsalabim (https://github.com/kostergroup/SIMsalabim)
    All the different case scenarios are described in the 'Test Configurations.pptx' or .pdf attached in the 'Test simulation' Folder.

    Parameters
    ----------
    path_scaps : list of str
        list containing the path to the scaps simulation results

    str_lst : list of str
        list of the string to run for SIMsalabim

    has_TL : bool
        If True the test include transport layers

    labels : list of str
        list containing the labels for the plot

    JV_files : list of str
        list containing the filenames for the JV_file output from SIMsalabim

    Var_files : list of str
        list containing the filenames for the Var_file output from SIMsalabim

    path2SIMsalabim : str
        string with the path to SIMsalabim

    Store_Folder : str
        name of the folder where the SIMsalabim output it stored 

    JVscaps_lst : list of str
        list containing the filenames of the output from scaps containing the JV data

    NPexp_lst : list of str
        list containing the filenames of the output from scaps containing the carrier density data

    JV_plot_filename : str
        filename for the saved image of the JV comparison figure
        
    np_plot_filename : str
        filename for the saved image of the carrier density comparison figure
    """    
    # General Inputs
    system = platform.system()                  # Operating system
    # Max number of parallel simulations (for number of CPU use: os.cpu_count() )
    max_jobs = os.cpu_count()-2
    if system == 'Windows':             # cannot easily do multiprocessing in Windows
        max_jobs = 1
        slash = '/'
        try:
            os.system('taskkill.exe /F /IM SIMsalabim.exe')
        except:
            pass
    else:
        slash = '/'

    # Simulation input
    run_simu = True                                         # Rerun simu?
 

    # Initialize
    curr_dir = os.getcwd()
    idx = 0
    start = time()
    lines = ['-', '--', '-.', ':']
    sys_lst, path_lst = [], []
    for i in str_lst:
        sys_lst.append(system)
        path_lst.append(curr_dir+slash+path2SIMsalabim)
        
    # Figures control
    size_fig = (16, 12)
    f_JVs = plt.figure(0, figsize=size_fig)
    f_dens = plt.figure(1, figsize=size_fig)
    density2plot = ['n', 'p']
    dens_patches = [mlines.Line2D([], [], color='k', linestyle=lines[0], marker='None', label='Electron'), mlines.Line2D([], [], color='k', linestyle=lines[1], marker='None', label='Hole')]

    # Color range for plotting
    colors = plt.cm.viridis(np.linspace(0, 1, max(len(str_lst), 4) + 1))

    # Run SIMsalabim
    if run_simu:
        p = multiprocessing.Pool(max_jobs)
        results = parmap.starmap(run_SIMsalabim, list(
            zip(str_lst, sys_lst, path_lst)), pm_pool=p, pm_processes=max_jobs, pm_pbar=True)
        p.close()
        p.join()
    ## Move output folder to new folder
    # Create directory if it does not exist
    if not os.path.exists(Store_Folder):
        os.makedirs(Store_Folder)
    # move file into the new folder
    for i in JV_files:
        if os.path.exists(os.path.join(curr_dir,path2SIMsalabim,i)):
            os.replace(os.path.join(curr_dir,path2SIMsalabim,i),Store_Folder/i)
    for i in Var_files:
        if os.path.exists(os.path.join(curr_dir,path2SIMsalabim,i)):
            os.replace(os.path.join(curr_dir,path2SIMsalabim,i),Store_Folder/i)
    # Start file import and plotting
    for Simu_str, exp_name, lbl, JV_file_name, Var_file_name, npexp_name in zip(str_lst, JVscaps_lst, labels, JV_files, Var_files, NPexp_lst):

        # Import files
        data_JV = make_df_JV(Path(Store_Folder / JV_file_name))
        if has_TL:
            names_scaps = ['v(V)', 'jtot(mA/cm2)', 'j_total_rec(mA/cm2)', 'j_total_gen(mA/cm2)', '   jbulk(mA/cm2)', 'jifr(mA/cm2)',
                           'jminor_left(mA/cm2)', 'jminor_right(mA/cm2)', 'j_SRH(mA/cm2)', 'j_Radiative(mA/cm2)', 'j_Auger(mA/cm2)']
            data_JVscaps = pd.read_csv(exp_name, names=names_scaps, engine="python", header=None, delim_whitespace=True,
                                       usecols=[0, 1], skiprows=30, skipfooter=10)
            names_densscaps = ['i', 'x(um)', 'y', 'Ec(eV)', 'Fn(eV)', 'Fp(eV)', 'Ev(eV)', 'n(/cm3)', 'p(/cm3)', 'rho(defect) (/cm3)', 'net doping (/cm3)', 'rho(/cm3)', 'E(V/cm)', 'jn(mA/cm2)', 'jp(mA/cm2)', 'jn_tunnel(mA/cm2)',
                               'jp_tunnel(mA/cm2)', 'jtot(mA/cm2)', 'generation(#/cm3.s)', 'recombination(#/cm3.s)', 'cumulative generation (mA/cm2)', 'cumulative recombination L to R (mA/cm2)', 'cumulative recombination R to L (mA/cm2)']
            data_densscaps = pd.read_csv(npexp_name, names=names_densscaps, engine="python",
                                         header=None, delim_whitespace=True, usecols=[1, 7, 8], skiprows=34, skipfooter=10)
        else:
            names_scaps = ['v(V)', 'jtot(mA/cm2)', 'j_total_rec(mA/cm2)', 'j_total_gen(mA/cm2)', '   jbulk(mA/cm2)', 'jifr(mA/cm2)',
                           'jminor_left(mA/cm2)', 'jminor_right(mA/cm2)', 'j_SRH(mA/cm2)', 'j_Radiative(mA/cm2)', 'j_Auger(mA/cm2)']
            data_JVscaps = pd.read_csv(exp_name, names=names_scaps, engine="python", header=None,
                                       delim_whitespace=True, usecols=[0, 1], skiprows=27, skipfooter=10)
            names_densscaps = ['i', 'x(um)', 'y', 'Ec(eV)', 'Fn(eV)', 'Fp(eV)', 'Ev(eV)', 'n(/cm3)', 'p(/cm3)', 'rho(defect) (/cm3)', 'net doping (/cm3)', 'rho(/cm3)', 'E(V/cm)', 'jn(mA/cm2)', 'jp(mA/cm2)', 'jn_tunnel(mA/cm2)',
                               'jp_tunnel(mA/cm2)', 'jtot(mA/cm2)', 'generation(#/cm3.s)', 'recombination(#/cm3.s)', 'cumulative generation (mA/cm2)', 'cumulative recombination L to R (mA/cm2)', 'cumulative recombination R to L (mA/cm2)']
            data_densscaps = pd.read_csv(npexp_name, names=names_densscaps, engine="python",
                                         header=None, delim_whitespace=True, usecols=[1, 7, 8], skiprows=34, skipfooter=10)

        data_var = make_df_Var(Path(Store_Folder / Var_file_name))
        # data_var['x'] = data_var['x'] * 1e9

        ########################################################
        ## Compare SIMsalabim and scaps JVs
        plt.figure(0)
        plt.plot(data_JVscaps['v(V)'], data_JVscaps['jtot(mA/cm2)'], 'o',
                markeredgecolor=colors[idx], markersize=10, markerfacecolor='None', markeredgewidth=3)
        SIMsalabim_JVs_plot(0, data_JV, x='Vext', y=['Jext'], xlimits=[], ylimits=[
        ], plot_type=0, labels=lbl, colors=colors[idx], legend=True, save_yes=True, pic_save_name=Store_Folder / JV_plot_filename)
        plt.legend(prop={"size": 20})

        ## Compare SIMsalabim and scaps electron and hole densities at Voc
        plt.figure(1)
        plt.semilogy(1000*data_densscaps['x(um)'], data_densscaps['n(/cm3)'], 'o',
                        markeredgecolor=colors[idx], markersize=5, markerfacecolor='None', markeredgewidth=3)
        plt.semilogy(1000*data_densscaps['x(um)'], data_densscaps['p(/cm3)'], 'o',
                        markeredgecolor=colors[idx], markersize=5, markerfacecolor='None', markeredgewidth=3)
        dens_patches.append(mlines.Line2D(
            [], [], color=colors[idx], linestyle='-', marker='None', label=lbl))
        plt.legend(handles=dens_patches, loc='best',
                    frameon=False, fontsize=30)
        SIMsalabim_dens_plot(1, data_var,y=['n', 'p'], xlimits=[], ylimits=[], plot_type=2, labels=lbl, colors=colors[idx],
                                line_type=lines, legend=False, save_yes=True, pic_save_name=Store_Folder / np_plot_filename)
        plt.legend(prop={"size": 15}, ncol=3, )
        plt.ylim(1e8, 1e19)

        idx = idx+1

    plt.show()


if __name__ == '__main__':

    # Inputs
    slash = '/'
    curr_dir = os.getcwd()                      # Current working directory
    path2SIMsalabim = 'Simulation_program/DDSuite_v418/SIMsalabim'+slash    # Path to SIMsalabim in curr_dir
    ext_save_pic = '.jpg'
    
    # Simulation types
    MIM_configuration = False
    Pin_nip_configuration = False
    Traps_configuration = False
    Interface_traps_configuration = True

    # First setup for MIM
    if MIM_configuration:
        print('\n')
        print('Start the MIM configuration comparison:')
        path_scaps = Path(curr_dir+slash+'Test simulation/MIM Configuration')
        Store_Folder = Path(curr_dir+slash+'Test simulation/MIM Configuration')
        str_lst = ['-Nc 5e24 -L 300e-9 -JV_file JV_MIM_ref.dat -Var_file Var_MIM_ref.dat',
                   '-Nc 5e24 -kdirect 1e-17 -JV_file JV_MIM_1.dat -Var_file Var_MIM_1.dat',
                   '-Nc 5e24 -mun_0 1e-7 -mup_0 1e-7 -JV_file JV_MIM_2.dat -Var_file Var_MIM_2.dat',
                   '-Nc 5e24 -mup_0 1e-7 -JV_file JV_MIM_3.dat -Var_file Var_MIM_3.dat',
                   '-Nc 5e24 -L 100E-9 -Gehp 1.3e28 -JV_file JV_MIM_4.dat -Var_file Var_MIM_4.dat',
                   '-Nc 5e24 -W_L 3.8 -W_R 5 -JV_file JV_MIM_5.dat -Var_file Var_MIM_5.dat',
                   '-Nc 5e24 -VB 4.8 -W_R 4.8 -JV_file JV_MIM_6.dat -Var_file Var_MIM_6.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 4', 'Test number 5',
                  'Test number 6']
        JV_files = ['JV_MIM_ref.dat', 'JV_MIM_1.dat','JV_MIM_2.dat','JV_MIM_3.dat','JV_MIM_4.dat','JV_MIM_5.dat','JV_MIM_6.dat']
        Var_files = ['Var_MIM_ref.dat', 'Var_MIM_1.dat','Var_MIM_2.dat','Var_MIM_3.dat','Var_MIM_4.dat','Var_MIM_5.dat','Var_MIM_6.dat']
        JVscaps_lst = [path_scaps / 'testref.iv', path_scaps / 'test1.iv', path_scaps / 'test2.iv', path_scaps / 'test3.iv',
                     path_scaps / 'test4.iv', path_scaps / 'test5.iv', path_scaps / 'test6.iv']
        NPexp_lst = [path_scaps / 'nptestref.eb', path_scaps / 'nptest1.eb', path_scaps / 'nptest2.eb', path_scaps / 'nptest3.eb',
                     path_scaps / 'nptest4.eb', path_scaps / 'nptest5.eb', path_scaps / 'nptest6.eb']
        JV_plot_filename = 'JV_MIM_configuration_SIMsalabim'+ext_save_pic 
        np_plot_filename = 'np_MIM_configuration_SIMsalabim'+ext_save_pic 
        has_TL = False
        SIMsalabim_vs_scaps(path_scaps, str_lst, has_TL, labels, JV_files, Var_files,path2SIMsalabim,Store_Folder, JVscaps_lst, NPexp_lst, JV_plot_filename, np_plot_filename)
        

    # Pin-nip simulation
    if Pin_nip_configuration:
        print('\n')
        print('Start the pin-nip configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Pin-nip Configuration')
        Store_Folder = Path(curr_dir + slash + 'Test simulation/Pin-nip Configuration')
        str_lst = ['-Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_TL_ref.dat -Var_file Var_TL_ref.dat',
                   '-Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -JV_file JV_TL_1.dat -Var_file Var_TL_1.dat',
                   '-Nc 1e24 -Gehp 4.2e27 -L 400e-9 -L_LTL 50e-9 -L_RTL 50e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -JV_file JV_TL_2.dat -Var_file Var_TL_2.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_TL_3.dat -Var_file Var_TL_3.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -JV_file JV_TL_4.dat -Var_file Var_TL_4.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -Nc 1e24 -Gehp 4.2e27 -L 400e-9 -L_LTL 50e-9 -L_RTL 50e-9 -mob_LTL 1e-7 -mob_RTL 1e-7 -JV_file JV_TL_5.dat -Var_file Var_TL_5.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2',
                  'Test number 3', 'Test number 4', 'Test number 5']
        JV_files = ['JV_TL_ref.dat','JV_TL_1.dat','JV_TL_2.dat','JV_TL_3.dat','JV_TL_4.dat','JV_TL_5.dat']
        Var_files = ['Var_TL_ref.dat','Var_TL_1.dat','Var_TL_2.dat','Var_TL_3.dat','Var_TL_4.dat','Var_TL_5.dat']
        JVscaps_lst = [path_scaps / 'TLtestref.iv', path_scaps / 'TLtest1.iv', path_scaps / 'TLtest2.iv',
                     path_scaps / 'TLtest3.iv', path_scaps / 'TLtest4.iv', path_scaps / 'TLtest5.iv']
        NPexp_lst = [path_scaps / 'npTLtestref.eb', path_scaps / 'npTLtest1.eb', path_scaps / 'npTLtest2.eb',
                     path_scaps / 'npTLtest3.eb', path_scaps / 'npTLtest4.eb', path_scaps / 'npTLtest5.eb']
        JV_plot_filename = 'JV_Pin_nip_configuration_SIMsalabim'+ext_save_pic 
        np_plot_filename = 'np_Pin_nip_configuration_SIMsalabim'+ext_save_pic 
        has_TL = True
        SIMsalabim_vs_scaps(path_scaps, str_lst, has_TL, labels, JV_files, Var_files,path2SIMsalabim,Store_Folder, JVscaps_lst, NPexp_lst, JV_plot_filename, np_plot_filename)

    # Traps simulation
    if Traps_configuration:
        print('\n')
        print('Start the bulk traps configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Traps Configuration')
        Store_Folder = Path(curr_dir + slash + 'Test simulation/Traps Configuration')
        str_lst = ['-Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_ref.dat -Var_file Var_Traps_ref.dat -Tr_type_B 0',
                   '-Bulk_tr 1e20 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_1.dat -Var_file Var_Traps_1.dat -Tr_type_B 0',
                   '-Tr_type_B 1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_2.dat -Var_file Var_Traps_2.dat',
                   '-Tr_type_B 1 -Bulk_tr 1e20 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_3.dat -Var_file Var_Traps_3.dat',
                   '-Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_4.dat -Var_file Var_Traps_4.dat',
                   '-Cn 1e-12 -Cp 1e-12 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.68e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_5.dat -Var_file Var_Traps_5.dat',
                   '-Tr_type_B -1 -Bulk_tr 1e20 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_6.dat -Var_file Var_Traps_6.dat',
                   '-Etrap 3.7 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_7.dat -Var_file Var_Traps_7.dat',
                   '-Etrap 3.9 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_8.dat -Var_file Var_Traps_8.dat',
                   '-Etrap 5.1 -Tr_type_B -1 -Bulk_tr 1e21 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_9.dat -Var_file Var_Traps_9.dat',
                   '-Etrap 5.1 -Tr_type_B -1 -Bulk_tr 1e22 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_10.dat -Var_file Var_Traps_10.dat',
                   '-Etrap 5.1 -Tr_type_B 1 -Bulk_tr 1e22 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_Traps_11.dat -Var_file Var_Traps_11.dat', ]
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 4', 'Test number 5',
                  'Test number 6', 'Test number 7', 'Test number 8', 'Test number 9', 'Test number 10', 'Test number 11', ]
        JV_files = ['JV_Traps_ref.dat','JV_Traps_1.dat','JV_Traps_2.dat','JV_Traps_3.dat','JV_Traps_4.dat','JV_Traps_5.dat','JV_Traps_6.dat','JV_Traps_7.dat','JV_Traps_8.dat','JV_Traps_9.dat','JV_Traps_10.dat','JV_Traps_11.dat', ]
        Var_files = ['Var_Traps_ref.dat','Var_Traps_1.dat','Var_Traps_2.dat','Var_Traps_3.dat','Var_Traps_4.dat','Var_Traps_5.dat','Var_Traps_6.dat','Var_Traps_7.dat','Var_Traps_8.dat','Var_Traps_9.dat','Var_Traps_10.dat','Var_Traps_11.dat', ]
        JVscaps_lst = [path_scaps / 'traptestref.iv', path_scaps / 'traptest1.iv', path_scaps / 'traptest2.iv',
                     path_scaps / 'traptest3.iv',path_scaps / 'traptest4.iv', path_scaps / 'traptest5.iv',
                     path_scaps / 'traptest6.iv',path_scaps / 'traptest7.iv', path_scaps / 'traptest8.iv',
                     path_scaps / 'traptest9.iv', path_scaps / 'traptest10.iv', path_scaps / 'traptest11.iv']
        NPexp_lst = [path_scaps / 'nptraptestref.eb', path_scaps / 'nptraptest1.eb', path_scaps / 'nptraptest2.eb',
                     path_scaps / 'nptraptest3.eb',path_scaps / 'nptraptest4.eb', path_scaps / 'nptraptest5.eb',
                     path_scaps / 'nptraptest6.eb',path_scaps / 'nptraptest7.eb', path_scaps / 'nptraptest8.eb',
                     path_scaps / 'nptraptest9.eb', path_scaps / 'nptraptest10.eb', path_scaps / 'nptraptest11.eb']
        JV_plot_filename = 'JV_Traps_configuration_SIMsalabim'+ext_save_pic 
        np_plot_filename = 'np_Traps_configuration_SIMsalabim'+ext_save_pic 
        has_TL = True
        SIMsalabim_vs_scaps(path_scaps, str_lst, has_TL, labels, JV_files, Var_files,path2SIMsalabim,Store_Folder, JVscaps_lst, NPexp_lst, JV_plot_filename, np_plot_filename)

    # Interface traps simulation
    if Interface_traps_configuration:
        print('\n')
        print('Start the interface traps configuration comparison:')
        path_scaps = Path(curr_dir + slash + 'Test simulation/Interface Traps Configuration')
        Store_Folder = Path(curr_dir + slash + 'Test simulation/Interface Traps Configuration')
        str_lst = ['-St_L 1e14 -Tr_type_L 0 -St_R 1e14 -Tr_type_R 0 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_ref.dat -Var_file Var_InterfaceTraps_ref.dat',
                   '-St_L 1e13 -Tr_type_L 0 -St_R 1e13 -Tr_type_R 0 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_1.dat -Var_file Var_InterfaceTraps_1.dat',
                   '-St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_2.dat -Var_file Var_InterfaceTraps_2.dat',
                   '-St_L 1e14 -Tr_type_L -1 -St_R 1e12 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_3.dat -Var_file Var_InterfaceTraps_3.dat',
                   '-St_L 1e14 -Tr_type_L 1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_4.dat -Var_file Var_InterfaceTraps_4.dat',
                   '-Cn 1e-12 -Cp 1e-12 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_5.dat -Var_file Var_InterfaceTraps_5.dat',
                   '-St_L 1e14 -Tr_type_L 1 -St_R 1e14 -Tr_type_R -1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_6.dat -Var_file Var_InterfaceTraps_6.dat',
                   '-Etrap 3.7 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_7.dat -Var_file Var_InterfaceTraps_7.dat',
                   '-Etrap 3.9 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_8.dat -Var_file Var_InterfaceTraps_8.dat',
                   '-Etrap 5.1 -St_L 1e14 -Tr_type_L -1 -St_R 1e14 -Tr_type_R 1 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_9.dat -Var_file Var_InterfaceTraps_9.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -St_L 1e14 -Tr_type_L 0 -St_R 1e14 -Tr_type_R 0 -Nc 1e24 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_10.dat -Var_file Var_InterfaceTraps_10.dat',
                   '-Nc_LTL 1e26 -Nc_RTL 1e26 -St_L 1e13 -Tr_type_L 0 -St_R 1e13 -Tr_type_R 0 -Nc 1e24 -Gehp 4.2e27 -L 340e-9 -L_LTL 20e-9 -L_RTL 20e-9 -JV_file JV_InterfaceTraps_11.dat -Var_file Var_InterfaceTraps_11.dat']
        labels = ['Reference', 'Test number 1', 'Test number 2', 'Test number 3', 'Test number 4', 'Test number 5',
                  'Test number 6', 'Test number 7', 'Test number 8', 'Test number 9', 'Test number 10',
                  'Test number 11', ]
        JV_files = ['JV_InterfaceTraps_ref.dat','JV_InterfaceTraps_1.dat','JV_InterfaceTraps_2.dat','JV_InterfaceTraps_3.dat','JV_InterfaceTraps_4.dat','JV_InterfaceTraps_5.dat','JV_InterfaceTraps_6.dat','JV_InterfaceTraps_7.dat','JV_InterfaceTraps_8.dat','JV_InterfaceTraps_9.dat','JV_InterfaceTraps_10.dat','JV_InterfaceTraps_11.dat', ]
        Var_files = ['Var_InterfaceTraps_ref.dat','Var_InterfaceTraps_1.dat','Var_InterfaceTraps_2.dat','Var_InterfaceTraps_3.dat','Var_InterfaceTraps_4.dat','Var_InterfaceTraps_5.dat','Var_InterfaceTraps_6.dat','Var_InterfaceTraps_7.dat','Var_InterfaceTraps_8.dat','Var_InterfaceTraps_9.dat','Var_InterfaceTraps_10.dat','Var_InterfaceTraps_11.dat', ]
        JVscaps_lst = [path_scaps / 'iftestref.iv', path_scaps / 'iftest1.iv', path_scaps / 'iftest2.iv',
                     path_scaps / 'iftest3.iv',path_scaps / 'iftest4.iv', path_scaps / 'iftest5.iv',
                     path_scaps / 'iftest6.iv',path_scaps / 'iftest7.iv', path_scaps / 'iftest8.iv',
                     path_scaps / 'iftest9.iv', path_scaps / 'iftest10.iv', path_scaps / 'iftest11.iv']
        NPexp_lst = [path_scaps / 'npiftestref.eb', path_scaps / 'npiftest1.eb', path_scaps / 'npiftest2.eb',
                     path_scaps / 'npiftest3.eb',path_scaps / 'npiftest4.eb', path_scaps / 'npiftest5.eb',
                     path_scaps / 'npiftest6.eb',path_scaps / 'npiftest7.eb', path_scaps / 'npiftest8.eb',
                     path_scaps / 'npiftest9.eb', path_scaps / 'npiftest10.eb', path_scaps / 'npiftest11.eb']
        JV_plot_filename = 'JV_Interface_traps_configuration_SIMsalabim'+ext_save_pic 
        np_plot_filename = 'np_Interface_traps_configuration_SIMsalabim'+ext_save_pic 
        has_TL = True
        SIMsalabim_vs_scaps(path_scaps, str_lst, has_TL, labels, JV_files, Var_files,path2SIMsalabim,Store_Folder, JVscaps_lst, NPexp_lst, JV_plot_filename, np_plot_filename)
