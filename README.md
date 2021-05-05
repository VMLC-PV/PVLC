# PVLC codes 

## Author
Vincent M. Le Corre (Main)\
Lena Lindenmeier

## Description
All these codes are written to be used with the open-source drift-diffusion suit written by Prof. L. Jan Anton Koster from University of Groningen. (See [GitHub repository](https://github.com/kostergroup)) They can be used to run the simulation, plot and do analysis of the output from the simulations for both steady-state with [SIMsalabim](https://github.com/kostergroup/SIMsalabim) and transient with ZimT (not open-source yet).\
Codes can be ran on Windows and Linux. However, the running simulations in parallel is not possible on Windows yet. 

## Steady-state simulation codes
- "run_plot_SIMsalabim.py" can be used to run SIMsalabim and plot the output from the JV_file or from the Var_file.
- "SIMsalabim_testing.py" can be used to test new version of SIMsalabim and compare it with [SCAPS](http://scaps.elis.ugent.be/) output for different case scenarios defined in *SIMsalabim_vs_SCAPS.pptx*. (Note: In future versions, SIMsalabim will also be compared to simple physical model which can be solved analytically)

## Transient simulation codes
- "run_plot_zimt.py" can be used to run SIMsalabim and plot the output from the JV_file or from the Var_file.
- "TDCF_zimt.py", "BACE_zimt.py", "TPV_zimt.py", "TPC_zimt.py", "Impedance_zimt.py" and "JV_Hyst_zimt.py" can be used to run zimt and analyze/plot the results from the simulation of the corresponding experiment.
- "tVG_gen.py" contains all the functions to generate the input tVG_file for zimt for the wanted experiment.

## Functions package
- "VLC_useful_func.py" contains all the necessary functions used by the different scripts.
- "plot_settings_screen.py" control the default setting for the plotting functions (Font size, line thickness...)
