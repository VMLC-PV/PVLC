################################################################
############### Function to format tVG files ###################
################################################################
# by Vincent M. Le Corre
# Package import
from unicodedata import name
import pandas as pds
import matplotlib.pyplot as plt
import numpy as np 
import math
from pathlib import Path
from scipy import integrate

def gaussian_pulse(t, tpulse, width, Imax):
    """Returns a gaussian pulse

    Parameters
    ----------
    t : 1-D sequence of floats
        t time axis (unit: s)
    tpulse : float
        tpulse center of the pulse (unit: s)
    width : float
        width of the pulse (unit: s)
    Imax : float
        Imax maximum of the pulse

    Returns
    -------
    1-D sequence of floats
        Vector containing the gaussian pulse
    """    
    return Imax *np.exp(-np.power(t - tpulse, 2.) / (2 * np.power(width, 2.)))


def zimt_light_decay(tmin,tmax,Gstart,Gfinal,Va,tstep,trf = 20e-9,time_exp =False,tVG_name='tVG.txt'):
    """Make tVG file for light decay experiment

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Gstart : float
        initial generation rate, i.e. light intensity (steady-state) (unit: m^-3 s^-1)
    Gfinal : float
        final generation rate (unit: m^-3 s^-1)
    Va : float
        applied voltage (unit: V)
    tstep : float
        time step if lin or number of point if time_exp = True (unit: s)
    trf : float, optional
        LED/laser fall/rise time, by default 20e-9 (unit: s)
    time_exp : bool, optional
        If True exponential time step is used, else linear time step, by default False
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    

    # Chose between exponential or linear time step
    if time_exp == True:
        t = np.geomspace(tmin,tmax,tstep)
    else :
        t = np.arange(tmin,tmax,tstep)

    t=np.insert(t,0,0)
    V,G = [],[]
    # Calculate the light decay assuming and exponential decay of the generation rate
    # with a lifetime (trf) of the LED/laser fall/rise time
    for i in t:
        G.append((Gstart-Gfinal)*np.exp(-i/trf)+Gfinal)
        V.append(Va)

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')

def zimt_voltage_step(tmin,tmax,Vstart,Vfinal,Gen,steps=100,trf = 10e-9,time_exp =False,tVG_name='tVG.txt'):
    """Make tVG file for Voltage decay experiment

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Vstart : float
        initial applied voltage (steady-state) (unit: V)
    Vfinal : float
        final applied voltage (unit: V)
    Gen : float
        constant generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    steps : int, optional
        number of time steps, by default 100
    trf : float, optional
        Voltage pulse fall/rise time , by default 10e-9 (unit: s)
    time_exp : bool, optional
        If True exponential time step is used, else linear time step, by default False
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """ 

    # Chose between exponential or linear time step
    if time_exp == True:
        t = np.geomspace(tmin,tmax,num=steps)
    else :
        t = np.linspace(0,tmax,int(steps),endpoint=True)

    t=np.insert(t,0,0)
    V,G = [],[]
    # Calculate the voltage decay assuming and exponential decay of the voltage
    # with a lifetime (trf) of the voltage pulse fall/rise time
    for i in t:
        V.append((Vstart-Vfinal)*np.exp(-i/trf)+Vfinal)
        G.append(Gen)

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')


def zimt_JV_sweep(Vstart,Vfinal,scan_speed,Gen,steps,time_exp =False,tVG_name='tVG.txt'):
    """Make tVG file for one JV sweep experiment

    Parameters
    ----------
    Vstart : float
        initial applied voltage (steady-state) (unit: V)
    Vfinal : float
        final applied voltage (unit: V)
    scan_speed : float
        scan speed (unit: V/s)
    Gen : float
        constant generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    steps : int
        number of JV points
    time_exp : bool, optional
        If True exponential time step is used, else linear time step, by default False
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    
    # Calculate duration of the JV sweep depending on the scan speed
    tmax = abs(Vfinal - Vstart)/scan_speed

    # Chose between exponential or linear time step
    if time_exp == True:
        t = np.geomspace(0,tmax,int(steps))
        t=np.insert(t,0,0)
    else :
        t = np.linspace(0,tmax,int(steps),endpoint=True)
    
    V,G = [],[]
    for i in t:
        V.append(scan_speed*i + Vstart)
        G.append(Gen)

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')

def zimt_JV_double_sweep(Vstart,Vfinal,scan_speed,Gen,steps,time_exp =False,tVG_name='tVG.txt'):
    """Make tVG file for double JV sweep experiment
    Scan voltage back and forth

    Parameters
    ----------
    Vstart : float
        initial applied voltage (steady-state) (unit: V)
    Vfinal : float
        final applied voltage (unit: V)
    scan_speed : float
        scan speed (unit: V/s)
    Gen : float
        constant generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    steps : int
        number of JV points 
    time_exp : bool, optional
        If True exponential time step is used, else linear time step, by default False
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """
    # Calculate duration of one sweep, the total experiment duration will be 2*tmax  
    tmax = abs((Vfinal - Vstart)/scan_speed)
    
    # Chose between exponential or linear time step
    if time_exp == True:
        t = np.geomspace(0,2*tmax,int(steps),endpoint=True)
        t=np.insert(t,0,0)
    else :
        t1 = np.linspace(0,tmax,int(steps/2),endpoint=True)
        t2 = np.linspace(tmax,2*tmax,int(steps/2),endpoint=True)
        t2 = np.delete(t2,[0])
        t = np.append(t1,t2)
    
    V,G = [],[]
    for i in t:
        if i <= tmax:
            V.append(np.sign(Vfinal-Vstart)*scan_speed*i + Vstart)
        else:
            V.append(-np.sign(Vfinal-Vstart)*(scan_speed*i) + np.sign(Vfinal-Vstart)*Vstart +2*Vfinal)
        G.append(Gen)

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')

def zimt_tdcf(tmin,tmax,Vpre,Vcol,Gen,tpulse,tstep,tdelay,width_pulse = 2e-9,tVp = 10e-9,time_exp=False,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for time-delayed collection field (TDCF)

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Vpre : float
        initial applied voltage (steady-state) or pre-bias in typical TDCF language (unit: V)
    Vcol : float
        final applied voltage or collection bias in typical TDCF language (unit: V)
    Gen : float
        Total number of carrier generated by the gaussian pulse (unit: m^-3)
    tpulse : float
        middle of the gaussian pulse (unit: s)
    tstep : float
        time step for the linear regime (unit: s)
    tdelay : float
        delay between middle of laser pulse and voltage switch (unit: s)
    width_pulse : float, optional
        width of the light pulse (unit: s), by default 2e-9
    tVp : float, optional
        Voltage pulse fall/rise time (unit: s), by default 10e-9
    time_exp : bool, optional
        if True chose exponential time step else keep time step linear, by default False
    steps : int, optional
        if time_exp = True number of exponential time step after voltage switch, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    
    
    tswitch = width_pulse + tdelay
    V,G = [],[]

    # Define laser pulse by a Gaussian and always keep tstep to 1e-10 for during the pulse to ensure that the amount of photon is consistent when changing tstep
    t = np.arange(tmin,width_pulse*1.5,1e-9)
    
    for i in t:
        if i < tswitch:
            V.append(Vpre)
        else:
            V.append((Vpre-Vcol)*np.exp(-i/tVp)+Vcol)
    G = gaussian_pulse(t,1.5*width_pulse/2,width_pulse,1)

    if width_pulse< tswitch:
        # time step before voltage delay 
        t1 = np.arange(width_pulse*1.5,tswitch-tstep,tstep)

        for i in t1:
            if i < tswitch:
                V.append(Vpre)
            else:
                V.append((Vpre-Vcol)*np.exp(-i/tVp)+Vcol)
            G=np.append(G,0)

        # Begin of the voltage voltage delay 
        if time_exp == True:
            t2 = np.geomspace(tswitch,tmax,num=steps)
        else :
            t2 = np.arange(tswitch,tmax,tstep)

        for i in t2:
                V.append((Vpre-Vcol)*np.exp(-i/tVp)+Vcol)
                G=np.append(G,0)
        
        t = np.append(t,t1)
        t = np.append(t,t2)
    else:
        # Begin of the voltage voltage delay 
        if time_exp == True:
            t2 = np.geomspace(1.5*width_pulse,tmax,num=steps)
        else :
            t2 = np.arange(1.5*width_pulse,tmax,tstep)

        for i in t2:
                V.append((Vpre-Vcol)*np.exp(-i/tVp)+Vcol)
                G=np.append(G,0)
        
        t = np.append(t,t2)

    # Insert initial conditions
    t = np.insert(t,0,0)
    V = np.asarray(V)
    V = np.insert(V,0,Vpre)
    G = np.asarray(G)
    G = np.insert(G,0,0)
    if Gen > 0:
        # ensure that the total number of generated charges is equal to Gen
        int_G  = integrate.cumtrapz(G, t, initial=0)
        G = G*Gen/int_G[-1]
    else:
        G = 0*G

    tVG = pds.DataFrame(np.stack([t,V,G]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e') 


def zimt_BACE(tmin,tmax,Gen,Vpre,Vextr,tstep,tLp = 20e-9,tVp = 10e-9,time_exp=False,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for bias-assisted charge extraction (BACE)

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Gen : float
        initial generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    Vpre : float
        initial applied voltage (or pre-bias) (unit: V)
    Vextr : float
        extraction voltage (unit: V)
    tstep : float
        time step for the linear regime (unit: s)
    tLp : float, optional
        LED pulse fall/rise time (unit: s), by default 20e-9
    tVp : float, optional
        Voltage pulse fall/rise time (unit: s), by default 10e-9
    time_exp : bool, optional
        if True chose exponential time step else keep time step linear, by default False
    steps : int, optional
        if time_exp = True number of exponential time step after voltage switch, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    

    V,G = [],[]
    
    if time_exp == True:
        t = np.geomspace(tmin,tmax,num=steps)
    else :
        t = np.arange(tmin,tmax,tstep)
    t=np.insert(t,0,0)
    for i in t:
            V.append((Vpre-Vextr)*np.exp(-i/tVp)+Vextr)
            G=np.append(G,(Gen)*np.exp(-i/tLp))

    # t=np.insert(t,0,0)
    # V=np.insert(t,0,Vpre)
    # G=np.insert(t,0,Gen)

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e') 

def zimt_TPV(tmin,tmax,Gen_pulse,G0,tstep,tpulse,width_pulse = 2e-9,time_exp =False,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for transient photovoltage (TPV) experiment

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Gen_pulse : float
        max generation rate (i.e. light intensity) of the gaussian pulse (unit: m^-3 s^-1)
    G0 : float
        background generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    tstep : float
        time step for the linear regime (unit: s)
    tpulse : float
        middle of the gaussian pulse (unit: s)
    width_pulse : float, optional
        width of the light pulse (unit: s), by default 2e-9
    time_exp : bool, optional
        if True chose exponential time step else keep time step linear, by default False
    steps : int, optional
        if time_exp = True number of exponential time step after voltage switch, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    

    # Define laser pulse by a Gaussian and always keep tstep to 1e-9 for during the pulse to ensure that the amount of photon is consistent when changing tstep
    t = np.arange(tmin,tpulse + width_pulse*3 - 1e-9,1e-9)
    t=np.insert(t,0,0)
    V,G = [],[]

    for i in t:
        V.append('oc')   
    G = gaussian_pulse(t,tpulse,width_pulse,Gen_pulse)


    if time_exp == True:
        t1 = np.geomspace(tpulse + width_pulse*3,tmax,steps)
    else :
        t1 = np.arange(tpulse + width_pulse*3,tmax,tstep)

    for i in t1:
        G=np.append(G,G0)
        V.append('oc')
    
    t = np.append(t,t1)

    for i in range(len(G)): #set G = 0 when G is too small (for stability in ZimT)
        if G[i] < G0:
            G[i] = G0

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')


def zimt_TPC(tmin,tmax,Gen_pulse,G0,tstep,tpulse,width_pulse = 2e-9,time_exp =False,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for transient photocurrent (TPC) experiment

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Gen_pulse : float
        max generation rate (i.e. light intensity) of the gaussian pulse (unit: m^-3 s^-1)
    G0 : float
        background generation rate (i.e. light intensity) (unit: m^-3 s^-1)
    tstep : float
        time step for the linear regime (unit: s)
    tpulse : float
        middle of the gaussian pulse (unit: s)
    width_pulse : float, optional
        width of the light pulse (unit: s), by default 2e-9
    time_exp : bool, optional
        if True chose exponential time step else keep time step linear, by default False
    steps : int, optional
        if time_exp = True number of exponential time step after voltage switch, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """   
    # Define laser pulse by a Gaussian and always keep tstep to 1e-9 for during the pulse to ensure that the amount of photon is consitent when changing tstep
    t = np.arange(tmin,tpulse + width_pulse*3 - 1e-9,1e-9)
    t=np.insert(t,0,0)
    V,G = [],[]

    for i in t:
        V.append(0)   
    G = gaussian_pulse(t,tpulse,width_pulse,Gen_pulse)


    if time_exp == True:
        t1 = np.geomspace(tpulse + width_pulse*3,tmax,steps)
    else :
        t1 = np.arange(tpulse + width_pulse*3,tmax,tstep)

    for i in t1:
        G=np.append(G,G0)
        V.append(0)
    
    t = np.append(t,t1)

    for i in range(len(G)): #set G = G0 
        if G[i] < G0:
            G[i] = G0
        if G[i] < 1:
            G[i] = 0

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e')

def zimt_CELIV(tmin,tmax,Voffset,slopeV,Vpulse,Gen,tpulse,tstep,tdelay,width_pulse = 2e-9,time_exp=False,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for charge extraction by linearly increasing voltage (CELIV) experiment

    Parameters
    ----------
    tmin : float
        first time step after 0 (unit: s)
    tmax : float
        final time step (unit: s)
    Voffset : float
        initial applied voltage (steady-state) (unit: V)
    slopeV : float
        slope of the applied voltage increase (V/s)
    Vpulse : float
        Voltage pulse duration (unit: s)
    Gen : float
        max generation rate (i.e. light intensity) of the gaussian pulse (unit: m^-3 s^-1)
    tpulse : float
        middle of the gaussian pulse (unit: s)
    tstep : float
        time step for the linear regime (unit: s)
    tdelay : float
        delay between middle of laser pulse and voltage switch (unit: s)
    width_pulse : float, optional
        width of the light pulse (unit: s), by default 2e-9
    time_exp : bool, optional
        if True chose exponential time step else keep time step linear, by default False
    steps : int, optional
        if time_exp = True number of exponential time step after voltage switch, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    
    
    tswitch = tpulse + tdelay
    V,G = [],[]

    # Define laser pulse by a Gaussian and always keep tstep to 1e-9 for during the pulse to ensure that the amount of photon is consitent when changing tstep
    t = np.arange(tmin,tpulse + width_pulse*3 - 1e-9,1e-9)
    t=np.insert(t,0,0)

    for i in t:
        if i < tswitch:
            V.append(Voffset)
        else:
            V.append(slopeV*(i-tswitch)+Voffset)
    

    if tpulse + width_pulse*3 < tswitch:
        # time step before voltage delay 
        t1 = np.arange(tpulse + width_pulse*3,tswitch-tstep,tstep)

        for i in t1:
            if i < tswitch:
                V.append(Voffset)
            else:
                V.append(slopeV*(i-tswitch)+Voffset)
            # G=np.append(G,0)

        # Begin of the voltage voltage delay 
        if time_exp == True:
            t2 = np.geomspace(tswitch,tmax,num=steps)
        else :
            t2 = np.arange(tswitch,tmax,tstep)

        for i in t2:
                V.append(slopeV*(i-tswitch)+Voffset)
                # G=np.append(G,0)
        
        t = np.append(t,t1)
        t = np.append(t,t2)
    else:
        # Begin of the voltage voltage delay 
        if time_exp == True:
            t2 = np.geomspace(tpulse + width_pulse*3,tmax,num=steps)
        else :
            t2 = np.arange(tpulse + width_pulse*3,tmax,tstep)

        for i in t2:
                V.append(slopeV*(i-tswitch)+Voffset)
                # G=np.append(G,0)
        
        t = np.append(t,t2)

    G = gaussian_pulse(t,tpulse,width_pulse,Gen)

    for i in range(len(G)): #set G = 0 when G is too small (for stability in ZimT)
        if G[i] < 1:
            G[i] = min(0,G[i])

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e') 


def zimt_impedance(Vapp,Vamp,freq,Gen,steps=100,tVG_name='tVG.txt'):
    """Make tVG file for impedance experiment


    Parameters
    ----------
    Vapp : [type]
        offset applied voltage (steady-state) (unit: V)
    Vamp : [type]
        amplitude of the voltage perturbation (unit: V) 
    freq : [type]
        frequency of the oscillation (unit: Hz)
    Gen : [type]
        max generation rate (i.e. light intensity) of the gaussian pulse (unit: m^-3 s^-1)
    steps : int, optional
        number of time step, by default 100
    tVG_name : str, optional
        tVG_file name, by default 'tVG.txt'
    """    

    w = 2*math.pi*freq
    t = np.linspace(0,3/freq,steps)

    G,V = [],[]
    for i in t:
        G.append(Gen)
        V.append(Vapp + Vamp * np.sin(w*i))
 

    tVG = pds.DataFrame(np.stack([t,np.asarray(V),np.asarray(G)]).T,columns=['t','Vext','Gehp'])

    tVG.to_csv(tVG_name,sep=' ',index=False,float_format='%.3e') 

def zimt_TID(Vfill, tfill, fill_steps, Vdrift, tdrift, Vamp, freq, drift_steps, Gen, tVG_path, tVG_name='tVG.txt'):
    """Make tVG file for tid experiment

    Parameters
    ----------
    Vfill : [float]
        Filling voltage (unit: V)
    tfill : [float]
        Time of filling pulse (unit: s)
    fill_steps : [float]
        Number of steps in the filling steps (unit: none)
    Vdrift : [float]
        Drift voltage (unit: V) 
    tdrift : [float]
        Time of drifting (unit: s)
    Vamp : [float]
        Amplitude of perturbation voltage (unit V)
    freq : [float]
        frequency of the ac perturbation (unit: Hz)
    drift_steps : [float]
        Number of steps in one perturbation period (unit: none)
    Gen : [type]
        max generation rate (i.e. light intensity) of the gaussian pulse (unit: m^-3 s^-1)
    tVG_name : str, optional
    """
    t_fill = np.linspace(0, tfill, fill_steps) 
    V_fill = np.empty(len(t_fill))
    V_fill.fill(Vfill)
    
    t_drift = np.linspace(tfill, tfill+tdrift, int(drift_steps*freq))[1:]
    V_drift = Vdrift + Vamp * np.sin(2*np.pi*freq*(t_drift-tdrift))

    

    t = np.concatenate((t_fill, t_drift))
    V = np.concatenate((V_fill, V_drift))
    G = np.empty(len(t))
    G.fill(Gen)

    tVG = pds.DataFrame(np.transpose([t, V, G]), columns=['t','Vext','Gehp'])

    tVG.to_csv(Path(tVG_path,tVG_name),sep=' ',index=False,float_format='%.3e') 

    
