import numpy as np
from scipy.integrate import odeint

x = f_params()
y0 = initial_values()

tspan = range(1801)
t = np.array(tspan)
condition = 3

sim_RP        = np.empty((len(tspan),condition))
sim_ShP       = np.empty((len(tspan),condition))
sim_PI3K_act  = np.empty((len(tspan),condition))
sim_Raf_act   = np.empty((len(tspan),condition))
sim_MEKPP     = np.empty((len(tspan),condition))
sim_ERKPP     = np.empty((len(tspan),condition))
sim_Akt_PI_PP = np.empty((len(tspan),condition))

for i in range(condition):
    if i==0:
        y0[HRG] = 330. # 10nM
    elif i==1:
        y0[HRG] = 33.  # 1nM
    elif i==2:
        y0[HRG] = 3.3  # 0.1nM

    Y = odeint(diffeq,y0,tspan,args=tuple(x))

    sim_RP[:,i] = 2*(Y[:,RP]+Y[:,R_PI3K]+Y[:,R_PI3K_act]+Y[:,R_ShGS]+Y[:,R_ShP]+Y[:,R_Shc])/y0[R]*100.
    sim_ShP[:,i] = Y[:,ShP]/y0[Shc]*100.
    sim_PI3K_act[:,i] = Y[:,PI3K_act]/y0[PI3K]*100.
    sim_Raf_act[:,i] = Y[:,Raf_act]/y0[Raf]*100.
    sim_MEKPP[:,i] = Y[:,MEKPP]/y0[MEK]*100.
    sim_ERKPP[:,i] = Y[:,ERKPP]/y0[ERK]*100.
    sim_Akt_PI_PP[:,i] = Y[:,Akt_PI_PP]/y0[Akt]*100.