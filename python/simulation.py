import numpy as np
from scipy.integrate import odeint

from model.name2idx import f_parameter as C
from model.name2idx import f_variable as V
from model.param_const import f_params
from model.initial_condition import initial_values
from model.differential_equation import diffeq

class Simulation(object):

    x = f_params()
    y0 = initial_values()

    tspan = range(1801)
    t = np.array(tspan)
    condition = 3

    RP        = np.empty((len(tspan),condition))
    ShP       = np.empty((len(tspan),condition))
    PI3K_act  = np.empty((len(tspan),condition))
    Raf_act   = np.empty((len(tspan),condition))
    MEKPP     = np.empty((len(tspan),condition))
    ERKPP     = np.empty((len(tspan),condition))
    Akt_PI_PP = np.empty((len(tspan),condition))

    for i in range(condition):
        if i==0:
            y0[V.HRG] = 330. # 10nM
        elif i==1:
            y0[V.HRG] = 33.  # 1nM
        elif i==2:
            y0[V.HRG] = 3.3  # 0.1nM

        Y = odeint(diffeq,y0,tspan,args=tuple(x))

        RP[:,i] = 2*(Y[:,V.RP]+Y[:,V.R_PI3K]+Y[:,V.R_PI3K_act]+Y[:,V.R_ShGS]+Y[:,V.R_ShP]+Y[:,V.R_Shc])/y0[V.R]*100.
        ShP[:,i] = Y[:,V.ShP]/y0[V.Shc]*100.
        PI3K_act[:,i] = Y[:,V.PI3K_act]/y0[V.PI3K]*100.
        Raf_act[:,i] = Y[:,V.Raf_act]/y0[V.Raf]*100.
        MEKPP[:,i] = Y[:,V.MEKPP]/y0[V.MEK]*100.
        ERKPP[:,i] = Y[:,V.ERKPP]/y0[V.ERK]*100.
        Akt_PI_PP[:,i] = Y[:,V.Akt_PI_PP]/y0[V.Akt]*100.