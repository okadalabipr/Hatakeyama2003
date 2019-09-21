from .name2idx import variables as V

def initial_values():
    y0 = [0]*V.len_f_vars

    y0[V.R] = 80.
    y0[V.Shc] = 1000.
    y0[V.GS] = 10.
    y0[V.RasGDP] = 120.
    y0[V.Raf] = 100.
    y0[V.MEK] = 120.
    y0[V.ERK] = 1000.
    y0[V.PI3K] = 10.
    y0[V.PI] = 800.
    y0[V.Akt] = 10.
    y0[V.E] = 7.
    y0[V.MKP3] = 2.4
    y0[V.PP2A] = 11.4

    return y0