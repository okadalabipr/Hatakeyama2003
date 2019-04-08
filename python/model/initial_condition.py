def initial_values():
    y0 = [0]*len(F_V)

    y0[R] = 80.
    y0[Shc] = 1000.
    y0[GS] = 10.
    y0[RasGDP] = 120.
    y0[Raf] = 100.
    y0[MEK] = 120.
    y0[ERK] = 1000.
    y0[PI3K] = 10.
    y0[PI] = 800.
    y0[Akt] = 10.
    y0[E] = 7.
    y0[MKP3] = 2.4
    y0[PP2A] = 11.4

    return y0