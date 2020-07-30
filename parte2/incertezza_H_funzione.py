import numpy as np
from math import *
import uncertainties
from libphysics import numpify
from uncertainties import ufloat
from uncertainties.umath import *

def incertezza_H(Cin, Cout, t_schermo, freqs):
    Cin_abs = abs(Cin)
    Cout_abs = abs(Cout)
    Cin_fase = float(np.angle(Cin))
    Cout_fase = float(np.angle(Cin))
    
    dVout = 2*Cout_abs*1.5/100
    dVin = 2*Cin_abs*1.5/100
    
    dt = 8e-3*t_schermo

    dCout_abs = dVout
    dCout_fase = numpify(2*pi*freqs*dt)*180/pi
    dCin_abs = dVin
    dCin_fase = numpify(2*pi*freqs*dt)*180/pi

    Cin_abs = ufloat(Cin_abs, dCin_abs)
    Cin_fase = ufloat(Cin_fase, dCin_fase)
    Cout_abs = ufloat(Cout_abs, dCout_abs)
    Cout_fase = ufloat(Cout_fase, dCout_fase)
    
    H_abs = Cout_abs/Cin_abs
    H_fase = Cout_fase-Cin_fase


    return {"abs":H_abs.s, "arg":H_fase.s}
