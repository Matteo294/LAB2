import numpy as np
import math
import uncertainties
from uncertainties import ufloat
from uncertainties.umath import *

def incertezza_H(A_in, B_in, A_out, B_out, dA_in, dB_in, dA_out, dB_out):
    A_in = ufloat(A_in, dA_in)
    B_in = ufloat(B_in, dB_in)
    A_out = ufloat(A_out, dA_out)
    B_out = ufloat(B_out, dB_out)

    C_in_re = A_in
    C_in_im = -B_in
    C_out_re = A_out
    C_out_im = B_out

    H_re = (C_out_re*C_in_re + C_in_im*C_out_re)/(C_in_re**2 + C_in_im**2)
    H_im = (C_out_im*C_in_re - C_out_re*C_in_im)/(C_in_re**2 + C_in_im**2)
    H_abs = H_re**2 + H_im **2
    H_arg = atan(H_im/H_re)

    return {"abs":H_abs.s, "arg":H_arg.s}
