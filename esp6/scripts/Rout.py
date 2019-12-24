import math
import numpy as np
from matplotlib import pyplot as plt
from labbclass import Analisi
import sys 
import os

file = '../Misure/zener.csv'

Rout = Analisi()
Rout.Rl = Rout.leggi_colonna(file, 0)
Rout.Rl = np.delete(Rout.Rl, -1)
Rout.sigmaRl = Rout.Rl *1/100
Rout.Vout = Rout.leggi_colonna(file, 1)
Rout.sigmaVout = Rout.leggi_colonna(file, 2) *24/200
Rout.Vca = float(Rout.Vout[-1])
Rout.sigmaVca = float(Rout.sigmaVout[-1])
Rout.Vout = np.delete(Rout.Vout, -1)
Rout.sigmaVout = np.delete(Rout.sigmaVout, -1)


Rout.Req = Rout.Rl * (Rout.Vca/Rout.Vout - 1)
Rout.sigmaReq = np.sqrt((Rout.Rl/Rout.Vout * Rout.sigmaVca)**2 + (Rout.Rl*Rout.Vca/(Rout.Vout)**2 * Rout.sigmaVout)**2 + (Rout.sigmaRl*(Rout.Vca/Rout.Vout -1))**2)
Rout.meanReq = np.average(Rout.Req, weights=Rout.sigmaReq)

for Req, sigmaReq in zip(Rout.Req, Rout.sigmaReq):
    print("Rout = {0:0.4f} +- {1:0.4f}".format(Req, sigmaReq))
print(Rout.meanReq)