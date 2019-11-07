from labbclass import LinearFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# PROVA LETTURA FILE
tau_prova = LinearFit()
tau_prova.leggiDati('Misure/Ctot/sistemati/r1_1.csv')

plt.plot(tau_prova.xdata, tau_prova.ydata, '.')
plt.show()

'''n_ripetute = 5 # numero di set per ogni misura
tau_tot = np.asarray([])

for tau in tau_tot:
    set_misure = [LinearFit() for i in range(n_ripetute)] # Array di set: per ogni resistenza (per ogni tau) ho 5 set
    tau = set_misure # tau è un array di oggetti della LinearFit: sarà il risultato dell'analisi di questi set'''