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
tau = np.ndarray((5,5), dtype=LinearFit)     # matrice una colonna hai i 5 dati fatti con la stessa resistenza; colonne diverse per resistenze diverse
partial_path = 'Misure/Ctot/sistemati/r'
for i in range(5):
    suffisso_1 = str(i+1) + '_'
    for j in range(5):
        suffisso_2 = str(j+1) + '.csv'
        suffisso_2_txt = str(j+1) + '.txt'
        total_path = partial_path + suffisso_1 + suffisso_2
        total_path_txt = partial_path + suffisso_1 + suffisso_2_txt
        tau[i,j] = LinearFit()
        tau[i,j].leggiDati(total_path)
        # le righe che seguono servono a beccare dal txt la risoluzione dell'oscilloscopio
        with open(total_path_txt, 'r') as filetxt:
            lines = filetxt.readlines()
            numero = lines[1][11:14]
            ris = numero[0] + numero[1] + numero[2]
            ris = int(ris)*1e-3
        tau[i,j].add_sigmas(sigmay=ris)


# plt.plot(tau[0,0].xdata, tau[0,0].ydata, '.')
# plt.show()

# REGRESSIONE PER TROVARE TAU SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        tau[i,j].reg_lin(cambiaVariabili=True, y=np.log(tau[i,j].ydata), x=tau[i,j].xdata, sigma=tau[i,j].sigmay/tau[i,j].ydata)

print(1/tau[0,0].B)
'''n_ripetute = 5 # numero di set per ogni misura
tau_tot = np.asarray([])

for tau in tau_tot:
    set_misure = [LinearFit() for i in range(n_ripetute)] # Array di set: per ogni resistenza (per ogni tau) ho 5 set
    tau = set_misure # tau è un array di oggetti della LinearFit: sarà il risultato dell'analisi di questi set'''