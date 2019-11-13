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
scariche = np.ndarray((5,5), dtype=LinearFit)     # matrice: in una colonna hai i 5 scariche fatte con la stessa resistenza; colonne diverse per resistenze diverse
partial_path = 'Misure/Ctot/sistemati/r'
for i in range(5):
    suffisso_1 = str(i+1) + '_'
    min = 0
    epsilon = 1e-10
    for j in range(5):
        suffisso_2 = str(j+1) + '.csv'
        suffisso_2_txt = str(j+1) + '.txt'
        total_path = partial_path + suffisso_1 + suffisso_2
        total_path_txt = partial_path + suffisso_1 + suffisso_2_txt
        scariche[i,j] = LinearFit()
        scariche[i,j].leggiDati(total_path)
    #traslo in su se trovo valori negativi
        for value in scariche[i,j].ydata:
            if value < min:
                min = value
        scariche[i,j].ydata += -min + epsilon
    #le righe che seguono servono a beccare dal txt la risoluzione dell'oscilloscopio
        with open(total_path_txt, 'r') as filetxt:
            lines = filetxt.readlines()
            numero = lines[1][11:14]
            ris = numero[0] + numero[1] + numero[2]
            ris = int(ris)*1e-3
        scariche[i,j].add_sigmas(sigmay=ris)


# REGRESSIONE PER TROVARE scariche SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        scariche[i,j].reg_lin(cambiaVariabili=True, y=np.log(scariche[i,j].ydata), x=scariche[i,j].xdata, sigma=scariche[i,j].sigmay/scariche[i,j].ydata, trasferisci=True)


tau_R = LinearFit()       # per fare la regressione tra 1/tau e 1/R
R = np.array([1,1,1,1,1]) # inserire valori resistenze usate davvero
sigmaR = 1
for i in range(5):      # per ogni colonna i, calcolo la media delle B (B=1/tau)
    B_colonna = np.empty(5)   # vettore delle B per ciascun elemento della colonna
    for j in range(5):         # sulla colonna, calcolo la media delle B
            B_colonna[j] = -scariche[j, i].B
    sigma_B = np.std(B_colonna)/(5)
    tau_R.ydata = np.append(tau_R.ydata, np.mean(B_colonna))
    tau_R.xdata = np.append(tau_R.xdata, 1/R[i])
    tau_R.sigmay = np.append(tau_R.sigmay, np.std(B_colonna)/(5))
    print(np.mean(B_colonna))
    tau_R.sigmax = np.append(tau_R.sigmax, 1/R[i]*sigmaR)        # non sono sicura del /n, o /(n-1), controllare

print(tau_R.sigmay)
tau_R.reg_lin(trasferisci=True)
print(1/tau_R.B)

# print(1/scariche[0,0].B, scariche[0,0].sigma_B)
'''n_ripetute = 5 # numero di set per ogni misura
scariche_tot = np.asarray([])

for scariche in scariche_tot:
    set_misure = [LinearFit() for i in range(n_ripetute)] # Array di set: per ogni resistenza (per ogni scariche) ho 5 set
    scariche = set_misure # scariche è un array di oggetti della LinearFit: sarà il risultato dell'analisi di questi set'''