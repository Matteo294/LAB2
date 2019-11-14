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
            dV_line = 1
            dV_pos = [i for i in range(11,14)]
            dV_mul = lines[dV_line][14]
            dV = []
            for idx in dV_pos:
                dV += lines[dV_line][idx]
            dV = float(''.join(dV))
            if dV_mul=='m':
                dV *= 1e-3

            dT_line = 9
            dT_pos = [i for i in range(36,41)]
            dT = []
            for idx in dT_pos:
                dT += lines[dT_line][idx]
            dT = float(''.join(dT))
            dT_mul = lines[9][41]
            if dT_mul=='u':
                dT *= 1e-6
            if dT_mul == 'm':
                dT *= 1e-3

        scariche[i,j].add_sigmas(sigmay=dV, sigmax=dT)


# REGRESSIONE PER TROVARE tau SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        scariche[i,j].reg_lin(trasferisci=True, logy=True)


tau_R = LinearFit()       # per fare la regressione tra 1/tau e 1/R
R = np.array([1e3,99.57e3,21.73e3,39.36e3,9.94e3]) # inserire valori resistenze usate davvero
sigmaR = 1      #temporanea
for i in range(5):      # per ogni colonna i, calcolo la media delle B (B=1/tau)
    B_colonna = np.empty(5)   # vettore delle B per ciascun elemento della colonna
    for j in range(5):         # sulla colonna, calcolo la media delle B
            B_colonna[j] = -scariche[j, i].B
    tau_R.ydata = np.append(tau_R.ydata, np.mean(B_colonna))       #ricorda -B = 1/tau
    tau_R.xdata = np.append(tau_R.xdata, 1/R[i])
    tau_R.sigmay = np.append(tau_R.sigmay, np.std(B_colonna)/(5))       
    tau_R.sigma_reg = tau_R.sigmay
    tau_R.sigmax = np.append(tau_R.sigmax, 1/R[i]*sigmaR)        


tau_R.reg_lin()

print(1/tau_R.B)

# for tau in scariche:
#     for r in tau:
#         y = 1 / r.sigmay * np.log(r.ydata)
#         A = np.matrix.transpose(1/r.sigmay * np.asarray([np.ones(r.xdata.size), r.xdata, 1/r.ydata]))
#         B = np.dot(np.matrix.transpose(A), A)
#         r.offset_params = np.dot( np.dot(np.linalg.inv(B), np.matrix.transpose(A)), y )
#         print(r.offset_params.size)
#         print("Parametri dell'analisi offset (A, B, C):")
#         print(r.offset_params)
 