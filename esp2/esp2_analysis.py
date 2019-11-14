from labbclass import LinearFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

enable_offset_printing = True
C_teo = 35.5e-9
R_teoriche = [1.001, 99.570, 21.73, 39.36, 9.94]

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
        scariche[i,j].R_teorica = R_teoriche[i]
    #traslo in su se trovo valori negativi
        '''for value in scariche[i,j].ydata:
            if value < min:
                min = value
        scariche[i,j].ydata += -min + epsilon'''
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

# Ritaglio i primi 200 elementi dell'array e li metto un array provvisorio in _scariche
_scariche = np.ndarray((5,5), dtype=LinearFit)
for tau, i in zip(scariche, range(scariche.size)):
    _scariche[i] = [LinearFit() for R in tau]

for tau, i in zip(scariche, range(scariche.size)):
    for r, j in zip(tau, range(tau.size)):
        _scariche[i,j].xdata = np.delete(r.xdata, np.arange(2000, r.xdata.size, 1))
        _scariche[i,j].sigmax = np.delete(r.sigmax, np.arange(2000, r.xdata.size, 1))
        _scariche[i,j].ydata = np.delete(r.ydata, np.arange(2000, r.xdata.size, 1))
        _scariche[i,j].sigmay = np.delete(r.sigmay, np.arange(2000, r.xdata.size, 1))
        _scariche[i,j].R_teorica = r.R_teorica

# Risolvo con il metodo di Cramer la soluzione alle equazioni che massimizzano la likelihood (o minimizzano il chi2)
# d(chi2)/da = 0, dove a è uno dei coefficienti da ricavare
for tau in _scariche:
    print('\n------------------------------------------------------------------------------------------\n')
    for r in tau:

        # y = A + Bt + Cx   con x = 1/Vm  e  y = log(Vm)    
        x = 1/r.ydata
        t = r.xdata
        y = np.log(r.ydata)

        w = 1/r.sigmay**2

        M = np.asarray([
            [np.sum(w), np.sum(w*t), np.sum(w*x)],
            [np.sum(w*t), np.sum(w*(t**2)), np.sum(w*x*t)],
            [np.sum(w*x), np.sum(w*x*t), np.sum(w*(x**2))]
        ])

        delta = np.linalg.det(M)

        MA = np.copy(M)
        MB = np.copy(M)
        MC = np.copy(M)
        MA[:,0] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*x*y)]
        MB[:,1] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*y*x)]
        MC[:,2] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*y*x)]

        A = np.linalg.det(MA) / delta
        B = np.linalg.det(MB) / delta
        C = np.linalg.det(MC) / delta
        if enable_offset_printing:
            print("Risultati: \t A = {0:.4f} \t B = {1:.4f} \t C = {2:.4f} \t resistenza reg: {3:.4f} \t resistenza teo: {4:.4f}".format(A, B, C, -1/B/C_teo/1000, r.R_teorica))


 
# print(1/scariche[0,0].B, scariche[0,0].sigma_B)
'''n_ripetute = 5 # numero di set per ogni misura
scariche_tot = np.asarray([])

for scariche in scariche_tot:
    set_misure = [LinearFit() for i in range(n_ripetute)] # Array di set: per ogni resistenza (per ogni scariche) ho 5 set
    scariche = set_misure # scariche è un array di oggetti della LinearFit: sarà il risultato dell'analisi di questi set'''
