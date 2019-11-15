from labbclass import LinearFit
from labbclass import PolynomialFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

enable_offset_printing = True
enable_tau_printing = False
enable_plots = True
ndati = 2000 # !!!!!!!!!!!! VALE PER TUTTO IL PROGRAMMA: ovunque ci sia _scariche invece di scariche !!!!!!!!!!!!!!

C_teo = 35.5e-9
R_teoriche = np.asarray([1.001, 99.570, 21.73, 39.36, 9.94])

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# LETTURA FILE
scariche = np.ndarray((5,5), dtype=LinearFit)     # matrice: in una colonna hai i 5 scariche fatte con la stessa resistenza; colonne diverse per resistenze diverse
partial_path = 'Misure/Ctot/sistemati/r'
for j in range(5):
    suffisso_1 = str(j+1) + '_'
    min = 0
    epsilon = 1e-10
    for i in range(5):
        suffisso_2 = str(i+1) + '.csv'
        suffisso_2_txt = str(i+1) + '.txt'
        total_path = partial_path + suffisso_1 + suffisso_2
        total_path_txt = partial_path + suffisso_1 + suffisso_2_txt
        scariche[i,j] = LinearFit()
        scariche[i,j].leggiDati(total_path)
        scariche[i,j].R_teorica = R_teoriche[j]
    #traslo in su se trovo valori negativi
        '''for value in scariche[i,j].ydata:
            if value < min:
                min = value
        scariche[i,j].ydata += -min + epsilon'''
    #le righe che seguono servono a beccare dal txt la risoluzione dell'oscilloscopio
        with open(total_path_txt, 'r') as filetxt:
            lines = filetxt.readlines()
            # dV
            dV_line = 1
            dV_pos = [i for i in range(11,14)]
            dV_mul = lines[dV_line][14]
            dV = []
            for idx in dV_pos:
                dV += lines[dV_line][idx]
            dV = float(''.join(dV))
            if dV_mul=='m':
                dV *= 1e-3
            # dT
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


# Ritaglio i primi 200 elementi dell'array e li metto un array provvisorio in _scariche
_scariche = np.ndarray((5,5), dtype=LinearFit)
for tau, i in zip(scariche, range(scariche.size)):
    _scariche[i] = [LinearFit() for R in tau]

for tau, i in zip(scariche, range(scariche.size)):
    for r, j in zip(tau, range(tau.size)):
        _scariche[i,j].xdata = np.delete(r.xdata, np.arange(ndati, r.xdata.size, 1))
        _scariche[i,j].sigmax = np.delete(r.sigmax, np.arange(ndati, r.xdata.size, 1))
        _scariche[i,j].ydata = np.delete(r.ydata, np.arange(ndati, r.xdata.size, 1))
        _scariche[i,j].sigmay = np.delete(r.sigmay, np.arange(ndati, r.xdata.size, 1))
        _scariche[i,j].R_teorica = r.R_teorica

# I coefficienti che minimizzano il chi quadro (massimizzano la Likelihood) sono
# x = inv(A' A) A' y
for tau, i in zip(_scariche, range(scariche.size)):
    for r, j in zip(tau, range(scariche[i].size)):
        # y = A + Bt + Cx   con x = 1/Vm  e  y = log(Vm)    
        x = 1/r.ydata
        t = r.xdata
        y = np.log(r.ydata)
        M = np.transpose(np.asarray([np.ones(r.xdata.size), t, x]) / r.sigmay)
        offset = PolynomialFit(M, y)
        offset.minimizza()
        if enable_offset_printing:
            print(offset)  
        _scariche[i,j].ydata += offset.coefficienti[2]  # Traslo i dati

#scariche[2,4] = scariche[3,4]


# REGRESSIONE PER TROVARE tau SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        _scariche[i,j].reg_lin(trasferisci=True, logy=True)
        _scariche[i,j].chi_quadro(logy=True)
#print(scariche[0,0].chi_ridotto)

for i in range(5):
    for j in range(5):
        if enable_tau_printing:
            print(1/_scariche[i,j].B)

tau_R = LinearFit()       # per fare la regressione tra 1/tau e 1/R
R = np.array([1e3, 99.57e3, 21.73e3, 39.36e3, 9.94e3]) # inserire valori resistenze usate davvero
sigmaR = 1      #temporanea
for j in range(5):      # per ogni colonna j, calcolo la media delle B (B=1/tau)
    B_colonna = np.empty(5)   # vettore delle B per ciascun elemento della colonna
    for i in range(5):         # sulla colonna, calcolo la media delle B
            B_colonna[i] = -_scariche[i, j].B
    tau_R.ydata = np.append(tau_R.ydata, np.mean(B_colonna))       #ricorda -B = 1/tau
    tau_R.xdata = np.append(tau_R.xdata, 1/R[j])
    tau_R.sigmay = np.append(tau_R.sigmay, np.std(B_colonna)/(5))       
    tau_R.sigma_reg = tau_R.sigmay
    tau_R.sigmax = np.append(tau_R.sigmax, 1/R[j]*sigmaR)        
    if enable_tau_printing:
        print(1/tau_R.ydata[j])

tau_R.reg_lin()
tau_R.chi_quadro()
if enable_tau_printing:
    print(tau_R)

if enable_plots:
    tau_R.plotData()
    plt.show()
    tau_R.residui()
    plt.show()
