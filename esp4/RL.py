from labbclass import LinearFit
from numpy.linalg import inv
from labbclass import PolynomialFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

# per non far apparire tremila scritte e grafici
enable_offset_printing = False
enable_tau_printing = True
enable_plots = True
ndati = 5000 # !!!!!!!!!!!! VALE PER TUTTO IL PROGRAMMA: ovunque ci sia _scariche invece di scariche !!!!!!!!!!!!!!

# parametri noti
R_gen = 50
R_dmm = np.asarray([99.5, 46.75, 99.5+46.75, 99.5*46.75/(99.5+46.75) + 99.13, 99.13+99.5])*1000        
sigmaR = np.array([11, 10, 15, 6000, 16])
n_scariche = 5
n_resistenze = 5
L_teorico = 2e-3 #!!!! mettere

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# LETTURA FILE
scariche = np.ndarray((5,5), dtype=LinearFit)     # matrice: in una colonna hai i 5 scariche fatte con la stessa resistenza; colonne diverse per resistenze diverse
partial_path = 'dati/RL/sistemati/r'
for j in range(n_resistenze):
    suffisso_1 = str(j+1) + '_'
    min = 0
    epsilon = 1e-10
    for i in range(n_scariche):
        suffisso_2 = str(i+1) + '.csv'
        suffisso_2_txt = str(i+1) + '.txt'
        total_path = partial_path + suffisso_1 + suffisso_2
        total_path_txt = partial_path + suffisso_1 + suffisso_2_txt
        scariche[i,j] = LinearFit()
        scariche[i,j].leggiDati(total_path)
        scariche[i,j].R_teorica = R_dmm[j]
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
                dV *= 1e-3 / math.sqrt(12)
                #dV = (dV*0.1 + 0.03*8*dV + 2e-3 ) / math.sqrt(12)
            # dT
            dT_line = 9
            dT_pos = [i for i in range(36,41)]
            dT = []
            for idx in dT_pos:
                dT += lines[dT_line][idx]
            dT = float(''.join(dT))
            dT_mul = lines[9][41]
            if dT_mul=='u':
                dT *= 1e-6 / math.sqrt(12)
            if dT_mul == 'm':
                dT *= 1e-3 / math.sqrt(12)
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
    if enable_offset_printing:
        print("------------------------------------------------------------------------------------------------------------------------- \n")
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

# REGRESSIONE PER TROVARE 1/tau SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        _scariche[i,j].reg_lin(trasferisci=True, logy=True)
        _scariche[i,j].chi_quadro(logy=True)

tau_R = LinearFit()                         # per fare la regressione tra 1/tau e 1/R

for j in range(n_resistenze):               # per ogni colonna j, calcolo la media delle B (B=-1/tau)
    B_colonna = np.empty(n_scariche)        # vettore delle B per ciascun elemento della colonna
    for i in range(5):                      # sulla colonna, calcolo la media delle B
            B_colonna[i] = -_scariche[i, j].B
    tau_R.ydata = np.append(tau_R.ydata, np.mean(B_colonna))       
    tau_R.xdata = np.append(tau_R.xdata, R_dmm[j])
    tau_R.sigmay = np.append(tau_R.sigmay, np.std(B_colonna)/(n_scariche))       
    tau_R.sigma_reg = tau_R.sigmay
    tau_R.sigmax = np.append(tau_R.sigmax, sigmaR[j])        

tau_R.reg_lin()
tau_R.chi_quadro()

if enable_tau_printing:
    print("L = ", 1/tau_R.B, " +- ", 1/tau_R.B**2 * tau_R.sigmaB)
    print("L teorico = ", L_teorico)
    print("R_gen + R_L = ", tau_R.A / tau_R.B, " (teorico 50.5 ohm circa)")

if enable_plots:
    tau_R.plotData()
    plt.show()
    tau_R.residui()
    plt.show()
