from labbclass import LinearFit
from numpy.linalg import inv
from labbclass import PolynomialFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

# per non far apparire tremila scritte e grafici
enable_plots = True # Scegliere se visualizzare i grafici
ndati = 5000 # Utile per selezionare una porzione di dati

# parametri noti
R_gen = 50
R_dmm = np.array([99.5, 46.75, 99.5+46.75, 99.5*46.75/(99.5+46.75) + 99.13, 99.13+99.5])  
print(R_dmm)
sigmaR = np.array([0.01, 0.009, 0.01, 0.01, 0.01]) 
n_scariche = 5
n_resistenze = 5
L_teorico = 2e-3 #!!!! mettere

passo = [150, 200, 200, 150, 200] # utile per la derivata numerica (valori scelti osservando i grafici della derivata)

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
        scariche[j,i] = LinearFit()
        scariche[j,i].leggiDati(total_path)
        scariche[j,i].R_teorica = R_dmm[j]
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
        scariche[j,i].add_sigmas(sigmay=dV, sigmax=dT)

# Solo valori positivi
_scariche = np.ndarray((5,5), dtype=LinearFit)
for tau, i in zip(scariche, range(scariche.size)):
    _scariche[i] = [LinearFit() for R in tau]
for tau, i in zip(scariche, range(scariche.size)):
    for r, j in zip(tau, range(tau.size)):
        indici = np.where(r.ydata <= 0)
        _scariche[i,j].xdata = np.delete(r.xdata, indici)
        _scariche[i,j].sigmax = np.delete(r.sigmax, indici)
        _scariche[i,j].ydata = np.delete(r.ydata, indici)
        _scariche[i,j].sigmay = np.delete(r.sigmay, indici)
        _scariche[i,j].R_teorica = r.R_teorica

# Ritaglio solo i primi dati
dati = np.ndarray((5,5), dtype=LinearFit)
for tau, i in zip(_scariche, range(_scariche.size)):
    dati[i] = [LinearFit() for R in tau]
for tau, i in zip(_scariche, range(_scariche.size)):
    for r, j in zip(tau, range(tau.size)):
        indici = np.arange(ndati, r.xdata.size, 1)
        dati[i,j].xdata = np.delete(r.xdata, indici)
        dati[i,j].sigmax = np.delete(r.sigmax, indici)
        dati[i,j].ydata = np.delete(r.ydata, indici)
        dati[i,j].sigmay = np.delete(r.sigmay, indici)
        dati[i,j].R_teorica = r.R_teorica

tau_R = LinearFit() # Qui tengo la media delle tau (in realtà -1/tau) per ogni resisteza. Poi faro' una regressione tra i vari tau e le corrispettive resistenze

# Calcolo derivata numerica dV/dt
sigma_y = np.array([])
for t, i in zip(dati, range(dati.size)):
    tau = [LinearFit() for i in range(t.size)]
    _da_mediare = np.array([])
    for r, j in zip(t, range(t.size)):
        valori, sigma = r.derivata_numerica(r.xdata, r.ydata, sigmax=r.sigmax, sigmay=r.sigmay, passo=passo[i])
        tau[j].ydata = np.append(tau[j].ydata, valori)
        tau[j].xdata = np.append(tau[j].xdata, r.ydata[:-passo[i]])
        tau[j].add_sigmas(sigmax=r.sigmay[:-passo[i]], sigmay=sigma)
        tau[j].reg_lin()
        tau[j].chi_quadro()  
        _da_mediare = np.append(_da_mediare, tau[j].B)

    tau_R.ydata = np.append(tau_R.ydata, -np.mean(_da_mediare))
    sigma_y = np.append(sigma_y, np.std(_da_mediare) / np.sqrt(5))
    tau_R.xdata = np.append(tau_R.xdata, R_dmm[i])
    del tau 
    del _da_mediare

sigma_x = np.copy(sigmaR)
tau_R.add_sigmas(sigmax=sigma_x, sigmay=sigma_y)
tau_R.reg_lin()
tau_R.chi_quadro()

if enable_plots:

    print(tau_R.xdata)

    ax1 = plt.subplot(2,1,1)
    tau_R.plotData(xlabel=r"Resistenza [$\Omega$]", ylabel=r"$\frac{1}{\tau}$  $[\frac{\Omega}{H}]$")
    tau_R.data_plot.set_markersize(6)

    plt.subplot(2,1,2, sharex=ax1)
    tau_R.residui(xlabel=r"Resistenza [$\Omega$]", ylabel="Residui")
    tau_R.residui_plot.set_markersize(6)

    plt.xlim(35, 210)
    plt.show()

print(f"L = {-1/tau_R.B} \u00B1 {tau_R.sigmaB/tau_R.B**2}")
print(f"Rgen + RL = {tau_R.A/tau_R.B} \u00B1 {math.sqrt( (tau_R.sigmaA/tau_R.B)**2 + (tau_R.A*tau_R.sigmaB/tau_R.B**2)**2 )}")
print(tau_R)