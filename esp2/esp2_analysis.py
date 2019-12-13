from labbclass import LinearFit
from labbclass import PolynomialFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

enable_offset_printing = False
enable_tau_printing = True
enable_plots = True
ndati = 3000 # !!!!!!!!!!!! VALE PER TUTTO IL PROGRAMMA: ovunque ci sia dati invece di scariche !!!!!!!!!!!!!!

passo = 30 # utile per la derivata numerica

C_teo = 35.5e-9
R_dmm = np.asarray([1.001e3, 99.570e3, 21.73e3, 39.36e3, 9.94e3])

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
        scariche[i,j].R_teorica = R_dmm[j]
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

# Solo valori positivi
_scariche = np.ndarray((5,5), dtype=LinearFit)
for tau, i in zip(scariche, range(scariche.size)):
    _scariche[i] = [LinearFit() for R in tau]
for tau, i in zip(scariche, range(scariche.size)):
    for r, j in zip(tau, range(tau.size)):
        indici = np.where(r.ydata <= 0.1)
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
for t, i in zip(dati, range(dati.size)):
    tau = [LinearFit() for i in range(t.size)]
    print("\n----------------------------------------------------------------------------------------------------\n")
    for r, j in zip(t, range(t.size)):
        valori, sigma = r.derivata_numerica(r.xdata, r.ydata, sigmax=r.sigmax, sigmay=r.sigmay, passo=passo)
        tau[j].ydata = np.append(tau[j].ydata, valori)
        tau[j].xdata = np.append(tau[j].xdata, r.ydata[:-passo])
        tau[j].add_sigmas(sigmax=r.sigmay[:-passo], sigmay=sigma)
        tau[j].reg_lin()
        tau[j].chi_quadro()
        print("Tau = ", -1/tau[j].B, " \t Ctot = ", -1/(tau[j].B*r.R_teorica), " \t R teorica = ", r.R_teorica)
    del tau    

'''
if enable_plots:
    tau_R.plotData()
    plt.show()
    tau_R.residui()
    plt.show()
'''