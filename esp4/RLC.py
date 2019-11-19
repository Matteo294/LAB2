from labbclass import FDT
from labbclass import Misura
from matplotlib import pyplot as plt 
import os
import sys
import numpy as np 
import csv

enable_plots = True
enable_simulation_plot = False
# E' possibile indicare da cmd quali grafici visualizzare 
# indicando un numero da 0 a n_resistenze - 1. Mettere x per non plottare
grafici_da_plottare = []
if len(sys.argv) > 1:
    for plot in sys.argv[1:]:
        if plot == 'x':
            enable_plots = False
        else:
            grafici_da_plottare.append(int(plot))

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

resistenze_files = ['dati/R1.csv', 'dati/R2.csv', 'dati/R3.csv'] # un file per ogni resistenza
resistenze = [Misura(9939, 2), Misura(5.10353e3, 2), Misura(467.34, 0.2)] # da cambiare !!
L = Misura(2.2e-3, 0.01) # da cambiare !!
C = Misura(48.2e-9, 0.2e-9) # incertezza da cambiare

# PLOT SIMULAZIONE
file_lettura = "simulazione.csv"
myfile = os.path.join(file_lettura)
freq = np.array([])
ampl = np.array([])
if os.path.isfile(myfile):
    with open(file_lettura, 'r') as csvFile:
        reader = csv.reader(csvFile)
        # Per ogni riga del file, leggo i valori frequenza e ampiezza
        for r in reader:
            row = [float(r[0]), float(r[1])] 
            freq = np.append(freq, row[0])
            ampl = np.append(ampl, row[1])
else:
    print("Problema: non trovo il file " + file_lettura)
if enable_simulation_plot:
    plt.plot(freq, ampl, '-')
    plt.show()

# Creo l'oggetto dalla classe FDT e leggo i dati dal file nel formato: freq Vout fase Vin
for R, f in zip(resistenze, resistenze_files):
    
    rlc = FDT()
    rlc.leggiDati(f)
    rlc.fdt_teorica(numeratore=[L.valore/R.valore, 0], denominatore=[L.valore*C.valore, L.valore/R.valore, 1]) # Plotto la fdt teorica per un circuito RLC ideale
    rlc.f_ris = 1 / (2*np.pi * np.sqrt(L.valore*C.valore))
    
    if enable_plots and resistenze.index(R) in grafici_da_plottare:
        rlc.plot_teorica_ampiezza()
        plt.semilogx(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), '.', markersize=10)
        plt.plot([rlc.f_ris, rlc.f_ris], [min(rlc._ampiezza_teo), max(rlc._ampiezza_teo)], '--', linewidth=1.8, color='red', label="Frequenza di risonanza")
        plt.legend()
        plt.grid()
        plt.show()

        rlc.plot_teorica_fase()
        plt.semilogx(rlc.freq, rlc.fase, '.', markersize=10)
        plt.plot([rlc.f_ris, rlc.f_ris], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color='red', label="Frequenza di risonanza")
        plt.legend()
        plt.grid()
        plt.show()
