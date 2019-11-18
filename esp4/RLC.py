from labbclass import FDT
from labbclass import Misura
from matplotlib import pyplot as plt 
import os
import sys
import numpy as np 

enable_plots = True

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

resistenze_files = ['dati/R1.csv'] # un file per ogni resistenza
resistenze = [Misura(160, 0.02)] # da cambiare !!
L = Misura(10, 0.01) # da cambiare !!
C = Misura(35.5e-9, 0.2e-9)

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
        plt.semilogx(rlc.freq, 20*np.log10(rlc.Vout / rlc.fase), '.', markersize=10)
        plt.plot([rlc.f_ris, rlc.f_ris], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color='red', label="Frequenza di risonanza")
        plt.legend()
        plt.grid()
        plt.show()