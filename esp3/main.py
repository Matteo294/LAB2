import numpy as np
from matplotlib import pyplot as plt
from labbclass import FDT
from labbclass import Misura
import os
import sys

plot_grafici = True # Se True stampa tutti i grafici

# E' possibile indicare da cmd quali visualizzare
grafici_da_plottare = []
if len(sys.argv) > 1:
    for plot in sys.argv:
        grafici_da_plottare.append(int(plot))
else:
    grafici_da_plottare = [2]

# Resistenze e capacità usate: valore misurato e incertezza (vedi classe Misura in labbclass.py)
R = [Misura(146.27, 0.02), Misura(1600, 0.2), Misura(16000, 2)]
C = Misura(35.5e-9, 0.2e-9)

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Ciascun elemento dell'array contiene i file relativi alla corrispettiva resistenza
files = ['dati/R1.csv', 'dati/R2.csv', 'dati/R3.csv']

# Creo un array con le varie funzioni di trasferimento e leggo i dati dai file
fdt = np.array([])
for r,i in zip(R, range(len(R))):
    fdt = np.append(fdt, FDT())
    if i == 2:
        fdt[i].leggiDati(files[i])
    else:
        fdt[i].leggiDati(files[i], scale_f = 10e3)
    fdt[i].fase = -fdt[i].fase
    fdt[i].tau = r.valore*C.valore

# Calcolo le fdt teoriche
for H in fdt:
    H.fdt_teorica(numeratore=[1], denominatore=[H.tau, 1])

# Plotto solo le fdt richieste
if plot_grafici is not False:
    for i in grafici_da_plottare:

        # Plot ampiezze misurate sul grafico teorico
        fdt[i].plot_teorica_ampiezza()
        plt.semilogx(fdt[i].freq, 20*np.log10(fdt[i].Vout / fdt[i].Vin), '.')
        plt.grid()
        plt.show()

        # Plot fasi misurate sul grafico teorico
        fdt[i].plot_teorica_fase()
        plt.semilogx(fdt[i].freq, fdt[i].fase, '.')
        plt.grid()
        plt.show()
    
    
