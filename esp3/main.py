import numpy as np
from matplotlib import pyplot as plt
from labbclass import FDT
from labbclass import Misura
import os
import sys

plot_grafici = False

# Quando si runna da riga di comando si può decidere se plottare o meno aggiungendo True/False
if len(sys.argv > 1):
    plot_grafici = sys.argv[1]

R_da_plottare = [0] # Comando brutto per scagliere quali plot visualizzare
Vpp1 = [1.51, 1.757, 1.953, 1.917, 1.592, 1.544]

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Ciascun array contiene i file relativi alla corrispettiva resistenza
files = ['R1.csv', 'R2.csv', 'R3.csv']
cartella = 'dati/'

# Resistenze e capacità usate: valore misurato e incertezza
R = [Misura(146.27, 0.02), Misura(1600, 0.2), Misura(16000, 2)]
C = Misura(35.5e-9, 0.2e-9)

# Creo un array con le varie funzioni di trasferimento e leggo i dati dai file
fdt = np.asarray([])
for r,i in zip(R, range(len(R))):
    fdt = np.append(fdt, FDT())
    for f in files[i]:
        fdt[i].leggiDati(cartella+f, scale_x=1e3) # Cambio udm da kOhm a Ohm
        fdt[i].setVin(Vpp1)
    fdt[i].tau = r.val*C.val

# Calcolo le fdt teoriche
for H in fdt:
    H.fdt_teorica(numeratore=[1], denominatore=[H.tau, 1])

# Plotto solo le fdt richieste
if plot_grafici is not False:
    for i in R_da_plottare:
        fdt[i].plot_teorica_ampiezza()
        plt.semilogx(fdt[i].omega, 20*np.log10(fdt[i].Vout / fdt[i].Vin), '.')
        plt.grid()
        plt.show()
        fdt[i].plot_teorica_fase()
        plt.semilogx(fdt[i].omega, fdt[i].fase, '.')
        plt.show()
    
    
