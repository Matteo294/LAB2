import numpy as np
from matplotlib import pyplot as plt
from labbclass import FDT
import os

R_da_plottare = [0] # Comando brutto per scagliere quali plot visualizzare

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Ciascun array contiene i file relativi alla corrispettiva resistenza
files = [['R1_1.csv', 'R1_2.csv'], ['R2_1.csv', 'R2_2.csv'], ['R3_1.csv', 'R3_1.csv']]
cartella = 'dati/'

# Resistenze e capacità usate
R = [53000, 5300, 530]
C = 10e-8

# Creo un array con le varie funzioni di trasferimento e leggo i dati dai file
fdt = np.asarray([])
for r,i in zip(R, range(len(R))):
    fdt = np.append(fdt, FDT())
    for f in files[i]:
        fdt[i].leggiDati(cartella+f)
    fdt[i].tau = r*C

# Calcolo le fdt teoriche
for H in fdt:
    H.fdt_teorica(numeratore=[1], denominatore=[H.tau, 1])

# Plotto solo le fdt richieste
for i in R_da_plottare:
    fdt[i].plot_teorica_ampiezza()
    plt.semilogx(fdt[i].omega, fdt[i].ydata / fdt[i].xdata)
    plt.grid()
    plt.show()
    
    
