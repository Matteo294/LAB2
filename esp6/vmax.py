''' Questo script serve per calcolare le tensioni di ripple e produrre un grafico per spiegare come si sono calcolate '''

#from sympy import Eq, plot, symbols, solveset, sin, init_printing, Interval
import math
import numpy as np
from matplotlib import pyplot as plt
from labbclass import Analisi
import sys 
import os

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Da terminale si può selezionare se plottare i grafici: 0 -> non plottare, 1 -> plotta
enable_plots = 0
if len(sys.argv) > 1:
    enable_plots = int(argv[1])

# Impostazioni per i grafici
plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=22)

# File delle misure
file_soloC = 'Misure/soloC.csv'
file_dmm = 'Misure/dmm.csv'
file_zener = 'Misure/zener.csv'

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
Vd = 0.7 # Tensione su ciascun diodo
Vp = V0 - 2*Vd # Valore di picco dopo il ponte (caduta su due diodi)
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6
t0 = np.arcsin(2*Vd/V0) / w

graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
graetz.Vout = graetz.leggi_colonna(file_soloC, 2)
graetz.I = graetz.Vout / graetz.resistenze

# !!!!!!!!!! Caratteristica del diodo !!!!!!!!!!
# # Vdiodo = lambda i: e^i....                          <----------------- # DA INSERIRE!!: vanno passati anche gli altri parametri ma non so quali siano
# Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])

graetz.Vmax = V0 - 2*Vdiodo_vettorizzata(graetz.I)

for Vmax, Vout in zip(graetz.Vout, graetz.Vmax):
    print(f"Vmax: {Vmax} \t Vout: {Vout}") # Aggiungere incertezze

# Vanno messe le barre d'errore
if enable_plots:
    plt.plot(graetz.resistenze, graetz.Vout, '.', label="Vout misurata", linewidth=2.0)
    plt.plot(graetz.resistenze, graetz.Vmax, '.', label="Vmax calcolata", linewidth=2.0)
    plt.xlabel(r"Resistenza [$\Omega$]")
    plt.ylabel("Vmax [V]")
    plt.title("Tensione massima in uscita dal ponte di Graetz")
    plt.label()
    plt.grid()
    plt.show()

