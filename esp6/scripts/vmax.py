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
enable_plots = 1
if len(sys.argv) > 1:
    enable_plots = int(sys.argv[1])

# Impostazioni per i grafici
plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=18)

# File delle misure
file_soloC = 'Misure/soloC.csv'
file_dmm = 'Misure/dmm.csv'
file_zener = 'Misure/zener.csv'

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6

graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
graetz.Vout = graetz.leggi_colonna(file_soloC, 1)
graetz.I = graetz.Vout / graetz.resistenze

# Caratteristica del diodo
def Vdiodo(i):
    if i > 0:
        return 0.9096 + 0.04596*np.log(i)
    elif i == 0:
        return 0
Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])

# Funzione ausiliaria
func = lambda V, R: V0 - 2*Vdiodo_vettorizzata(V/R) - V
func_vettorizzata = np.vectorize(func, [float])

graetz.Vmax = np.array([])
for R in graetz.resistenze:
    graetz.Vmax = np.append(graetz.Vmax, graetz.risolvi_numericamente(func, 5, 12, nsteps=10000, param1=R))
for R, Vout, Vmax in zip(graetz.resistenze, graetz.Vout, graetz.Vmax):
    print("R: {0:.1f}   \t Vmax: {1:.4f} \t Vout: {2:.4f}".format(R, Vmax, Vout)) 
print(np.array([graetz.resistenze, graetz.Vmax]))
graetz.scriviDati(file_scrittura='Vmax_teorica_no_zener.csv', dati=np.array([graetz.resistenze, graetz.Vmax]))

# Vanno messe le barre d'errore
if enable_plots:
    plt.plot(graetz.Vout/graetz.resistenze, graetz.Vout, '.', label=r"$V^{max}$ misurata", linewidth=2.0)
    plt.plot(graetz.Vout/graetz.resistenze, graetz.Vmax, '.', label=r"$V^{max}$ calcolata", linewidth=2.0)
    plt.xlabel(r"$i_L$ [$A$]")
    plt.ylabel(r"$V^{max}$ [V]")
    plt.title("Tensione massima in uscita dal circuito")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

