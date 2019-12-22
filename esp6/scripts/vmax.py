''' Questo script serve per calcolare le tensioni di ripple e produrre un grafico per spiegare come si sono calcolate '''

#from sympy import Eq, plot, symbols, solveset, sin, init_printing, Interval
import math
import numpy as np
from matplotlib import pyplot as plt
from labbclass import Analisi
from labbclass import Misura
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
file_soloC = '../Misure/soloC.csv'
file_dmm = '../Misure/dmm.csv'
file_zener = '../Misure/zener.csv'

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6
A = Misura(9096e-4, 3e-4)
B = Misura(4596e-5, 6e-5)

graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
print(graetz.resistenze)
graetz.sigma_resistenze = graetz.resistenze/100
print(graetz.sigma_resistenze)
graetz.Vout = graetz.leggi_colonna(file_soloC, 1)
graetz.sigmaVout = graetz.leggi_colonna(file_soloC,2)
graetz.sigmaVout = graetz.sigmaVout * 24 / 100

# Caratteristica del diodo
def Vdiodo(i):
    if i > 0:
        return A.valore + B.valore*np.log(i)
    elif i == 0:
        return 0
Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])

# Funzione ausiliaria
func = lambda V, R: V0 - 2*Vdiodo_vettorizzata(V/R) - V
func_vettorizzata = np.vectorize(func, [float])

graetz.Vmax = np.array([])
graetz.sigma_Vmax = np.array([])
for R, dR in zip(graetz.resistenze, graetz.sigma_resistenze):
    VM = graetz.risolvi_numericamente(func, 5, 12, nsteps=10000, param1=R)
    graetz.Vmax = np.append(graetz.Vmax, VM)
    graetz.sigma_Vmax = np.append(graetz.sigma_Vmax, 2 * np.sqrt( (A.sigma)**2 + (np.log(VM/R)*B.sigma)**2 + (B.valore/R*dR)**2))
for R, dR, Vout, dVout, Vmax, dVmax in zip(graetz.resistenze, graetz.sigma_resistenze, graetz.Vout, graetz.sigmaVout, graetz.Vmax, graetz.sigma_Vmax):
    print("R: {0:.1f} \u00B1 {1:.1f}  \t Vmax: {2:.4f} \u00B1 {3:.4f} \t Vout: {4:.4f} \u00B1 {5:.4f}".format(R, dR, Vmax, dVmax, Vout, dVout)) 
#print(np.array([graetz.resistenze, graetz.Vmax]))


# Vanno messe le barre d'errore
if enable_plots:
    plt.errorbar(graetz.Vout/graetz.resistenze, graetz.Vout, yerr=graetz.sigmaVout, marker='.', color = 'royalblue', ecolor = 'lightgray', linestyle='', label=r"$V^{max}$ misurata",linewidth=2.0, markersize=8)
    plt.errorbar(graetz.Vout/graetz.resistenze, graetz.Vmax, yerr=graetz.sigma_Vmax, marker='.', color = 'orange', linestyle='', label=r"$V^{max}$ calcolata", markersize=8)
    plt.xlabel(r"$i_L$ [$A$]")
    plt.ylabel(r"$V^{max}$ [V]")
    plt.title("Tensione massima in uscita dal circuito")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

