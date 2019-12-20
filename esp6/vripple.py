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

# Da terminale si possono selezionare i grafici da plottare, con numeri da 0 a 9
da_plottare = []
for nplot in sys.argv[1:]:
    da_plottare.append(int(nplot))

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
Vd = 0.65 # Tensione su ciascun diodo
Vp = V0 - 2*Vd # Vmax in uscita
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6
t0 = np.arcsin(2*Vd/V0) / w

Vin = lambda t: V0*np.cos(w*t) # Tensione in ingresso al ponte

def Vc(t, R, Vmax):
    V = Vmax*np.exp(-t/(R*C)) # Scarica del condensatore
    return V

# Curva caratteristica del diodo ricostruita sperimentalmente
def Vdiodo(i):
    if i > 0:
        return 0.9096 + 0.04596*np.log(i)
    elif i == 0:
        return 0
# Funzione ausiliaria per la linea di carico
# Parametro 1: R, parametro 2: Vmax
def func(t, param1, param2): 
    R = param1
    Vmax = param2
    print(Vmax)
    return Vin_vettorizzata(t) - 2*Vdiodo_vettorizzata(Vc_vettorizzata(t, R)/R) - Vc_vettorizzata(t, R, Vmax)/R

# Vettorizzazione le funzioni delle tensioni: utile quando si devono applicare ad un array di valori
Vc_vettorizzata = np.vectorize(Vc, [float])
Vin_vettorizzata = np.vectorize(Vin, [float])
Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])
func_vettorizzata = np.vectorize(func, [float])

#--------------------------- Inizio analisi -------------------------------------
graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
print("Resistenze: ", graetz.resistenze, '\n\n')
graetz.ripple = graetz.leggi_colonna(file_soloC, 3)
graetz.ripple_teo = np.array([])
graetz.Vmax = graetz.leggi_colonna(file_soloC, 1)

for R, i, Vmax in zip(graetz.resistenze, range(len(graetz.resistenze)), graetz.Vmax):

    tau_C = R*C
    # Risoluzione numerica dell'equazione
    t0 = graetz.risolvi_numericamente(func_vettorizzata, 3/2*np.pi/w, 2*np.pi/w, nsteps=10000, param1=R, param2=Vmax)
    Vmin = Vc(t0, R, Vmax)

    # Print dei risultati
    print("Resistenza: ", R)
    print("Tensione a fine scarica: ",  Vmin)
    print("Ripple calcolato: ", Vc(t0, R, Vmax) - Vmin,  "\t Ripple misurato: ", graetz.ripple[i]) # Aggiungere incertezze

    graetz.ripple_teo = np.append(graetz.ripple_teo, Vc(t0, R, Vmax) - Vmin)
    
    # Vanno messe le barre d'errore
    if i in da_plottare:
        # Grafico cos'è successo
        t = np.linspace(np.pi/2/w + t0, 2*np.pi/w, 1000)
        #plt.plot(t, Vponte_vettorizzata(t), label="$V_{ponte}$", linewidth=2, color=[0, 0, 0.9], alpha=0.7)
        plt.plot(t, Vc_vettorizzata(t, tau_C), label="$V_c$", linewidth=2, color=[1, 0.5, 0], alpha=0.7)
        plt.plot(t, Vin_vettorizzata(t), label="$V_{in}$", linewidth=2, alpha=0.7, color="gray")
        plt.plot([t, t], [-V0, Vmin], '--', linewidth=2, color="gray")
        plt.title("Grafico tensioni")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Tensione [V]")
        plt.legend()
        plt.grid()
        plt.show()
    
    print("\n")

plt.semilogx(graetz.resistenze, graetz.ripple, '.',  markersize=16, label="Ripple misurato")
plt.semilogx(graetz.resistenze, graetz.ripple_teo, '.', markersize=16, label="Ripple calcolato")
plt.xlabel(r"Resistenza [$\Omega$]")
plt.ylabel("Ripple [V]")
plt.title("Tensioni di ripple")
plt.legend()
plt.grid()
plt.show()
