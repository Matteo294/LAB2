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
print_risultati = True
plot_finale_ripple = True


# Impostazioni per i grafici
plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=18)

# File delle misure
file_soloC = '../Misure/soloC.csv'
file_dmm = '../Misure/dmm.csv'
file_zener = '../Misure/zener.csv'
file_Vmax_no_zener = 'Vmax_teorica_no_zener.csv'

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6

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
    return np.abs(Vin_vettorizzata(t)) - 2*Vdiodo_vettorizzata(Vc_vettorizzata(t, R, Vmax)/R) - Vc_vettorizzata(t, R, Vmax)

# Vettorizzazione le funzioni delle tensioni: utile quando si devono applicare ad un array di valori
Vc_vettorizzata = np.vectorize(Vc, [float])
Vin_vettorizzata = np.vectorize(Vin, [float])
Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])
func_vettorizzata = np.vectorize(func, [float])

#--------------------------- Inizio analisi -------------------------------------
graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
graetz.ripple = graetz.leggi_colonna(file_soloC, 3)
graetz.Vmax = graetz.leggi_colonna(file_Vmax_no_zener, 1)
graetz.Vmax_sperimentale = graetz.leggi_colonna(file_soloC, 1)
graetz.sigmaVripple = graetz.leggi_colonna(file_soloC, 2)
graetz.sigmaVripple = graetz.sigmaVripple * 24/100

# Inizializzo gli array di storage
graetz.ripple_teo = np.array([])
graetz.dV_teo = np.array([])
graetz.dV_sperimentale = np.array([])

for R, i, Vmax in zip(graetz.resistenze, range(len(graetz.resistenze)), graetz.Vmax):

    tau_C = R*C
    # Risoluzione numerica dell'equazione
    t0 = graetz.risolvi_numericamente(func_vettorizzata, 1/2*np.pi/w, np.pi/w, nsteps=10000, param1=R, param2=Vmax)
    Vmin = Vc(t0, R, Vmax)

    # Aggiungo i risultati agli array di storage
    graetz.ripple_teo = np.append(graetz.ripple_teo, Vmax - Vmin)
    graetz.dV_sperimentale = np.append(graetz.dV_sperimentale, graetz.ripple[i] / graetz.Vmax_sperimentale[i] * 100)
    graetz.dV_teo = np.append(graetz.dV_teo,  (Vmax - Vmin) / Vmax * 100)

    # Print dei risultati
    if print_risultati:
        print("-----------------------------------------------------------------------------------------------------------------------")
        print("Resistenza: ", R)
        print("Vmax teorica: ", Vmax)
        print("tempo di scarica: ", t0)
        print("Tensione teorica a fine scarica: ",  Vmin)
        print("Ripple calcolato: ", Vmax - Vmin,  "\t Ripple misurato: ", graetz.ripple[i]) # Aggiungere incertezze
        print("Ripple percentuale calcolato: ", graetz.dV_teo[i], "%\tRipple percentuale sperimentale: ", graetz.dV_sperimentale[i], "%")
        print("-----------------------------------------------------------------------------------------------------------------------")
    
    # Vanno messe le barre d'errore
    if i in da_plottare:
        # Grafico cos'è successo
        t = np.linspace(0, np.pi/w, 1000)
        plt.plot(t, Vc_vettorizzata(t, R, Vmax), label="$V_c$", linewidth=2, color=[1, 0.5, 0], alpha=0.7)
        plt.plot(t, np.abs(Vin_vettorizzata(t)) - 2*Vdiodo_vettorizzata(Vc_vettorizzata(t, R, Vmax)/R), label=r"$\left| V_{in}\right| - 2V_D(i(t))$", linewidth=2, alpha=0.7, color="firebrick")
        plt.plot(t0, Vmin, 'o', fillstyle='full', markersize = 8, color="gray", label=r'$V^{min}$')
        plt.title("Grafico tensioni")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Tensione [V]")
        plt.legend(loc = 'lower left')
        plt.grid()
        plt.show()
        print("\n")

if plot_finale_ripple:
    print(graetz.sigmaVripple)
    print(graetz.ripple)
    plt.errorbar(graetz.ripple/graetz.resistenze, graetz.ripple, yerr=graetz.sigmaVripple, marker = '.',  markersize=8, ecolor ='lightgray', color='royalblue', linestyle = '', label="Ripple misurato")
    plt.semilogx(graetz.ripple/graetz.resistenze, graetz.ripple_teo, '.', markersize=8, color = 'orange', label="Ripple calcolato")
    plt.xlabel(r"$i_L$ [A]")
    plt.ylabel("Ripple [V]")
    plt.title("Tensioni di ripple")
    plt.legend(loc='upper left')
    plt.grid()
    plt.show()

