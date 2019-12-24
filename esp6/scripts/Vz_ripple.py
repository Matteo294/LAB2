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
print_risultati = 0
plot_finale_ripple = True

# Impostazioni per i grafici
plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=18)

# File delle misure
file_zener = '../Misure/zener.csv'
file_Vmax_zener = 'Vmax_teorica_zener.csv'

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6
Rz = 100
# Coefficienti zener
Az = -1.81110
Bz = 0.35762
# Coefficienti diodo
Ad = 0.9096
Bd = 0.04596

Vin = lambda t: V0*np.cos(w*t) # Tensione in ingresso al ponte
Vin_vettorizzata = np.vectorize(Vin, [float])

# Coefficiente comodo
R0 = lambda RL: Bz*Rz*RL + Rz + RL

# Tau scarica condensatore
tau = lambda RL: C*R0(RL)/(1+Bz*RL)

def Vc(t, RL, Vc0):
    V = Vc0*np.exp(-t/tau(RL)) - Az*RL/(1 + Bz*RL)*(1 - np.exp(-t/tau(RL)))
    return V
Vc_vettorizzata = np.vectorize(Vc, [float])

def Vout(t, RL):
    return RL/R0(RL) * Vc(t, RL, Vc0) - Rz*RL/R0(RL)*Az
Vout_vettorizzata = np.vectorize(Vout, [float])

# Caratteristica del diodo
def Vdiodo(i):
    if i > 0:
        return Ad + Bd*np.log(i)
    elif i == 0:
        return 0
Vdiodo_vettorizzata = np.vectorize(Vdiodo, [float])

# Caratteristica zener
def iZener(V):
    return Bz*V - Az
iZener_vettorizzata = np.vectorize(iZener, [float])

def Vzener(iz):
    return (iz + Az)/Bz

# Funzione ausiliaria per la linea di carico
# Parametro 1: R, parametro 2: Vmax
def func(t, param1, param2): 
    RL = param1
    Vc0 = param2
    return np.abs(Vin_vettorizzata(t)) - 2*Vdiodo_vettorizzata(Vout_vettorizzata(t, RL)/RL + iZener(Vout_vettorizzata(t, RL))) - Vc_vettorizzata(t, RL, Vc0)
func_vettorizzata = np.vectorize(func, [float])

#--------------------------- Inizio analisi -------------------------------------
graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_zener, 0)
graetz.ripple = graetz.leggi_colonna(file_zener, 5)
graetz.Vc0 = graetz.leggi_colonna(file_zener, 3)
graetz.Vmax = graetz.leggi_colonna(file_Vmax_zener, 1)
graetz.sigmaVripple = graetz.leggi_colonna(file_zener, 6)
graetz.sigmaVripple = graetz.sigmaVripple * 32/100 *math.sqrt(2)

# Inizializzo gli array di storage
graetz.ripple_teo = np.array([])

for RL, i, Vc0, Vmax in zip(graetz.resistenze, range(len(graetz.resistenze)), graetz.Vc0, graetz.Vmax):

    # Risoluzione numerica dell'equazione
    if math.isinf(RL):
        Vmin = Vmax
    else:
        t0 = graetz.risolvi_numericamente(func_vettorizzata, 1/2*np.pi/w + 0.0005, np.pi/w, nsteps=10000, param1=RL, param2=Vc0)
        Vmin = Vout(t0, RL)

    # Aggiungo i risultati agli array di storage
    graetz.ripple_teo = np.append(graetz.ripple_teo, Vmax - Vmin)

    # Print dei risultati
    if print_risultati:
        print("-----------------------------------------------------------------------------------------------------------------------")
        print("Resistenza: ", RL)
        print("Vmax teorica: ", Vmax)
        print("tempo di scarica: ", t0)
        print("Tensione teorica a fine scarica: ",  Vmin)
        print("Ripple calcolato: ", Vmax - Vmin,  "\t Ripple misurato: ", graetz.ripple[i], "\u00B1", graetz.ripple[i]*0.3) # Aggiungere incertezze
        #print("Ripple percentuale calcolato: ", graetz.dV_teo[i], "%\tRipple percentuale sperimentale: ", graetz.dV_sperimentale[i], "%")
        print("-----------------------------------------------------------------------------------------------------------------------")
    
    # Vanno messe le barre d'errore
    if i in da_plottare:
        # Grafico cos'è successo
        t = np.linspace(0, np.pi/w, 1000)
        plt.plot(t, Vc_vettorizzata(t, RL, Vc0), label="$V_c$", linewidth=2, color=[1, 0.5, 0], alpha=0.7)
        plt.plot(t, np.abs(Vin_vettorizzata(t)) - 2*Vdiodo_vettorizzata(Vout_vettorizzata(t, RL)/RL + iZener(Vout_vettorizzata(t, RL))), label=r"$\left| V_{in}\right| - 2V_D(i(t))$", linewidth=2, alpha=0.7, color="firebrick")
        plt.plot(t0, Vmin, 'o', fillstyle='full', markersize = 8, color="gray", label=r'$V^{min}$')
        plt.title("Grafico tensioni")
        plt.xlabel("Tempo [s]")
        plt.ylabel("Tensione [V]")
        plt.legend(loc = 'lower left')
        plt.grid()
        plt.show()
        print("\n")

if plot_finale_ripple:
    plt.errorbar(graetz.resistenze, graetz.ripple, yerr=graetz.sigmaVripple, marker = '.',  markersize=8, ecolor ='gray', color='royalblue', linestyle = '', label="Ripple misurato")
    plt.semilogx(graetz.resistenze, graetz.ripple_teo, '.', markersize=8, color = 'orange', label="Ripple calcolato")
    plt.xlabel(r"$R_L~[\Omega]$")
    plt.ylabel("Ripple [V]")
    plt.title("Tensioni di ripple")
    plt.legend(loc='upper right')
    plt.grid()
    plt.show()

print(graetz.sigmaVripple)