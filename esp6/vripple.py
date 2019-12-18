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
Vd = 0.7 # Tensione su ciascun diodo
Vp = V0 - 2*Vd # Valore di picco dopo il ponte (caduta su due diodi)
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
C = 200e-6
t0 = np.arcsin(2*Vd/V0) / w

Vin = lambda t: V0*np.sin(w*t) # Tensione in ingresso al ponte

# Tensione dopo il ponte. Quando Vin < 2Vd, la tensione in uscita dal ponte è 0
def Vponte(t):
    V = 0
    if t < t0:
        V = 0
    elif t < (T/2 - t0):
        V = np.abs(Vin(t)) - 2*Vd
    elif t < (T/2 + t0):
        V = 0
    elif t < (T - t0):
        V = np.abs(Vin(t)) - 2*Vd
    return V

def Vc(t, tau):
    if t < np.pi/2/w:
        V = Vponte(t)
    else:
        V = Vp*np.exp(-(t-(np.pi/2/w))/tau) # Scarica del condensatore
    return V

# Funzioni ausiliarie comode per la risoluzione numerica delle equazioni
def func(t, tau):
    return Vponte(t) - Vc(t, tau)
def func2(t):
    return Vin(t) - 2*Vd

V_vettorizzata = np.vectorize(func, [float]) # Questa funzione vettorizza la chiamata della funzione qualora la si dovesse usare per un array di parametri.
Vponte_vettorizzata = np.vectorize(Vponte, [float])
Vc_vettorizzata = np.vectorize(Vc, [float])
Vin_vettorizzata = np.vectorize(Vin, [float])
deltaV_vettorizzata = np.vectorize(func2, [float])


#--------------------------- Inizio analisi -------------------------------------
graetz = Analisi()
graetz.resistenze = graetz.leggi_colonna(file_soloC, 0)
print("Resistenze: ", graetz.resistenze, '\n\n')
graetz.ripple = graetz.leggi_colonna(file_soloC, 3)
graetz.ripple_teo = np.array([])

for R, i in zip(graetz.resistenze, range(len(graetz.resistenze))):

    tau_C = R*C

    # Risoluzione numerica dell'equazione
    #t0 = graetz.risolvi_numericamente(deltaV_vettorizzata, 0, np.pi/4/w, nsteps=10000) # Tempo al quale Vin(t) = 2Vd
    #print("t0: ", t0)

    # Risoluzione numerica dell'equazione
    delta_t = graetz.risolvi_numericamente(V_vettorizzata, np.pi/w, 3/2*np.pi/w, nsteps=10000, params=tau_C)

    # Print dei risultati
    print("Resistenza: ", R)
    print("Tempo di scarica: ", delta_t)
    print("Tensione a fine scarica: ",  Vc(delta_t, tau_C))
    print("Ripple calcolato: ", Vc(np.pi/2/w + t0, tau_C) - Vc(delta_t, tau_C), "\t Ripple misurato: ", graetz.ripple[i]) # Aggiungere incertezze

    graetz.ripple_teo = np.append(graetz.ripple_teo, Vc(np.pi/2/w + t0, tau_C) - Vc(delta_t, tau_C))
    
    # Vanno messe le barre d'errore
    if i in da_plottare:
        # Grafico cos'è successo
        t = np.linspace(np.pi/2/w + t0, 2*np.pi/w, 1000)
        plt.plot(t, Vponte_vettorizzata(t), label="$V_{ponte}$", linewidth=2, color=[0, 0, 0.9], alpha=0.7)
        plt.plot(t, Vc_vettorizzata(t, tau_C), label="$V_c$", linewidth=2, color=[1, 0.5, 0], alpha=0.7)
        plt.plot(t, Vin_vettorizzata(t), label="$V_{in}$", linewidth=2, alpha=0.7, color="gray")
        plt.plot([delta_t, delta_t], [-V0, Vc(delta_t, tau_C)], '--', linewidth=2, color="gray")
        plt.plot([delta_t, delta_t], [Vponte(delta_t), V0 - 2*Vd], color="Black", linewidth=2, label='Ripple')
        #plt.plot([0, 2*np.pi/w], [Vponte(np.pi/2/w + t0), Vponte(np.pi/2/w + t0)])
        #plt.plot([0, 2*np.pi/w], [Vc(delta_t, tau_C), Vc(delta_t, tau_C)])
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
