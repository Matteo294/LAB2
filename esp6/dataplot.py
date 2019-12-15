''' Questo script serve per costruire il plot necessario per spiegare le modalità con cui è stata calcolata
la tensione di ripple, nel circuito del ponte di Graetz, senza diodo zener '''

#from sympy import Eq, plot, symbols, solveset, sin, init_printing, Interval
import math
import numpy as np
from matplotlib import pyplot as plt
from labbclass import Analisi

plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=22)

# Costanti
Vrms = 7.5 # V efficace all'uscita dal trasformatore
V0 = Vrms * math.sqrt(2) # Valore di picco
Vd = 0.65 # Tensione su ciascun diodo
Vp = V0 - 2*Vd # Valore di picco dopo il ponte (caduta su due diodi)
f = 50 # Frequenza
w = 2*math.pi*f # Pulsazione
T = 1/f # Periodo
t0 = 3.9112e-4 # Tempo al quale Vin(t) = 2Vd
C = 200e-6
R = 1020
tau_C = R*C

# Tensione dopo il ponte. Quando Vin < 2Vd, la tensione in uscita dal ponte è 0
def Vponte(t):
    V = 0
    if t < t0:
        V = 0
    elif t < (T/2 - t0):
        V = Vp*np.abs(np.sin(w*t))
    elif t < (T/2 + t0):
        V = 0
    elif t < (T - t0):
        V = Vp*np.abs(np.sin(w*t))
    return V

def Vc(t):
    if t < np.pi/2/w:
        V = Vponte(t)
    else:
        V = Vp*np.exp(-(t-(np.pi/2/w))/tau_C) # Scarica del condensatore
    return V
    
Vin = lambda t: V0*np.sin(w*t) # Tensione in ingresso al ponte

def func(t):
    return Vponte(t) - Vc(t)

V_vettorizzata = np.vectorize(func, [float]) # Questa funzione vettorizza la chiamata della funzione qualora si dovesse usare per un array di parametri.
Vponte_vettorizzata = np.vectorize(Vponte, [float])
Vc_vettorizzata = np.vectorize(Vc, [float])
Vin_vettorizzata = np.vectorize(Vin, [float])

# Risoluzione numerica dell'equazione 
graetz = Analisi()
delta_t = graetz.risolvi_numericamente(V_vettorizzata, np.pi/w, 3/2*np.pi/w)

# Print dei risultati
print("Tempo di scarica: ", delta_t)
print("Tensione a fine scarica: ",  Vc(delta_t))
print("Tensione di ripple: ", Vc(np.pi/2/w) - Vc(delta_t))

# Grafico cos'è successo
t = np.linspace(np.pi/2/w + t0, 2*np.pi/w, 1000)
plt.plot(t, Vponte_vettorizzata(t), label="$V_{ponte}$", linewidth=1.8)
plt.plot(t, Vc_vettorizzata(t), label="$V_c$", linewidth=1.8)
plt.plot(t, Vin_vettorizzata(t), label="$V_{in}$", linewidth=1.8)
plt.plot([delta_t, delta_t], [0, Vc(delta_t)], '--', linewidth=1.8)
plt.plot([0, 2*np.pi/w], [Vc(np.pi/2/w), Vc(np.pi/2/w)])
plt.plot([0, 2*np.pi/w], [Vc(delta_t), Vc(delta_t)])
plt.title("Grafico tensioni")
plt.xlabel("Tempo [s]")
plt.ylabel("Tensione [V]")
plt.legend()
plt.grid()
plt.show()
