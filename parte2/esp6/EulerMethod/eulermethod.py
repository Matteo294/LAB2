''' Metodo di Eulero per la risoluzione di equazioni del tipo ydot = f(t, y) 
Guida per l'uso:
- Creare una funzione (func in questo esempio) che rappresenti f(t, y)
- Lanciare lo script da terminale specificando, nell'ordine, i parametri: tempo di arresto, numero di punti, y(t=0)
- I risultati saranno disponibili file "risultati.csv" all'interno della cartella
'''

import numpy as np 
from labbclass import Analisi
from matplotlib import pyplot as plt 
import sys 
import csv
import os


# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

filename_base = 'risultati'

# Impostazioni per i grafici
plt.rc('text', usetex=True) 
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
plt.rc('font', size=18)

# Costanti del problema
A = 0.9096
B = 0.04596
C = 220e-6
R = 100
C1 = lambda RL: C * (R + RL + RL*B*R)
tau = lambda RL: C*(R + RL + RL*B*R) / (1 + B*RL)

file_resistenze = '../Misure/zener.csv'
zener = Analisi()
zener.resistenze = zener.leggi_colonna(file_resistenze, colonna=0)
zener.Vc_max = zener.leggi_colonna(file_resistenze, colonna=3)

# Parametri passati da terminale
tstop = float(sys.argv[1])
npunti = int(sys.argv[2])
#y0 = float(sys.argv[3])

# Step temporale e array dei tempi
dt = tstop/npunti 
t = np.arange(0, tstop, dt)

# Funzione f(t, y) membro di destra dell'equazione differenziale
func = lambda t, y, RL: (-y - RL*B*y - RL*A)/C1(RL)

# Avanza di uno step nel tempo
def avanza(f, t, dt, y0, param=None):
    if param is None:
        y1 = y0 + dt*f(t, y0)
    else:
        y1 = y0 + dt*f(t, y0, param)
    return y1

def salva_su_file(fname, dati):
    myfile = os.path.join(fname)
    with open(myfile, 'w', newline='\n') as csvfile:
        mywriter = csv.writer(csvfile, delimiter = ',')
        mywriter.writerows(np.transpose(dati))

def main():
    for RL, Vmax, i in zip(zener.resistenze, zener.Vc_max, range(len(zener.Vc_max))):
        y0 = Vmax
        y = np.array([y0])
        for n in range(npunti):
            y = np.append(y, avanza(func, t[n], dt, y[n], param=RL))

        filename = filename_base + str(i) + '.csv'
        salva_su_file(filename, [t, y[:-1]])

        analitico = Vmax*np.exp(-t/tau(RL)) + A*RL/(1 + B*RL) * (np.exp(-t/tau(RL)) - 1)

        # Grafico di y(t)
        plt.plot(t, y[:-1], linewidth=2.0, color='black', alpha=0.8, label="y(t) numerico")
        plt.plot(t, analitico, linewidth=2.0, alpha=0.8, label="y(t) analitico")
        plt.xlabel('t')
        plt.ylabel('y(t)')
        plt.title(f"Soluzione dell'equazione differenziale con condizione iniziale $y_0 = {y0}$")
        plt.legend()
        plt.grid()
        plt.show()

        # Grafico differenza
        plt.plot(t, (y[:-1] - analitico) / analitico, linewidth=2.0, color='black', alpha=0.8, label="y(t) numerico")
        plt.xlabel('t')
        plt.ylabel('err')
        plt.title(f"Differenza tra risulato numerico e analitico con {npunti} punti")
        plt.legend()
        plt.grid()
        plt.show()

if __name__ == "__main__":
    main()