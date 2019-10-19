import scipy
import matplotlib
from matplotlib import pyplot as plt
import os
import csv

matplotlib.rcParams.update({'font.size': 22}) # Cambio dimensione font globale dei plot

nomefile = "lab1.csv" # Nome del file su cui salvare i dati (ed eventualmente caricare quelli esistenti)

# Array di storage dove tenere tutti i dati (che poi andranno anche salvati su un csv)
xdata = []
ydata = []
sigmadata = []

# Funzione di plotting, la faccio una volta per tutte
def plotxy(x , y):
    plt.plot(x, y, 'go--', linewidth=2, markersize=12, color=[0.01, 0.01, 0.01], alpha=0.5, markerfacecolor=[0, 0, 1])
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("X-Y Plot")
    plt.show()

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

#-------- Cerco se esiste già un file con questo nome ed eventualmente carico quei dati ----------#
myfile = os.path.join(nomefile)
if os.path.isfile(myfile):
    with open(nomefile, 'r') as csvFile:
        reader = csv.reader(csvFile)
        # Per ogni riga del file, aggiungo la tripletta all'array di storage
        for row in reader:
            xdata.append(row[0])
            ydata.append(row[1])
            sigmadata.append(row[2])
    plotxy(xdata, ydata)
#-------------------------------------------------------------------------------------------------#

while(True): # Arresto il ciclo quando vengono inseriti dei caratteri
    print("Inserire una tripla (x, y, sigma) nel formato 'dato1 dato2 dato3' per aggiungerla al grafico.")
    nuovo_dato = input()

    if(len(nuovo_dato.split()) == 3):
        x, y, sigma = nuovo_dato.split() # divido i valori in tre variabili diverse
        # Aggiungo i nuovi dati agli array di storage
        xdata.append(x)
        ydata.append(y)
        sigmadata.append(sigma)
        # Grafico
        plotxy(xdata, ydata)
    else:
        print("Input non riconosciuto")
        break