''' 
Questo script viene usato durante le sessioni di laboratorio per graficare in tempo reale i dati.
E' prevista una funzione per cui i dati inseriti nel plot vengono automaticamente salvati su un file csv in questa cartella.
Il programma inizialmente carica i dati presenti nel file di lettura indicato (se esiste) e li plotta subito.
'''

import scipy
import matplotlib
from matplotlib import pyplot as plt
import os
import csv

matplotlib.rcParams.update({'font.size': 22}) # Cambio dimensione font globale dei plot

file_lettura = "misure2.csv" # Nome del file da cui caricare quelli esistenti.
file_scrittura = "misure2.csv" # Nome del file dove salvare tutti i dati

# Array di storage dove tenere tutti i dati (che poi andranno anche salvati su un csv)
xdata = []
ydata = [] 
sigmadata = []

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Funzione di plotting, la faccio una volta per tutte
def plotxy(x , y):
    #print(type(y[-1]))
    plt.plot(x, y, 'go', linewidth=2, markersize=12, color=[0.01, 0.01, 0.01], alpha=0.5, markerfacecolor=[0, 0, 1])
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("X-Y Plot")
    plt.show()

# Funzione di lettura dati gia presenti
def leggiFile(nomefile):
    myfile = os.path.join(file_lettura)
    if os.path.isfile(myfile):
        with open(file_lettura, 'r') as csvFile:
            reader = csv.reader(csvFile)
            # Per ogni riga del file, aggiungo la tripletta all'array di storage
            for r in reader:
                row = [float(r[0]), float(r[1]), float(r[2])]
                xdata.append(row[0])
                ydata.append(row[1])
                sigmadata.append(row[2])
                print(row)
                if(nomefile != file_scrittura):
                    aggiornaFile(file_scrittura, row) # Funzione che scrive in modalita "append" sul file di scrittura

# Funzione di scrittura su file
def aggiornaFile(nomefile, xysigma):
    with open(nomefile, 'a') as filew:
        writer = csv.writer(filew, delimiter=',')
        writer.writerow(xysigma)




'''
Inizio main routine
'''

leggiFile(file_lettura)
print(xdata, ydata, sigmadata)
plotxy(xdata, ydata)

while(True):
    print("Inserire una tripla (x, y, sigma) nel formato 'dato1 dato2 dato3' per aggiungerla al grafico.")
    nuovo_dato = input()

    if(len(nuovo_dato.split()) == 3):
        x, y, sigma = nuovo_dato.split() # divido i valori in tre variabili diverse
        # Aggiungo i nuovi dati agli array di storage
        xdata.append(float(x))
        ydata.append(float(y))
        sigmadata.append(float(sigma))
        # Scrivo su file
        aggiornaFile(file_scrittura, [x, y, sigma])
        # Grafico
        plotxy(xdata, ydata)
    else:
        print("Input non riconosciuto")