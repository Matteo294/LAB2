import numpy as np
import scipy
import math
from matplotlib import pyplot as plt
import os 
import csv

# nome del file da cui si prenderanno i dati per la regressione
nomefile = "monte10V_5mA.csv"

# inizializzo i miei array 
xdata = np.array([])
ydata = np.array([])
sigma_x = np.array([])
sigma_y = np.array([])


# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
myfile = os.path.join(nomefile)

# riempio i miei array coi dati dal file
if os.path.isfile(myfile):
    with open(nomefile, 'r') as csvFile:
        reader = csv.reader(csvFile)
        # Per ogni riga del file, aggiungo la tripletta all'array di storage
        for row in reader:
            xdata = np.append(xdata, row[0])
            ydata = np.append(ydata, row[1])
            sigma_x = np.append(sigma_x, row[2])
            sigma_y = np.append(sigma_y, row[3])
else :
    print("File di input non trovato")

# faccio un cast a float, altrimenti non posso fare operazioni sugli array
xdata = xdata.astype('float64')
ydata = ydata.astype('float64')
sigma_x = sigma_x.astype('float64')
sigma_y = sigma_y.astype('float64')


## E qui inizia la regressione
# ipotizzando che ydata = A*xdata + B, trovo A e B e le loro sigme
n_misure = xdata.size

# faccio una funzione che mi fa i calcolazzi sporchi
def regredisci(xdata, ydata, sigma):
    w = 1/sigma**2
    delta = sum(w)*sum(xdata**2/w) - (sum(xdata/w))**2
    A = 1/delta * (sum(xdata**2/w)*sum(ydata/w) - sum(xdata/w)*sum(xdata*ydata/w))
    B = 1/delta * (sum(w)*sum(xdata*ydata/w) - sum(xdata/w)*sum(ydata/w))
    sigma_A = math.sqrt(1/delta * sum(xdata**2/w))
    sigma_B = math.sqrt(1/delta * sum(w))
    return A, B, sigma_A, sigma_B

# prima la invoco una volta con sigma = sigma_y (tanto va fatta comunque), poi chiedo se contare anche il contributo di sigma_x
A, B, sigma_A, sigma_B = regredisci(xdata, ydata, sigma_y)
flag = bool(input("Vuoi contare anche il contributo di xdata? (0=no, 1=si) "))
if (flag):
    sigma_trasformata = abs(B)*sigma_x
    sigma_regressione = np.sqrt(sigma_y**2 + sigma_trasformata**2)
    A, B, sigma_A, sigma_B = regredisci(xdata, ydata, sigma_regressione)
print("A = ", A, "\nsigma_A = ", sigma_A, "\nB = ", B, "\nsigma_B = ", sigma_B)