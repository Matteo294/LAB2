'''
Classe esperimenti di laboratorio. Questa classe contiene tutte le possibili funzioni e commandi che possono tornare utili nell'analisi dati di laboratorio.
Ne rimane escluso il programma per plottare in live i dati che è stato mantenuto separato per comodità di utilizzo.
Le linee guida di utilizzo sono scritte nel file .txt in questa cartella (devo ancora scriverlo LoL)
C'è una superclasse (Analisi) in cui sono scritte tutte le funzioni comuni, ovvero quelle funzioni che vanno bene per ogni tipo di analisi.
Ci sono poi delle sottoclassi (per ora solo LinearFit) con delle funzioni particolari per l'analisi cercata (utili per il fit lineare in questo caso).
''' 

'''
Lista di cose da aggiungere:
chi quadro a 2 parametri, trasferimento incertezze, calcolo probabilità di eccedere il chi quadro, migliorare la funzione di print, 
funzione per i grafici (deve creare l'oggetto self.axes in modo da potervi accedere dall'esterno per modificare il grafico), grafico dei residui
'''
import numpy as np
import scipy
import math
import matplotlib
from numpy.linalg import inv
from numpy import matrix as mat
import os
import csv

class Analisi:
    
    # Superclasse: usata per le cose comuni a tutte le sottoclassi
    
    def __init__(self):
        self.xdata = np.asarray([])
        self.ydata = np.asarray([])
        self.sigmax = np.asarray([])
        self.sigmay = np.asarray([])

    # Lo swap permette di invertre X e Y, mentre lo scale permette di cambiare multiplo di udm ai dati
    def leggiDati(self, file_lettura, swap_xy=0, scale_x=1, scale_y=1):
        myfile = os.path.join(file_lettura)
        if os.path.isfile(myfile):
            with open(file_lettura, 'r') as csvFile:
                reader = csv.reader(csvFile)
                # Per ogni riga del file, aggiungo la tripletta all'array di storage
                for r in reader:
                    row = [float(r[swap_xy]), float(r[not swap_xy])] # Leggo le due colonne r[0] e r[1]. Se swap_xy = 1 scambio le due colonne (leggo prima r[1] e poi r[0])
                    self.xdata = np.append(self.xdata, row[0] * scale_x)
                    self.ydata = np.append(self.ydata, row[1] * scale_y)

    def add_sigmas(self, sigmax=0, sigmay=0):
        # Se la sigma passata per parametro è costante (int o float), creo un array di costanti, altrimenti copio l'array
        if isinstance(sigmax, (int, float)):
            self.sigmax = np.full(self.xdata.size, sigmax)
        else:
            self.sigmax = np.asarray(sigmax)
        if isinstance(sigmay, (int, float)):
            self.sigmay = np.full(self.ydata.size, sigmay)
        else:
            self.sigmay = np.asarray(sigmay)

    def __str__(self):
        # Questa funzione dev'essere implementata nelle sottoclassi. Se non lo si fa, lancio un errore
        raise NotImplementedError




    
class LinearFit(Analisi):

    # Ho copiato solo una parte dalla funzione, non quella in cui vengono trasferite le incertezze
    def reg_lin(self, trasferisci=False): # a + bx model
        # Guardo se esiste la variabile sigma_regressione. All'inizio non esisterà, quindi la setto uguale a sigmay e sigmax non influisce.
        # Se poi voglio trasferire sigmax, allora aggiorno sigma_regressione e poi rieseguo il codice
        # In questo modo posso reinvocare dopo lo stesso codice cambiando solo sigma_regressione
        try:
            sigma_regressione
        except NameError:
            sigma_regressione = self.sigmay

        w = 1/sigma_regressione
        delta = sum(w)*sum(self.xdata**2/w) - (sum(self.xdata/w))**2
        self.A = 1/delta * (sum(self.xdata**2/w)*sum(self.ydata/w) - sum(self.xdata/w)*sum(self.xdata*self.ydata/w))
        self.B = 1/delta * (sum(w)*sum(self.xdata*self.ydata/w) - sum(self.xdata/w)*sum(self.ydata/w))
        self.sigma_A = math.sqrt(1/delta * sum(self.xdata**2/w))
        self.sigma_B = math.sqrt(1/delta * sum(w))
        
        if (trasferisci==True):
            sigma_trasformata = abs(self.B)*self.sigmax
            sigma_regressione = np.sqrt(self.sigmay**2 + sigma_trasformata**2)
            self.reg_lin(trasferisci=False)


    def __str__(self):
        # Printo con 4 decimali i valori e con 1 il chi ridotto (troncamento, non approssimazione)
        # Il codice \u00B1 è per il +-
        return(u"\nIntercetta: {0:.4f} \u00B1 {1:.4f} \t Coefficiente angolare: {2:.4f} \u00B1 {3:.4f} \t Chi quadro ridotto: {4:.1f} ({5})\n".format(
            self.A, self.sigma_A, self.B, self.sigma_B, self.chi_ridotto, self.xdata.size))
    
    def chi_quadro(self, n_params=1):
        if n_params == 1:
            self.chi_q = sum((self.B*self.xdata - self.ydata)**2 / (self.sigmay)**2)
            self.chi_ridotto = self.chi_q / self.xdata.size