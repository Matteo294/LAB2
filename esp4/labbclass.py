'''
Classe esperimenti di laboratorio. Questa classe contiene le principali funzioni e commandi che possono tornare utili nell'analisi dati di laboratorio.
Ne rimane escluso il programma per plottare in live i dati che è stato mantenuto separato per comodità di utilizzo.
C'è una superclasse (Analisi) in cui sono scritte tutte le funzioni comuni, ovvero quelle funzioni che vanno bene per ogni tipo di analisi come
leggere i dati da un file, aggiungere le incertezze ai dati ecc.
Ci sono poi delle sottoclassi (per ora solo LinearFit) con delle funzioni particolari per l'analisi cercata (utili per il fit lineare in questo caso).
I parametri xscale e yscale che si trovano qua e la nelle funzioni servono per cambiare udm ai dati. In particolare si trovano nella funzione leggiDati
(i due parametri servono per memorizzare i dati negli array con un multiplo di udm) e di nuovo nelle funzioni di plot (per facilitare la visualizzazione
dei dati potrebbe essere comodo ri-cambiare udm)
''' 

import numpy as np
from scipy import stats
from scipy import signal
import math
from matplotlib import pyplot as plt
#from numpy.linalg import inv
#from numpy import matrix as mat
import os
import csv

class Analisi:
    
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
        else:
            print("Problema: non trovo il file " + file_lettura)

    def add_sigmas(self, sigmax=0, sigmay=0):
        # Se il parametro passato è una costante, creo un array di costanti. Altrimenti copio l'array.
        if isinstance(sigmax, (int, float)):
            self.sigmax = np.full(self.xdata.size, sigmax)
        else:
            self.sigmax = np.asarray(sigmax)
        if isinstance(sigmay, (int, float)):
            self.sigmay = np.full(self.ydata.size, sigmay)
        else:
            self.sigmay = np.asarray(sigmay)

    def add_resistenza_tester_ICE(self, input_amp, input_vol):
        ''' input_amp: 1=5A; 2=500mA; 3=50mA; 4=5mA; 5=500uA; 6=50uA
            input_vol: 1= 100mV; 2=2V ; 3=10V ; 4=50V ; 5=200V ; 6=500V ; 7=1000V '''
        # I commenti messi con le ''' subito sotto l'intestazione della funzione servono da documentazione su come funziona la funzione.
        # Appaiono come suggerimenti quando le usi su vs code, molto utili per sapere cosa mettere come parametri
        R_interne = [0.064, 0.576, 5.76, 57.6, 576, 5760]
        R_tot = sum(R_interne)
        R_a = 1600
        R_sh = 0
        # calcolo R shunt
        for i in range(input_amp):
            R_sh += R_interne[i]
        self.R_amp = R_sh*(R_tot - R_sh + R_a) / (R_tot + R_a)
        # il terminale 6 ha una resistenza in piu
        if input_amp == 6:
            self.R_amp += 720

        # calcolo resistenza del voltmetro
        R_parallelo_int = R_tot*R_a / (R_tot + R_a)         # parallelo interno completo R_tot || R_a
        R_terminali_vol = [32.2e3, 6.52e3, 160e3, 800e3, 6e6, 100e6]
        if input_vol == 1:
            self.R_vol = 720 + R_parallelo_int
        else:
            R_entrata = 0
            for i in range(input_vol):
                R_entrata += R_terminali_vol[i]
            self.R_vol = R_parallelo_int + R_entrata

    # Queste funzioni devono essere implementate nelle sottoclassi. Se non lo si fa, lancio un errore
    def __str__(self):
        raise NotImplementedError

    def plotData(self):
        raise NotImplementedError




    
class LinearFit(Analisi):
    # Uso l'incertezza sigma_reg per la regressione e nel chi quadro. 
    # Di default è solo sigmay, se dovrò trasferire sigmax, lo farò in reg_lin
    def add_sigmas(self, sigmax=0, sigmay=0):
        super(LinearFit, self).add_sigmas(sigmax, sigmay)
        self.sigma_reg = self.sigmay

    # Fit lineare con una funzione A + Bx
    def reg_lin(self, trasferisci=False):
        w = 1/self.sigma_reg
        delta = sum(w)*sum(self.xdata**2/w) - (sum(self.xdata/w))**2
        self.A = 1/delta * (sum(self.xdata**2/w)*sum(self.ydata/w) - sum(self.xdata/w)*sum(self.xdata*self.ydata/w))
        self.B = 1/delta * (sum(w)*sum(self.xdata*self.ydata/w) - sum(self.xdata/w)*sum(self.ydata/w))
        self.sigma_A = math.sqrt(sum(self.xdata**2 * w) / delta)
        self.sigma_B = math.sqrt(sum(w) / delta)
        
        if (trasferisci==True):
            sigma_trasformata = abs(self.B)*self.sigmax
            self.sigma_reg = np.sqrt(self.sigmay**2 + sigma_trasformata**2)
            self.reg_lin(trasferisci=False)

    def __str__(self):
        # Controllo se esistono le variabili e man mano le aggiungo alla frase di print
        frase = ""
        try:
            frase += "Intercetta: {0:.4f} \u00B1 {1:.4f} \t".format(self.A, self.sigma_A)
        except:
            pass
        try:
            frase += "Coefficiente angolare: {0:4f} \u00B1 {1:.4f} \t".format(self.B, self.sigma_B)
        except:
            pass
        try :
            frase += "Chi quadro ridotto: {0:.4f} pvalue={1:.4f} su {2} dati \t".format(self.chi_ridotto, self.probabilita_chi, self.xdata.size)
        except:
            pass
        return frase
    
    def chi_quadro(self, n_params=2):
        if n_params == 1:
            self.chi_q = sum((self.B*self.xdata - self.ydata)**2 / (self.sigma_reg)**2)
            self.chi_ridotto = self.chi_q / (self.xdata.size - 1)
        elif n_params == 2:
            self.chi_q = sum((self.B*self.xdata + self.A - self.ydata)**2 / (self.sigma_reg)**2)
            self.chi_ridotto = self.chi_q / (self.xdata.size - 2)
            self.probabilita_chi = stats.chi2.sf(self.chi_q, self.xdata.size - 2)


    def residui(self, n_params=2, xlabel='ID Misura', ylabel='Y', title=None, xscale=1, yscale=1):
        if n_params == 1:
            residui = self.B*self.xdata - self.ydata
        elif n_params == 2:
            residui = self.ydata - (self.B*self.xdata + self.A)

        params = plt.errorbar(self.xdata*xscale, residui, self.sigma_reg, np.zeros(self.xdata.size), 'o', ecolor='red')

        self.residui_plot = params[0]

        # Sistemo alcuni parametri
        self.residui_plot.set_color('black')
        self.residui_plot.set_alpha(0.8)
        self.residui_plot.set_markersize(8)
        
        plt.xlabel(xlabel, fontsize=18)
        plt.ylabel(ylabel, fontsize=18)

        if title is not None:
            plt.title(title, fontsize=24)
        
        plt.plot([0, max(self.xdata*xscale)], [0, 0], '--', color='gray', linewidth=1.8)

        plt.grid()

    def plotData(self, xlabel='X data', ylabel='Y', title=None, xscale=1, yscale=1):
        params = plt.errorbar(xscale*self.xdata, yscale*self.ydata, self.sigmay*yscale, self.sigmax*xscale, 'o', ecolor='red', linewidth=1.8, markersize=8, label='Dati misurati')
        
        # Memorizzo l'oggetto del grafico di modo da potervi accedere dall'esterno
        self.data_plot = params[0]
        self.ax_data = plt.gca()

        # Sistemo alcuni parametri
        self.data_plot.set_color('black')
        self.data_plot.set_alpha(0.8)

        plt.xlabel(xlabel, fontsize=18)
        plt.ylabel(ylabel, fontsize=18)

        if title is not None:
            plt.title(title, fontsize=24)   

        # Array di 1000 punti nel range xmin-xmax dove calcolare la funzione di regressione
        xrange = np.linspace(min(self.xdata), max(self.xdata), 1000)
        self.regression_plot, = plt.plot(xrange * xscale, (self.A + self.B*xrange) * yscale, label='Regressione lineare')
        
        plt.grid() # Griglia

# Classe per l'analisi di funzioni di trasferimento
# In questo caso con x si intende Vin e con y Vout
class FDT(Analisi):

    def __init__(self):
        self.Vout = np.array([])
        self.freq = np.array([])
        self.fase = np.array([])
        self.Vin = np.array([])

    def leggiDati(self, file_lettura, scale_f=1, scale_Vout=1, scale_Vin=1):
        myfile = os.path.join(file_lettura)
        if os.path.isfile(myfile):
            with open(file_lettura, 'r') as csvFile:
                reader = csv.reader(csvFile)
                # Per ogni riga del file, aggiungo all'array di storage.
                # Le colonne sono: freq, Vout, fase(°), Vin
                for r in reader:
                    row = [float(r[0]), float(r[1]), float(r[2]), float(r[3])]
                    self.freq = np.append(self.freq, row[0] * scale_f)
                    self.Vout = np.append(self.Vout, row[1] * scale_Vout)
                    self.fase = np.append(self.fase, row[2])
                    self.Vin = np.append(self.Vin, row[3] * scale_Vin)
        else:
            print("Problema: non trovo il file " + file_lettura)
    
    # Eventualmente è possibile fissare "manualmente" la Vin (comoda se è Vin=cost.)
    def setVin(self, value):
        if isinstance(value, (int, float)):
            self.Vin = np.full(self.freq.size, value)
        else:
            self.Vin = np.asarray(value)

    def plotFDT(self):
        self.dataplot_magnitude, = plt.semilogx(self.freq, self.ydata/self.xdata)

    # numeratore: coefficienti del polinomio a numeratore della fdt potenza decrescente.
    # denominatore: coefficienti del polinomio a denominatore della fdt potenza decrescente
    # Esemio: H(w) = (t1*(jw)^2 + t2*(jw) + 5) / (3*(jw) + 4) --> 
    # --> numeratore = [t1, t2, 5] e denominatore = [3, 4]
    # Se auto_f = True, il range di f viene basato sui dati misurati
    def fdt_teorica(self, numeratore=1, denominatore=1):
        fdt = signal.lti(numeratore, denominatore)
        _w_plot, self._ampiezza_teo, _fase_teo = signal.bode(fdt, n = 1000) 
        self._f_teo = _w_plot / (2 * np.pi)       
        self._fase_teo = _fase_teo # rad -> deg
    
    # I parametri f e ampiezza servono per plottare un eventuale fdt creata all'esterno
    def plot_teorica_ampiezza(self, f=None, ampiezza=None):
        if f is None:
            f = self._f_teo
        if ampiezza is None:
            ampiezza = self._ampiezza_teo

        self.teoplot_ampiezza, = plt.semilogx(f, ampiezza, linewidth=1.8, color=[0, 0, 0])
        plt.xlabel("Frequenza [Hz]", fontsize=18)
        plt.ylabel("Ampiezza [dB]", fontsize=18)
        plt.title("Diagramma dell'ampiezza della F.D.T.", fontsize=24)

    def plot_teorica_fase(self, f=None, fase=None):
        if f is None:
            f = self._f_teo
        if fase is None:
            fase = self._fase_teo

        self.teoplot_fase, = plt.semilogx(f, fase, linewidth=1.8, color=[0, 0, 0])
        plt.xlabel("Frequenza [Hz]", fontsize=18)
        plt.ylabel("Fase [°]", fontsize=18)
        plt.title("Diagramma della fase della F.D.T.", fontsize=24)



    def __str__(self):
        return "Print non ancora implementata"


# Classe molto semplice per memorizzare i dati nel formato misura-incertezza
class Misura:
    def __init__(self, value, sigma):
        if isinstance(value, (int, float)):
            self.valore = value
            self.sigma = sigma
        else:
            self.valore = np.array(value)
            if isinstance(sigma, (int, float)):
                self.sigma = np.full(self.valore.size, sigma)
            else:
                self.sigma = np.array(sigma)