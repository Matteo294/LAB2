from labbclass import FDT
from labbclass import Misura
from matplotlib import pyplot as plt 
import os
import sys
import numpy as np 
import csv
import sympy as sym
from sympy import I, re, im
from sympy.abc import w

enable_plots = True
enable_simulation_plot = False
# E' possibile indicare da cmd quali grafici visualizzare 
# indicando un numero da 0 a n_resistenze - 1. Mettere x per non plottare
grafici_da_plottare = []
if len(sys.argv) > 1:
    for plot in sys.argv[1:]:
        if plot == 'x':
            enable_plots = False
        else:
            grafici_da_plottare.append(int(plot))

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

resistenze_files = ['dati/RLC/R1.csv', 'dati/RLC/R2.csv', 'dati/RLC/R3.csv'] # un file per ogni resistenza
resistenze = [Misura(9939, 2), Misura(5.10353e3, 2), Misura(467.34, 0.2)] # da cambiare !!
L = Misura(2.1e-3, 0.01) # da cambiare !!
R_l = Misura(0.5, 0.01)     # mettere più cifre significative
C = Misura(48.2e-9, 0.2e-9) # incertezza da cambiare
C_l = Misura(13e-9, 0.01e-9)    # incertezza da cambiare
C_osc = Misura(128.63e-12, 0.02e-12)
R_osc = 1e5
# PLOT SIMULAZIONE
file_lettura = "simulazione.csv"
myfile = os.path.join(file_lettura)
freq = np.array([])
ampl = np.array([])
if os.path.isfile(myfile):
    with open(file_lettura, 'r') as csvFile:
        reader = csv.reader(csvFile)
        # Per ogni riga del file, leggo i valori frequenza e ampiezza
        for r in reader:
            row = [float(r[0]), float(r[1])] 
            freq = np.append(freq, row[0])
            ampl = np.append(ampl, row[1])
else:
    print("Problema: non trovo il file " + file_lettura)
if enable_simulation_plot:
    plt.plot(freq, ampl, '-')
    plt.show()


# calcolo l'impedenza totale formata dal parallelo di L, C, Osc
z_L = (I*w*L.valore + R_l.valore) / (1 - w**2 * L.valore * C_l.valore + I*w*R_l.valore*C_l.valore)
z_C = 1/(I*w*C.valore)
z_Cosc = 1/(I*w*C_osc.valore)
z_Rosc = R_osc
z_tot = 1/(1/z_L + 1/z_C + 1/z_Cosc + 1/z_Rosc)

# Creo l'oggetto dalla classe FDT e leggo i dati dal file nel formato: freq Vout fase Vin
for R, f in zip(resistenze, resistenze_files):
    
    # calcolo la funzione di trasferimento, la separo in numeratore e denominatore e li rendo polinomi
    num, den = sym.fraction(sym.cancel(H))
    num = sym.Poly(num, w)
    den = sym.Poly(den, w)
    # trovo i coefficienti del numeratore e del denominatore
    num = num.all_coeffs()
    num = [complex(num[i]) for i in range(len(num))]    # altrimenti sono complessi di "sympy" e non funzionano in fdt
    den = den.all_coeffs()
    den = [complex(den[i]) for i in range(len(den))]
    
    rlc = FDT()
    rlc.leggiDati(f)
    # do in pasto alla fdt il numeratore e denominatore trovati di H
    rlc.fdt_teorica(numeratore=num, denominatore=den) 

    rlc.f_ris = 1 / (2*np.pi * np.sqrt(L.valore*C.valore))
    
    
    if enable_plots and resistenze.index(R) in grafici_da_plottare:
        rlc.plot_teorica_ampiezza()

        plt.semilogx(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), '.', markersize=10)
        plt.plot([rlc.f_ris, rlc.f_ris], [min(rlc._ampiezza_teo), max(rlc._ampiezza_teo)], '--', linewidth=1.8, color='red', label="Frequenza di risonanza")
        plt.legend()
        plt.grid()
        plt.show()

        rlc.plot_teorica_fase()
        plt.semilogx(rlc.freq, rlc.fase, '.', markersize=10)
        plt.plot([rlc.f_ris, rlc.f_ris], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color='red', label="Frequenza di risonanza")
        plt.legend()
        plt.grid()
        plt.show()