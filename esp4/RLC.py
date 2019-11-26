from labbclass import FDT
from labbclass import Misura
from labbclass import regressione
from matplotlib import pyplot as plt 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
import sys
import numpy as np 
import csv
import math
import sympy as sym
from sympy import I, re, im

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
L = Misura(2.09e-3, 0.01) # da cambiare !!
R_l = Misura(0.5, 0.01)     # mettere più cifre significative
C = Misura(34.61e-9, 0.2e-9) # incertezza da cambiare
C_l = Misura(1.335e-10, 0.01e-9)    # incertezza da cambiare
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
jw = sym.Symbol('jw')   # poi mi serviranno i coefficienti di jw, la metto come variabile
z_L = jw*L.valore
z_Cl = 1/(jw*C_l.valore)
z_L_tot = 1/(1/(z_L + R_l.valore) + 1/z_Cl)
z_C = 1/(jw*C.valore)
z_Cosc = 1/(jw*C_osc.valore)
z_Rosc = R_osc
z_tot = 1/(1/z_L_tot + 1/z_C + 1/z_Cosc + 1/z_Rosc)

# Creo l'oggetto dalla classe FDT e leggo i dati dal file nel formato: freq Vout fase Vin
for R, f in zip(resistenze, resistenze_files):
    # calcolo la funzione di trasferimento, la separo in numeratore e denominatore e li rendo polinomi
    H = z_tot / (R.valore + z_tot)
    num, den = sym.fraction(sym.cancel(H))
    num = sym.Poly(num, jw)
    den = sym.Poly(den, jw)
    # trovo i coefficienti del numeratore e del denominatore
    num = num.all_coeffs()
    num = [complex(num[i]) for i in range(len(num))]    # altrimenti sono complessi di "sympy" e non funzionano in fdt
    den = den.all_coeffs()
    den = [complex(den[i]) for i in range(len(den))]
    
    rlc = FDT()
    rlc.leggiDati(f)
    # do in pasto alla fdt il numeratore e denominatore trovati di H
    rlc.fdt_teorica(numeratore=num, denominatore=den)
    rlc.f_ris_teo = rlc._f_teo[np.where(rlc._ampiezza_teo == max(rlc._ampiezza_teo))]
    
    # regressione lineare delle fasi attorno alla freq di risonanza
    # seleziono solo le fasi < 20
    fase_limite = 20 
    fasi_vicine_zero = np.array([rlc.fase[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    sigma_fasi_vicine_zero = np.array([rlc.sigmaVout[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    freq_fasi_vicine_zero = np.array([rlc.freq[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    sigma_freq_fasi_vicine_zero = np.array([rlc.sigmaT[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    A,B,sigma_A,sigma_B = regressione(fasi_vicine_zero, freq_fasi_vicine_zero, sigma_fasi_vicine_zero, sigma_freq_fasi_vicine_zero)
    
    rlc.f_ris_regressione = -A/B
    print(rlc.f_ris_regressione)
    
    if enable_plots and resistenze.index(R) in grafici_da_plottare:
        plt.subplot(1,2,1)
        rlc.plot_teorica_ampiezza()
        plt.semilogx(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), '.', color=[0.6350, 0.0780, 0.1840], markersize=10, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._ampiezza_teo), max(rlc._ampiezza_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza teorica")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._ampiezza_teo), max(rlc._ampiezza_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        
        plt.legend()
        plt.grid(which='both')


        plt.subplot(1,2,2)
        rlc.plot_teorica_fase()
        plt.plot(rlc.freq, rlc.fase, '.', color=[0.6350, 0.0780, 0.1840], markersize=10, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        plt.xlim(1, 1e8)
        plt.legend()
        plt.grid(which='both')

        zoom = inset_axes(plt.gca(), loc=1, width="30%", height="40%", borderpad=1)
        plt.sca(zoom)
        rlc.plot_teorica_fase(axislabel=False)
        plt.plot(rlc.freq, rlc.fase, '.', color=[0.6350, 0.0780, 0.1840], markersize=10, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        zoom.set_xlim(16.5e3, 21e3)
        zoom.set_ylim(-20,20)
        plt.grid(which='both')
        plt.show()