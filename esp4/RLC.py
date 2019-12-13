from labbclass import FDT
from labbclass import Misura
from labbclass import regressione
from labbclass import findWidth
from matplotlib import pyplot as plt 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import os
import sys
import numpy as np 
import csv
import math
import sympy as sym
from sympy import I, re, im

enable_plots = True
enable_simulation_plot = False
enable_data_printing = True
enable_correzione_non_ohmica = False
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
resistenze = [Misura(9939, 2), Misura(5.105e3, 2), Misura(467.34, 0.2)]
L = Misura(2.216e-3, 4e-6)      # risultato della regressione
R_l = Misura(1.606, 0.004)    
#C = Misura(34.61e-9, 0.2e-9)
C = Misura(38.43e-9, 0.0624e-9) # risultato della regressione
#C_l = Misura(1.26585e-10, 7.986e-12)    
C_l = Misura(2.7e-10, 7.986e-12)
C_osc = Misura(128.63e-12, 0.02e-12)    
R_osc = 1e5
C_tot = Misura(C_l.valore + C_osc.valore + C.valore, math.sqrt(C.sigma**2 + C_l.sigma**2 + C_osc.sigma**2))


# grafico bode teorico fatto da simulatore (in realtà non usato)
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
if enable_correzione_non_ohmica:
    d = 0.01
    z_L = jw*L.valore*(1-I*d)
else:
    z_L = jw*L.valore
z_Cl = 1/(jw*C_l.valore)
z_L_tot = 1/(1/(z_L + R_l.valore) + 1/z_Cl)
z_C = 1/(jw*C.valore)
z_Cosc = 1/(jw*C_osc.valore)
z_Rosc = R_osc
z_tot = 1/(1/z_L_tot + 1/z_C + 1/z_Cosc + 1/z_Rosc)

outputf = open("output_dati.csv", 'w')
outputf.write("Resistenza,\"$f_{ris,teo}$[Hz]\",\"$f_r{ris,exp}$[Hz]\",\\makecell{Valore di \\ picco sperimentale},Larghezza picco teorica [Hz],Larghezza picco sperimentale [Hz]\n")

# Ciclo per calcolare funzione di trasferimento, f_ris teorica e sperimentale, f3dB
for R, f, idx in zip(resistenze, resistenze_files, range(len(resistenze))):
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
    # creo la classe e leggo i dati
    rlc = FDT()
    rlc.leggiDati(f)
    # do in pasto alla fdt il numeratore e denominatore trovati di H
    rlc.fdt_teorica(numeratore=num, denominatore=den)
    rlc.f_ris_teo = 1/(2*math.pi*math.sqrt(L.valore*C_tot.valore))
    rlc.sigma_f_ris_teo = 1/(2*math.pi) * math.sqrt(1/(4*(L.valore*C_tot.valore)**3) * ((L.valore*C_tot.sigma)**2 + (C_tot.valore*L.sigma)**2))

    rlc.Q_teorico = R.valore * C_tot.valore * 2*math.pi * rlc.f_ris_teo
    rlc.sigma_Q_teorico = rlc.Q_teorico * 2*math.pi* math.sqrt((R.sigma/R.valore)**2 + (C_tot.sigma/C_tot.valore)**2 + (rlc.sigma_f_ris_teo/rlc.f_ris_teo)**2)
    
    # F RIS SPERIMENTALE
    # regressione lineare delle fasi attorno alla freq di risonanza
    # seleziono solo le fasi < 20
    fase_limite = 20 
    fasi_vicine_zero = np.array([rlc.fase[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    sigma_fasi_vicine_zero = np.array([rlc.sigmaVout[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    freq_fasi_vicine_zero = np.array([rlc.freq[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    sigma_freq_fasi_vicine_zero = np.array([rlc.sigmaT[i] for i in range(np.size(rlc.fase)) if (abs(rlc.fase[i])<=fase_limite and rlc.freq[i]> 1e3)])
    A,B,sigma_A,sigma_B = regressione(fasi_vicine_zero, freq_fasi_vicine_zero, sigma_fasi_vicine_zero, sigma_freq_fasi_vicine_zero)
    rlc.f_ris_regressione = -A/B
    rlc.sigma_f_ris_regressione = rlc.f_ris_regressione * math.sqrt((sigma_A/A)**2 + (sigma_B/B)**2)


    # F3DB PASSA ALTO E PASSA BASSO
    # Le f3db sono a H = max/sqrt(2), trovo i due dati più vicini e interpolo tra loro
    rlc.ampiezza = rlc.Vout/rlc.Vin
    rlc.picco = max(rlc.ampiezza)
    rlc.sigma_picco = 0
    rlc.sigma_ampiezza = rlc.ampiezza * np.sqrt((rlc.sigmaVout/rlc.Vout)**2 + (rlc.sigmaVin/rlc.Vin)**2)
    
    larghezza_exp, sigma_larghezza_exp = findWidth(rlc.ampiezza, rlc.freq,  rlc.f_ris_regressione)
    larghezza_teo, sigma_larghezza_teo = findWidth(10**(rlc._amp_teo/20), rlc._f_teo, rlc.f_ris_teo)

    
    

    # GRAFICI: CERCARE DI CAPIRE IL CODICE A PROPRIO RISCHIO E PERICOLO
    if enable_plots and not enable_correzione_non_ohmica and resistenze.index(R) in grafici_da_plottare:
        # Grafico ampiezza  
        fig = plt.subplot(2,1,1)
        if(resistenze[idx].sigma >= 1):
            plt.suptitle(r'R = %1.f $\Omega$ $\pm$ %1.f $\Omega$' %(resistenze[idx].valore, resistenze[idx].sigma), fontsize=24)
        else:
            plt.suptitle(r'R = %.1f $\Omega$ $\pm$ %.1f $\Omega$' %(resistenze[idx].valore, resistenze[idx].sigma), fontsize=24)
        rlc.plot_teorica_amp()
        plt.errorbar(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), rlc.sigma_amp_dB, rlc.sigmaFreq,'.', ecolor=[0.6350, 0.0780, 0.1840], color=[0.6350, 0.0780, 0.1840], markersize=10,linewidth=1.8,  label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza teorica")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        plt.grid(which='both')
        primary_ax = plt.gca()
        # Zoom ampiezza (differenzio tra alcuni grafici dalle dimensioni diverse)
        if (idx==2):
            zoom = inset_axes(primary_ax, loc='lower right', borderpad=3, width="40%", height="45%")
        else:
            zoom = inset_axes(primary_ax, loc='upper left', borderpad=3, width="40%", height="45%")
        plt.sca(zoom)
        rlc.plot_teorica_amp(axislabel=False)
        plt.errorbar(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), rlc.sigma_amp_dB, rlc.sigmaFreq,'.', ecolor=[0.6350, 0.0780, 0.1840], color=[0.6350, 0.0780, 0.1840], markersize=10, linewidth=1.5, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza teorica")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        if (idx==2):
            zoom.set_xlim(1e4, 3e4)
            zoom.set_ylim(-10, 2)
        else:
            zoom.set_xlim(1.5e4, 2e4)
            zoom.set_ylim(-15, 0.5)
        plt.grid(which='both')
        mark_inset(primary_ax, zoom, loc1=2, loc2=4, fc="none", ec="0.5")

        # Grafico fase
        plt.subplot(2,1,2)
        rlc.plot_teorica_fase()
        plt.errorbar(rlc.freq, rlc.fase, rlc.sigmaFase, rlc.sigmaFreq, '.', ecolor=[0.6350, 0.0780, 0.1840], color=[0.6350, 0.0780, 0.1840], markersize=10, linewidth=1.8, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.5, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        plt.grid(which='both')
        # zoom grafico fase
        primary_ax = plt.gca()
        zoom = inset_axes(primary_ax, loc=3, borderpad=3, width="30%", height="40%")
        plt.sca(zoom)
        rlc.plot_teorica_fase(axislabel=False)
        plt.errorbar(rlc.freq, rlc.fase, rlc.sigmaFase, rlc.sigmaFreq, '.', ecolor=[0.6350, 0.0780, 0.1840], color=[0.6350, 0.0780, 0.1840], markersize=10, linewidth=1.8, label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._fase_teo), max(rlc._fase_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        zoom.set_xlim(16.5e3, 21e3)
        zoom.set_ylim(-25,25)
        mark_inset(primary_ax, zoom, loc1=2, loc2=4, fc="none", ec="0.5")
        plt.grid(which='both')
        plt.show()

    if enable_plots and enable_correzione_non_ohmica and resistenze.index(R) in grafici_da_plottare:
        if(resistenze[idx].sigma >= 1):
            plt.suptitle(r'R = %1.f $\Omega$ $\pm$ %1.f $\Omega$' %(resistenze[idx].valore, resistenze[idx].sigma), fontsize=24)
        else:
            plt.suptitle(r'R = %.1f $\Omega$ $\pm$ %.1f $\Omega$' %(resistenze[idx].valore, resistenze[idx].sigma), fontsize=24)
        rlc.plot_teorica_amp()
        plt.errorbar(rlc.freq, 20*np.log10(rlc.Vout / rlc.Vin), rlc.sigma_amp_dB, rlc.sigmaFreq,'.', ecolor=[0.6350, 0.0780, 0.1840], color=[0.6350, 0.0780, 0.1840], markersize=10,linewidth=1.8,  label="Valori sperimentali")
        plt.plot([rlc.f_ris_teo, rlc.f_ris_teo], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0.4660, 0.6740, 0.1880], label="Frequenza di risonanza teorica")
        plt.plot([rlc.f_ris_regressione, rlc.f_ris_regressione], [min(rlc._amp_teo), max(rlc._amp_teo)], '--', linewidth=1.8, color=[0, 0.4470, 0.7410], label="Frequenza di risonanza per regressione")
        if (idx==2):
            plt.xlim(1.5e4, 2e4)
            plt.ylim(-10, 2)
        else:
            plt.xlim(1.5e4, 2e4)
            plt.ylim(-15, 0.5)
        plt.grid(which='both')
        plt.show()

    if enable_data_printing:
        print("\nResistenza %i \nF ris teorica = %f +- %f" %(idx+1, rlc.f_ris_teo, rlc.sigma_f_ris_teo))
        print("F ris sperimentale = %f +- %f" %(rlc.f_ris_regressione, rlc.sigma_f_ris_regressione))
        # print("F3db passa alto = %f +- %f " %(rlc.f3db_passaalto.valore, rlc.f3db_passaalto.sigma))
        # print("F3db passa basso = %f +- %f" %(rlc.f3db_passabasso.valore, rlc.f3db_passabasso.sigma))
        print("Larghezza picco exp = %f +- %f" %(larghezza_exp, sigma_larghezza_exp))
        print("Larghezza picco teorica = %f +- %f " %(larghezza_teo, sigma_larghezza_teo))
    
    outputf.write("${} \\pm{}$,${}\\pm{}$,${} \\pm {}$,${} \\pm {}$,${}\\pm {}$,${}\\pm {}$\n".format(R.valore, R.sigma, rlc.f_ris_teo, rlc.sigma_f_ris_teo, rlc.f_ris_regressione, rlc.sigma_f_ris_regressione, rlc.picco, rlc.sigma_picco, larghezza_teo, sigma_larghezza_teo, larghezza_exp, sigma_larghezza_exp))

outputf.close()
    





