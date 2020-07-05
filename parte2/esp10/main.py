''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np
import sys
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *

plot_flag = 0
print_flag = 0
if len(sys.argv) > 1:
    plot_flag = int(sys.argv[1])

fileimmagine = '../esp9/Immagini/Faraday.png'

def Gdiff(f):
    w = 2*pi*f
    # carico = oscilloscopio
    Cosc = 110e-12
    Rosc = 1e6
    Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
    # parametri che servono
    Rc = 9.830e3
    Re = 119.25
    re = 45
    # Gdiff
    G0_diff = Rc / (2*(re + Re))
    return abs((G0_diff*Zosc)/(Rc + Zosc))

def dGdiff(f):
    return Gdiff(f) / 100

# Misure e costanti
R_lim = 11
dR_lim = 0.001
n_samples = 5 # numero misure ripetute

# Parametri bobine
mu0 = 4*pi*1e-7
sigma1 = np.pi * (17.5e-3/2)**2
sigma2 = sigma1
n_S = 30
n_R = 28
M_dipolo = lambda d: 2 * mu0/(4*np.pi) * n_S * n_R * sigma1 * sigma2 / d**3

distanze = [1.35e-2, 2.3e-2, 4.6e-2, 10.5e-2, 4.4e-2, 1.8e-2] # Sistemare prima distanza
frequenze = [1e3, 50e3, 150e3]
omegas = [2*np.pi*f for f in frequenze]

''' Analisi accoppiamento '''
Z_ctrl = []
dZ_ctrl = []
for i, f in enumerate(frequenze):

    Z = []
    
    for j in range(n_samples):

        filename = "Data/h/f" + str(i+1) + "/" + str(j+1) + ".csv"
        s_in, s_out, ds_in, ds_out = estrazione_segnale(filename, f, showplots=False)
        [_, A_in, B_in] = s_in 
        [_, A_out, B_out] = s_out 
        [_, dA_in, dB_int] = ds_in
        [_, dA_out, dB_out] = ds_out

        # Calcolo ampiezza complessa del segnale in ingresso
        C_in = A_in - 1j*B_in

        # Calcolo ampiezza complessa segnale in uscita
        C_out = A_out - 1j*B_out

        # Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)
        Z.append(np.imag(C_out / C_in / Gdiff(f) * R_lim)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
    
    Z_ctrl.append(np.mean(Z))
    dZ_ctrl.append(np.std(Z, ddof=1)) # Non solo questa, anche 1% fondoscala

params = linreg(omegas, Z_ctrl, dZ_ctrl)
M_ctrl = params['m']
'''-----------------------------------'''

# Array induttanza mutua a diverse distanze
Mrs = [] 
dMrs = []

for d in range(len(distanze)):

    base_input_file = "Data/d" + str(d+1)
    Z_eff = []
    dZ_eff = []

    for i, f in enumerate(frequenze):

        Z = []
        for j in range(n_samples):

            filename = base_input_file + "/f" + str(i+1) + "/" + str(j+1) + '.csv'
            s_in, s_out, ds_in, ds_out = estrazione_segnale(filename, f, showplots=False)
            [_, A_in, B_in] = s_in 
            [_, A_out, B_out] = s_out 
            [_, dA_in, dB_int] = ds_in
            [_, dA_out, dB_out] = ds_out

            # Calcolo ampiezza complessa del segnale in ingresso
            C_in = A_in - 1j*B_in

            # Calcolo ampiezza complessa segnale in uscita
            C_out = A_out - 1j*B_out

            # Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)
            Z.append(np.imag(C_out / C_in / Gdiff(f) * R_lim)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
    
        Z_eff.append(np.mean(Z))
        dZ_eff.append(np.std(Z, ddof=1)) # Non solo questa, anche 1% fondoscala

    Ze = linreg(omegas, Z_eff, dZ_eff)
    if print_flag:
        print("d =", distanze[d], ": \t Z0 =", Ze['b'], "\u00B1", Ze['db'], " \t Zeff =", Ze['m'], "\u00B1", Ze['dm'], '\n')
    if plot_flag == 1:
        fig = bodeplot(frequenze, Amp=np.abs(Z_eff), Phase=np.angle(Z_eff))
        plt.show()

    Mrs.append(Ze['m'] - M_ctrl)
    dMrs.append(Ze['dm'])

# Leggo i valori dei modelli teorici
d, val, approx = readCSV('Mrs/induzione.csv')

# Cambio udm ai dati
d = [dist/1e-3 for dist in d]
distanze = [dist/1e-3 for dist in distanze]
approx = [a/1e-6 for a in approx]
val = [v/1e-6 for v in val]
Mrs = [m/1e-6 for m in Mrs]
dMrs = [d/1e-6 for d in dMrs]

# Plot dati, salvo nella cartella immagini
plt.semilogy(d, approx, label="Approssimazione dipolo", linewidth=1.8, color=[0, 1, 0])
plt.semilogy(d, val, label="Formula Neumann", linewidth=1.8, color='blue')
plt.errorbar(distanze, Mrs, yerr=np.real(dMrs), fmt='.', markersize=8, markerfacecolor='red', color='black', label="Dati sperimentali")
plt.xlabel(r'Distanza   [mm]')
plt.ylabel(r'$M_{RS}   [\mu H]$')
plt.ylim((5e-4, 1e1))
plt.title("Mutua induzione tra bobine", fontsize=20)
plt.legend()
plt.savefig(fileimmagine, bbox_inches='tight', dpi=1000)
plt.show()
