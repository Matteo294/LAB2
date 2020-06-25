''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np
import sys

plot_flag = 0
if len(sys.argv) > 1:
    plot_flag = int(sys.argv[1])

def Gdiff(f):
    ReG = -3.0955117 + 5.03032336*f + 1.1460822e-9*f**2 - 4.1463407e-15*f**3
    ImG = -3.14127842 + 3.16888061*f -2.01140302e-9*f**2 - 3.83960293e-15*f**3
    return ReG + 1j*ImG

# Misure e costanti
R_lim = 4.634
dR_lim = 0.001
n_samples = 5 # numero misure ripetute

# Parametri bobine
mu0 = 4*pi*1e-7
sigma1 = np.pi * (17.5e-3/2)**2
sigma2 = sigma1
NS = 0
NR = 0
M_dipolo = lambda d: 2 * mu0/(4*np.pi) * NS * NR * sigma1 * sigma2 / d**3

distanze = [0, 2.3e-2, 4.6e-2, 10.5e-2, 4.4e-2, 1.8e-2]
frequenze = [1e3, 50e3, 150e3]
dG_diff = [0, 0, 0] # Incertezza per i guadagni sopra

# Array induttanza mutua a diverse distanze
Mrs = [] 
dMrs = []

for d in distanze:

    base_input_file = "Data/d" + str(d+1)
    Z_eff = []
    dZ_eff = []

    for i, f in enumerate(frequenze):

        A_in, B_in, A_out, B_out = [[] for n in range(4)]

        for j in range(n_samples):
            filename = base_input_file + "/f" + str(i+1) + "/" + str(j+1) + '.csv'
            segnale = estrazione_segnale(filename, f, showplots=False)
            A_in.append(segnale["A_in"])
            B_in.append(segnale["B_in"])
            A_out.append(segnale["A_out"])
            B_out.append(segnale["B_out"])
        
        # Valori medi delle misure ripetute
        A_in_media = np.average(A_in)
        A_in_std = np.std(A_in, ddof=1)

        B_in_media = np.average(B_in)
        B_in_std = np.std(B_in, ddof=1)

        A_out_media = np.average(A_out)
        A_out_std = np.std(A_out, ddof=1)

        B_out_media = np.average(B_out)
        B_out_std = np.std(B_out, ddof=1)

        # Calcolo ampiezza complessa del segnale in ingresso
        C_in = A_in_media - 1j*B_in_media
        dC_in = np.sqrt(A_in_std**2 + B_in_std**2)

        # Calcolo ampiezza complessa segnale in uscita
        C_out = A_out_media - 1j*B_out_media
        dC_out = np.sqrt(A_out_std**2 + B_out_std**2)

        # Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)
        i_s = C_in / R_lim
        di_s = i_s * np.sqrt((dC_in/C_in)**2 + (dR_lim/R_lim)**2) 
        Z_eff.append(np.imag(C_out / Gdiff(f[i]) / i_s)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
        _dZ_eff = Z_eff * np.sqrt((dC_out/C_out)**2 + (dG_diff[i]/Gdiff(f[i]))**2 + (di_s/i_s)**2)
        dZ_eff.append(_dZ_eff * np.imag(Z_eff)/np.absolute(Z_eff))

    Ze = linreg(2*pi*frequenze, Z_eff, dZ_eff)
    print("Risultati d =", d, ": \t Z0 =", Ze['b'], "\u00B1", Ze['db'], " \t Zeff =", Ze['m'], "\u00B1", Ze['dm'])

    if plot_flag == 1:
        fig = bodeplot(frequenze, H=C_out)
        fig.show()

    Mrs.append(Ze['m'])
    dMrs.append(Ze['dm'])

plt.errorbar(distanze, Mrs, dMrs, '.')
plt.show()
