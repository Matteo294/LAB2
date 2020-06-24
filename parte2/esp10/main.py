''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np
import sys

plot_flag = 0
if len(sys.argv) > 1:
    plot_flag = int(sys.argv[1])

# Misure e costanti
Gdiff = 0 # da calcolare
Rc = 0
R1 = 0
R2 = 0
R_lim = 0
dR_lim = 0
n_samples = 0
Mrs = 0
# Parametri bobine
mu0 = 4*pi*1e-7
sigma1 = np.pi * (17.5e-3/2)**2
sigma2 = sigma1
NS = 0
NR = 0
M_dipolo = lambda d: 2 * mu0/(4*np.pi) * NS * NR * sigma1 * sigma2 / d**3

distanze = [0, 2.3, 4.6, 10.5, 4.4, 1.8] * 1e-2 # m
frequenze = [1, 50, 150] * 1e3
G_diff = [0, 0, 0] # Guadagno differenziale per frequenze scelte
dG_diff = [0, 0, 0] # Incertezza per i guadagni sopra

# Array induttanza mutua a diverse distanze
Mrs = [] 
dMrs = []

for d in distanze:

    base_input_file = "Faraday" + int(d)

    Z_eff = []
    dZ_eff = []

    for i, f in enumerate(frequenze):

        A_in, B_in, A_out, B_out = [[] for n in range(4)]

        for j in range(n_samples):
            filename = base_input_file + str(i+1) + '.csv'
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
        Z_eff.append(np.imag(C_out / G_diff[i] / i_s)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
        _dZ_eff = Z_eff * np.sqrt((dC_out/C_out)**2 + (dG_diff[i]/G_diff[i])**2 + (di_s/i_s)**2)
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
