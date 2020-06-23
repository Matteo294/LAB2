''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np

input_file = "data_Faraday.csv"

# Misure e costanti
Gdiff = 0 # da calcolare
Rc = 0
R1 = 0
R2 = 0
R_lim = 0
dR_lim = 0
n_samples = 0
Mrs = 0

frequenze = [10e3, 200e3] # fmax e fmin

Z_eff = []
dZ_eff = []

for i, f in enumerate(frequenze):

    A_in, B_in, A_out, B_out = [[] for n in range(4)]

    for j in range(n_samples):
        segnale = estrazione_segnale(input_file, f, showplots=False)
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
    Z_eff.append(np.imag(C_out / i_s)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
    _dZ_eff = Z_eff * np.sqrt((dC_out/C_out)**2 + (di_s/i_s)**2)
    dZ_eff.append(_dZ_eff * np.imag(Z_eff)/np.absolute(Z_eff))

Ze = linreg(2*pi*frequenze, Z_eff, dZ_eff)
print("Risultati: \t Z0 =", Ze['b'], "\u00B1", Ze['db'], " \t Zeff =", Ze['m'], "\u00B1", Ze['dm'])

fig = bodeplot(frequenze, H=C_out)


