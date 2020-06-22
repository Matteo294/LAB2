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
nsamples = 0

frequenze = [10e3, 200e3] # fmax e fmin

for i, f in enumerate(frequenze):

    A_in, B_in, C_in, D_in = [[] for n in range(4)]

    for j in range(n_samples):
        segnale = estrazione_segnale(input_file)
        A_in.append(segnale["A_in"])
        B_in.append(segnale["B_in"])
        C_in.append(segnale["C_in"])
        D_in.append(segnale["D_in"])
    
    # Valori medi delle misure ripetute
    A_in_media = np.average(A_in)
    A_int_std = np.std(A_in, ddof=1)

    B_in_media = np.average(B_in)
    B_int_std = np.std(B_in, ddof=1)

    C_in_media = np.average(C_in)
    C_int_std = np.std(C_in, ddof=1)

    D_in_media = np.average(D_in)
    D_int_std = np.std(D_in, ddof=1)

    # Calcolo ampiezza complessa del segnale in ingresso
    C_in = A_in_media - 1j*B_in_media
    dC_in = np.sqrt(A_in_std**2 + B_in_std**2)

    # Calcolo ampiezza complessa segnale in uscita
    C_in = A_in_media - 1j*B_in_media
    dC_in = np.sqrt(A_in_std**2 + B_in_std**2)

