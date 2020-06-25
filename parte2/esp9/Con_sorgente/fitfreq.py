from libphysics import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)

''' 
    Fit polinomio di grado 3 G(f) = A + Bf + Cf^2 + Df^3 per interpolare i punti del grafico
'''
input_file = "output_differential_mode_sorgente.csv"

freqs = [] 
H = []
with open(input_file) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        freqs.append(float(row[0]))
        H.append(complex(row[1]))

def fitfunc(f, a, b, c, d):
    return a + b*f + c*f**2 + d*f**3

params_real = curve_fit(fitfunc, freqs, np.real(H))
params_imag = curve_fit(fitfunc, freqs, np.imag(H))
print("Parametri Re(H):", params_real[0])
print("Parametri Im(H):", params_imag[0])
Ar, Br, Cr, Dr = params_real[0]
Ai, Bi, Ci, Di = params_imag[0]

esponenti = np.linspace(2, 6, 1000)
basi = np.full(1000, 10)
f = np.power(basi, esponenti)

index = np.where(f > 210e3)[0]
imax = min(index)

plt.plot(freqs, np.real(H), '.')
plt.plot(f[0:imax], fitfunc(f[0:imax], Ar, Br, Cr, Dr))
plt.xlabel('f [Hz]')
plt.ylabel(r'$\textit{Re}(H)$')
plt.title('Interpolazione Re(H)')
plt.show()

plt.plot(freqs, np.imag(H), '.')
plt.plot(f[0:imax], fitfunc(f[0:imax], Ai, Bi, Ci, Di))
plt.xlabel('f [Hz]')
plt.ylabel(r'$\textit{Im}(H)$')
plt.title('Interpolazione Im(H)')
plt.show()