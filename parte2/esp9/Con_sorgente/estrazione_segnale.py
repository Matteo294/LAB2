from libphysics import *

data_file = "data_serie_temporale.csv"
freq = 1e3
w = 2*pi*freq

[t, V] = readCSV(data_file, skiprows=10, untilrow=1900)
dV = numpify(1e-3, dim=len(V))

#funzioni per il fit
func1= numpify(np.ones(len(t)))     # costante
func2 = numpify(np.sin(w*t))        # seno

matrix = np.hstack([func1, func2])
fit = lsq_fit(V, matrix, dV)
coeffs = fit["fit_out"]

plt.plot(t, V)
plt.plot(t, coeffs[0]*func1 + coeffs[1]*func2)
plt.show()