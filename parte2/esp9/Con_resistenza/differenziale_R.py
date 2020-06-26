from libphysics import *
from estrazione_segnale_funzione import *
from incertezza_H_funzione import *
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *


# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_R.csv"
[Gcm, Gcm_fase, dGcm, dGcm_fase] = readCSV(file_gcm)

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3]
output_file = "output_common_mode_R.csv"
Rc = ufloat(9.830, 0.001)*1e3
Re = ufloat(119.25, 0.03)

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gdiff = []
dGdiff = []
Gdiff_fase = []
dGdiff_fase = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_diff/" + str(i+1) + ".csv"
    [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
    [_, A_in, B_in] = segnale_in
    [_, A_out, B_out] = segnale_out
    [_, dA_in, dB_in] = dsegnale_in
    [_, dA_out, dB_out] = dsegnale_out  

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H = C_out/C_in
    dH = incertezza_H(A_in, B_in, A_out, B_out, dA_in, dB_in, dA_out, dB_out)
    [dH_amp, dH_fase] = [dH["abs"], dH["arg"]]

    # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    Gdiff_complesso = H - Gcm[i]*np.exp(1j*Gcm_fase[i])
    Gdiff.append(float(abs(Gdiff_complesso)))
    Gdiff_fase.append(float(np.angle(Gdiff_complesso)))
    dGdiff.append(sqrt(dH_amp**2 + dGcm[i]**2))


Gdiff = numpify(Gdiff, column = False)
Gdiff_fase = numpify(Gdiff_fase, column = False)
Gcm = numpify(Gcm, column = False)
Gcm_fase = numpify(Gcm_fase, column = False)
freqs = numpify(freqs, column = False)

#stima Cmrr e re
Gdiff_u = unumpy.uarray(Gdiff, dGdiff)
Gcm_u = unumpy.uarray(Gcm, dGcm)
Cmrr = abs(Gdiff_u/Gcm_u)
re = Rc/(2*Gdiff_u) - Re
re_stima = np.mean(unumpy.nominal_values(re)[:6])
dre_stima = np.std(unumpy.nominal_values(re)[:6])

b1 = bodeplot(freqs, Amp=Gcm, Phase=Gcm_fase)
b2 = bodeplot(freqs, Amp=Gdiff, Phase=Gdiff_fase, color="blue", figure=b1)
[ax1, ax2] = b1.axes
ax1.legend(["Gcm", "Gdiff"])
ax2.legend(["Gcm", "Gdiff"])
plt.show()

print("Stima re: {} +- {}".format(re_stima, dre_stima))