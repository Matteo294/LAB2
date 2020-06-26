from libphysics import *
from estrazione_segnale_funzione import *
from incertezza_H_funzione import *


freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3]
output_file = "output_common_mode_R.csv"
Rc = 10e3
Re = 120

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gcm = []
dGcm = []
Gcm_fase = []
dGcm_fase = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_common_mode/" + str(i+1) + ".csv"
    [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
    [_, A_in, B_in] = segnale_in
    [_, A_out, B_out] = segnale_out
    [_, dA_in, dB_in] = dsegnale_in
    [_, dA_out, dB_out] = dsegnale_out

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H = C_out/C_in
    dH = incertezza_H(A_in, B_in, A_out, B_out, dA_in, dB_in, dA_out, dB_out)

    Gcm.append(float(abs(H)))
    dGcm.append(dH["abs"])
    Gcm_fase.append(float(np.angle(H)))
    dGcm_fase.append(dH["arg"])


Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
dGcm = numpify(dGcm)
dGcm_fase = numpify(dGcm_fase)
freqs = numpify(freqs)


# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([float(Gcm[i]), float(Gcm_fase[i]), float(dGcm[i]), float(dGcm_fase[i])])

# bodeplot
# b1 = bodeplot(freqs, Amp = Gcm, Phase=Gcm_fase, deg=False, err=True, Amperr = dGcm, Phaseerr = dGcm_fase)
# plt.show()

