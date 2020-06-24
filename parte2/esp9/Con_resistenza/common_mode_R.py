from libphysics import *
from estrazione_segnale_funzione import *


freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3]
output_file = "output_common_mode_R.csv"
Rc = 10e3
Re = 120

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gcm = []
Gcm_fase = []


for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_common_mode/" + str(i+1) + ".csv"
    segnale = estrazione_segnale(input_file, freq)
    A_in = segnale["A_in"]
    B_in = segnale["B_in"]
    A_out = segnale["A_out"]
    B_out = segnale["B_out"]    

    C_in = A_in - 1j*B_in

    C_out = A_out - 1j*B_out

    H = C_out/C_in

    Gcm.append(float(abs(H)))
    Gcm_fase.append(float(np.angle(H)))


Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
freqs = numpify(freqs)


# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([float(Gcm[i]), float(Gcm_fase[i])])

# bodeplot
# b1 = bodeplot(freqs, Amp = Gcm, Phase=Gcm_fase, deg=False)
# plt.show()

