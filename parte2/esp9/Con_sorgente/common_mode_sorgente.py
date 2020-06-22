from libphysics import *
from estrazione_segnale_funzione import *

freqs = [1e3, 2e3]
n_samples = 3 # cambiare in lab in base a quante misure ripetute si fanno!
output_file = "output_common_mode_sorgente.csv"
Rc = 10e3
Re = 120

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gcm = []
Gcm_fase = []
dGcm = []
dGcm_fase = []


for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    # preparo array per fit sinusoidale
    A_in = []
    B_in = []
    A_out = []
    B_out = []
    for j in range(n_samples):      # ciclo sulle misure ripetute
        input_file = "Serie_temporali/" + str(i+1) + "_" + str(j+1) + ".csv"
        # estraggo dal fit le ampiezze di un misura
        segnale = estrazione_segnale(input_file, freq)
        A_in.append(segnale["A_in"])
        B_in.append(segnale["B_in"])
        A_out.append(segnale["A_out"])
        B_out.append(segnale["B_out"])

    A_in_media = np.average(A_in)
    A_in_std = np.std(A_in)

    B_in_media = np.average(B_in)
    B_in_std = np.std(B_in)

    A_out_media = np.average(A_out)
    A_out_std = np.std(A_out)

    B_out_media = np.average(B_out)
    B_out_std = np.std(B_out)

    # calcolo ampiezze complesse per funzione di trasferimento e Gcm
    # chiamo C la ampiezza complessa del segnale (vedi slides)
    C_in = A_in_media - 1j*B_in_media
    dC_in = sqrt(A_in_std**2 + B_in_std**2)

    C_out = A_out_media - 1j*B_out_media
    dC_out = sqrt(A_out_std**2 + B_out_std**2)

    H = C_out/C_in
    dH = H * sqrt((dC_out/C_out)**2 + (dC_in/C_in)**2)
    Gcm.append(abs(H))
    dGcm.append(abs(dH))
    Gcm_fase.append(np.angle(H))
    dGcm_fase.append(abs(H))

# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([Gcm[i], Gcm_fase[i], dGcm[i], dGcm_fase[i]])
