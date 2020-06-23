from libphysics import *
from estrazione_segnale_funzione import *

freqs = [1e3, 2e3]
n_samples = 3 # cambiare in lab in base a quante misure ripetute si fanno!
Rc = 10e3
Re = 120


# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_sorgente.csv"
[Gcm, Gcm_fase, dGcm, dGcm_fase] = readCSV(file_gcm)

# old [freq, Vin, Vout, fase] = readCSV(data_file, skiprows=1)
# preparo gli array per Gdiff
Gdiff = []
Gdiff_fase = []
dGdiff = []
dGdiff_fase = []


# estrazione ampiezze da serie temporali e ricostruzione Gdiff
for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    # preparo array per fit sinusoidale
    A_in = []
    B_in = []
    A_out = []
    B_out = []
    for j in range(n_samples):      # ciclo sulle misure ripetute
        input_file = "Data_diff/" + str(i+1) + "_" + str(j+1) + ".csv"
        # estraggo dal fit le ampiezze di un misura
        segnale = estrazione_segnale(input_file, freq)
        A_in.append(segnale["A_in"])
        B_in.append(segnale["B_in"])
        A_out.append(segnale["A_out"])
        B_out.append(segnale["B_out"])

    A_in_media = np.average(A_in)
    A_in_std = np.std(A_in, ddof=1)

    B_in_media = np.average(B_in)
    B_in_std = np.std(B_in, ddof=1)

    A_out_media = np.average(A_out)
    A_out_std = np.std(A_out, ddof=1)

    B_out_media = np.average(B_out)
    B_out_std = np.std(B_out, ddof=1)

    # calcolo ampiezze complesse per funzione di trasferimento e Gcm
    # chiamo C la ampiezza complessa del segnale (vedi slides)
    C_in = A_in_media - 1j*B_in_media
    dC_in = sqrt(A_in_std**2 + B_in_std**2)

    C_out = A_out_media - 1j*B_out_media
    dC_out = sqrt(A_out_std**2 + B_out_std**2)

    H = C_out/C_in
    # provvisorio dH = H * sqrt((dC_out/C_out)**2 + (dC_in/C_in)**2)

    # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    Gdiff_complesso = H - Gcm[i]*np.exp(1j*Gcm_fase[i])
    Gdiff.append(abs(Gdiff_complesso))
    Gdiff_fase.append(np.angle(Gdiff_complesso))
    # AGGIUNGERE INCERTEZZE
    

#Cmrr = abs(Gdiff/Gcm)
#re_stima = Rc/(2*Gdiff) - Re

b1 = bodeplot(freqs, Amp=Gcm, Phase=Gcm_fase)
b2 = bodeplot(freqs, Amp=Gdiff, Phase=Gdiff_fase, figure=b1, color="blue")
[ax1, ax2] = b1.axes
ax1.legend(["Gcm", "Gdiff"])
ax2.legend(["Gcm", "Gdiff"])
plt.show()

