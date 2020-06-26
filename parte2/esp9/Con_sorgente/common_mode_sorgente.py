from libphysics import *
from estrazione_segnale_funzione import *

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 63.3e3, 88.6e3, 141.3e3, 200e3, 159e3]
n_samples = 5 # cambiare in lab in base a quante misure ripetute si fanno!
output_file = "output_common_mode_sorgente.csv"

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gcm = []
Gcm_fase = []
dGcm = []
dGcm_fase = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    # preparo array per fit sinusoidale
    # A_in = []
    # B_in = []
    # A_out = []
    # B_out = []
    Gcm_temp = []
    Gcm_fase_temp = []
    for j in range(n_samples):      # ciclo sulle misure ripetute
        input_file = "Data_common_mode/" + str(i+1) + "/" + str(j+1) + ".csv"
        # estraggo dal fit le ampiezze di un misura
        plots=False
        if(j==0):
            plots = True

        [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
        [_, A_in, B_in] = segnale_in
        [_, A_out, B_out] = segnale_out
        C_in = A_in - 1j*B_in
        C_out = A_out - 1j*B_out
        H = C_out / C_in

        Gcm_temp.append(float(abs(H)))
        Gcm_fase_temp.append(float(np.angle(H)))
        
    Gcm.append(np.average(Gcm_temp))
    Gcm_fase.append(np.average(Gcm_fase_temp))
    dGcm.append(np.std(Gcm_temp, ddof=1))
    dGcm_fase.append(np.std(Gcm_fase_temp, ddof=1))
    
# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([Gcm[i], Gcm_fase[i], dGcm[i], dGcm_fase[i]])
