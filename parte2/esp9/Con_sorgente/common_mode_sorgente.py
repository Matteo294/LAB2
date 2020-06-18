from libphysics import *
from estrazione_segnale_funzione import *

# cambiare nomi files in lab se diversi!
input_files = ["Serie_temporali/data_serie_temporale_"+ str(i+1)+".csv" for i in range(3)]
output_file = "output_common_mode_sorgente.csv"
Rc = 10e3
Re = 120

Gcm = []
Gcm_fase = [] # sarà in rad, non deg!!
freqs = [1e3, 2e3, 3e3]

for file, freq in zip(input_files, freqs):
    H = estrazione_segnale(file, freq)
    Gcm.append(float(abs(H)))
    Gcm_fase.append(float(np.angle(H)))

# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([Gcm[i], Gcm_fase[i]])
