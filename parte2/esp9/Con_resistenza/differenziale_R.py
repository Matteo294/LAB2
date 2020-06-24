from libphysics import *
from estrazione_segnale_funzione import *


# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_R.csv"
[Gcm, Gcm_fase] = readCSV(file_gcm)

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3]
output_file = "output_common_mode_R.csv"
Rc = 10e3
Re = 120

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gdiff = []
Gdiff_fase = []


for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_diff/" + str(i+1) + ".csv"
    segnale = estrazione_segnale(input_file, freq)
    A_in = segnale["A_in"]
    B_in = segnale["B_in"]
    A_out = segnale["A_out"]
    B_out = segnale["B_out"]    

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H = C_out/C_in
    # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    Gdiff_complesso = H - Gcm[i]*np.exp(1j*Gcm_fase[i])
    Gdiff.append(float(abs(Gdiff_complesso)))
    Gdiff_fase.append(float(np.angle(Gdiff_complesso)))

Gdiff = numpify(Gdiff)
Gdiff_fase = numpify(Gdiff_fase)
Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
freqs = numpify(freqs)
    
print

#Cmrr = abs(Gdiff/Gcm)
#re_stima = Rc/(2*Gdiff) - Re

b1 = bodeplot(freqs, Amp=Gcm, Phase=Gcm_fase)
b2 = bodeplot(freqs, Amp=Gdiff, Phase=Gdiff_fase, color="blue", figure=b1)
[ax1, ax2] = b1.axes
ax1.legend(["Gcm", "Gdiff"])
ax2.legend(["Gcm", "Gdiff"])
plt.show()
