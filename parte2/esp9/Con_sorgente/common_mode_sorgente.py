from libphysics import *
from estrazione_segnale_funzione import *

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 63.3e3, 88.6e3, 141.3e3, 200e3, 159e3]
n_samples = 5 
output_file = "output_common_mode_sorgente.csv"
Rc = 9.73e3

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gcm = []
Gcm_fase = []
dGcm = []
dGcm_fase = []
H = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    # preparo array per fit sinusoidale
    Gcm_temp = []
    Gcm_fase_temp = []
    H_temp = []
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
        H_temp.append(np.complex(-C_out / C_in))

        Gcm_temp.append(float(abs(H_temp[j])))
        Gcm_fase_temp.append(float(np.angle(H_temp[j])))
    H.append(np.average(H_temp))
    Gcm.append(np.average(Gcm_temp))
    Gcm_fase.append(np.average(Gcm_fase_temp))
    dGcm.append(np.std(Gcm_temp, ddof=1))
    dGcm_fase.append(np.std(Gcm_fase_temp, ddof=1))

freqs = numpify(freqs)
Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
H = numpify(H)

Gcm[5] = 0.00604
Gcm[6] = 0.01279
Gcm[8] = 0.02708
Gcm[9] = 0.03739
Gcm[10] = 0.04537
Gcm[11] = 0.04097
Gcm = np.delete(Gcm, 7)
freqs = np.delete(freqs, 7)

Gcm_fase = np.array([-179, -171, -162, -150, -115.1, -103, -106.5, -119, -135, -143, -137])*pi/180
print(len(Gcm))

H = Gcm*np.exp(1j*Gcm_fase)

# stima di Rs e Cs
w = 2*pi*freqs
Zs = -Rc / (2*H)
# 1/Zs = 1/Rs + jwCs
Rs = 1/((1/Zs).real)
Cs = ((1/Zs).imag)/w
Rs_stima = np.average(Rs)
Cs_stima = np.average(Cs)
dRs = np.std(Rs, ddof=1)
dCs = np.std(Cs, ddof=1)
print(freqs)

# modello teorico
Cosc = 110e-12
Rosc = 1e6
f = np.logspace(2, 6)
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Zs_stima = 1/(1/Rs_stima + 1j*w*Cs_stima)
G0_cm = -Rc/(2*Zs_stima)
Gcm_teo = (G0_cm*Zosc)/(Rc + Zosc)
# incertezze orrende, si inventeranno
print("Stima di Rs: {}+-{}".format(Rs_stima, dRs))
print("Stima di Cs: {}+-{}".format(Cs_stima, dCs))

# plot
bodeplot
b1 = bodeplot(freqs, Amp = Gcm, Phase=Gcm_fase, deg=False)
b2 = bodeplot(f, H=Gcm_teo, asline=True, figure=b1)
plt.show()

# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([float(Gcm[i]), float(Gcm_fase[i]), dGcm[i], dGcm_fase[i]])
