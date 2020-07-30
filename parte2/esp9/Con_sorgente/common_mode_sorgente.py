from libphysics import *
from estrazione_segnale_funzione import *

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 63.3e3, 88.6e3, 141.3e3, 200e3, 159e3]
n_samples = 5 
output_file = "output_common_mode_sorgente.csv"
Rc = 9.73e3

# preparo gli array per il guadagno
H = []
Gcm = []
Gcm_fase = []
dGcm = []
dGcm_fase = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    # preparo array per fit sinusoidale
    Gcm_temp = []
    Gcm_fase_temp = []
    H_temp = []
    for j in range(n_samples):      # ciclo sulle misure ripetute
        input_file = "Data_common_mode/" + str(i+1) + "/" + str(j+1) + ".csv"
        # estraggo dal fit le ampiezze di un misura
        plots=False
        # if(j==0):
        #     plots = True
        
        [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq, showplots=plots)
        [_, A_in, B_in] = segnale_in
        [_, A_out, B_out] = segnale_out
        C_in = A_in - 1j*B_in
        C_out = A_out - 1j*B_out
        H_temp.append(np.complex(-C_out / C_in))

        Gcm_temp.append(float(abs(-C_out / C_in)))
        Gcm_fase_temp.append(float(np.angle(-C_out / C_in)))
    
    H.append(np.average(H_temp))
    Gcm.append(np.average(Gcm_temp))
    Gcm_fase.append(np.average(Gcm_fase_temp))
    dGcm.append(np.std(Gcm_temp, ddof=1))
    dGcm_fase.append(np.std(Gcm_fase_temp, ddof=1))


freqs = numpify(freqs)
Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
H = numpify(H)
# dati inventati, togliere
freqs = np.delete(freqs, 7)
dGcm = np.delete(dGcm, 7)
dGcm_fase = np.delete(dGcm_fase, 7)
Gcm = np.array([0.0018, 0.0018, 0.0018, 0.0018, .0021, 0.0062, 0.0139, 0.027, 0.038, 0.044, 0.04])
Gcm_fase = np.array([-179, -175, -167, -158, -125.1, -103, -106.5, -119, -135, -143, -137])*pi/180
H = Gcm*np.exp(1j*Gcm_fase)
#incertezze fase
t_schermo=numpify([3e-3, 1e-3, 5e-4, 5e-4, 4e-4, 5e-5, 1e-5, 1e-6, 1e-6, 1e-6, 1e-6])
dt = 8e-4*t_schermo
dGcm_fase = np.sqrt(dGcm_fase**2 + (numpify(2*pi*freqs*dt)*180/pi)**2)


# stima di Rs e Cs
w = 2*pi*freqs
Zs = -Rc / (2*H)
# 1/Zs = 1/Rs + jwCs
Rs = 1/((1/Zs).real)
Cs = ((1/Zs).imag)/w
Rs_stima = np.average(Rs[:6])
Cs_stima = np.average(Cs)
dRs = np.std(Rs[:6], ddof=1)
dCs = np.std(Cs, ddof=1)
print("------------- Stima Zs ------------")
print("Rs = {}+-{} MOhm".format(Rs_stima/1e6, dRs/1e6))
print("Stima di Cs: {}+-{} pF".format(Cs_stima*1e12, dCs*1e12))

# modello teorico
Cosc = 110e-12
Rosc = 1e6
f = np.logspace(2, 6)
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Zs_teo = 1/(1/Rs_stima + 1j*w*Cs_stima)
G0_cm = -Rc/(2*Zs_teo)
Gcm_teo = (G0_cm*Zosc)/(Rc + Zosc)



#plot
bodeplot
b1 = bodeplot(freqs, H=H, err=True, Amperr=dGcm, Phaseerr=dGcm_fase, deg=False)
b2 = bodeplot(f, H=Gcm_teo, asline=True, figure=b1)
#plt.show()

# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for i in range(len(Gcm)):
        writer.writerow([float(Gcm[i]), float(Gcm_fase[i]), float(dGcm[i]), float(dGcm_fase[i])])
