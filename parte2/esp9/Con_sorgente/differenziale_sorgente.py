from libphysics import *
from estrazione_segnale_funzione import *
from incertezza_H_funzione import *
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *

plt.rc('font', family='serif', size=18)
freqs = numpify([980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3, 159e3, 100, 250, 520])
Rc = 9.830e3
Re = 119.25

###########################
# PRESA DATI COMMON MODE E RIORDINAMENTO
###########################
# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_sorgente.csv"
file_gdelta = "output_differential_mode_sorgente.csv"
[Gcm, Gcm_fase, dGcm, dGcm_fase] = readCSV(file_gcm)
Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
dGcm = numpify(dGcm)
dGcm_fase = numpify(dGcm_fase)
H_cm = Gcm*np.exp(1j*Gcm_fase)

################################
# ANALISI DATI DIFFERENZIALE
################################

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gdiff = []
dGdiff = []
Gdiff_fase = []
dGdiff_fase = []
Gdiff_complesso = []
H = []
dH_amp = []
dH_fase = []

t_schermo=numpify([5e-4, 5e-4, 5e-5, 5e-5, 1e-6, 1e-6, 1e-6, 1e-6, 3e-3, 3e-3, 5e-4])

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_diff/" + str(i+1) + ".csv"
    [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
    [_, A_in, B_in] = segnale_in
    [_, A_out, B_out] = segnale_out
    [_, dA_in, dB_in] = dsegnale_in
    [_, dA_out, dB_out] = dsegnale_out    

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H.append(C_out/C_in)
    dH = incertezza_H(C_in, C_out, dVin_schermo=30e-3, dVout_schermo=400e-3, t_schermo=t_schermo[i], freqs=freqs[i])
    dH_amp.append(dH["abs"])
    dH_fase.append(dH["arg"])

    # # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    # Gdiff_complesso.append(H[i] - Gcm[i]/2*np.exp(1j*Gcm_fase[i]/2))
    # Gdiff.append(float(abs(Gdiff_complesso[i])))
    # Gdiff_fase.append(float(np.angle(Gdiff_complesso[i])))
    # dGdiff.append(sqrt(dH_amp[i]**2 + dGcm[i]**2))

H = numpify(H)
dH_amp = numpify(dH_amp)
dH_fase = numpify(dH_fase)
# !! i dati differenziale erano presi in ordine sbagliato, riordino
offset=8
freqs = np.concatenate((freqs[offset:], freqs[:offset]))
freqs = numpify(freqs)
H = np.concatenate((H[offset:], H[:offset]))
H = numpify(H)

Gdiff_complesso = H - Gcm/2*np.exp(1j*Gcm_fase/2)
Gdiff = np.abs(Gdiff_complesso)
Gdiff_fase = np.angle(Gdiff_complesso)
dGdiff = np.sqrt(dH_amp**2 + dGcm**2)


###########################
# STIME PARAMETRI E MODELLI TEORICI
###########################

#stima Cmrr e re
Gdiff_u = unumpy.uarray(Gdiff, dGdiff)
Gcm_u = unumpy.uarray(Gcm, dGcm)
Cmrr = np.average(unumpy.nominal_values(Gdiff_u[:4]/Gcm_u[:4]))
dCmrr = np.std(unumpy.nominal_values(Gdiff_u/Gcm_u))
re = Rc/(2*Gdiff_u) - Re
re_stima = np.mean(unumpy.nominal_values(re)[:4])
dre_stima = np.std(unumpy.nominal_values(re)[:4])
print("----------- Stime parametri------------")
print("re = {}+-{} Ohm".format(re_stima, dre_stima))
print("cmrr = {}+-{}".format(Cmrr, dCmrr))

# stima di Rs e Cs
w = 2*pi*freqs
Zs = -Rc / (2*H_cm)
# 1/Zs = 1/Rs + jwCs
Rs = 1/((1/Zs).real)
Cs = ((1/Zs).imag)/w
Rs_stima = np.average(Rs[:6])
Cs_stima = np.average(Cs[6:])
dRs = np.std(Rs[:6], ddof=1)
dCs = np.std(Cs[6:], ddof=1)
print("Rs = {}+-{} MOhm".format(Rs_stima/1e6, dRs/1e6))
print("Stima di Cs: {}+-{} pF\n\n".format(Cs_stima*1e12, dCs*1e12))

# modello teorico differenziale
Cosc = 110e-12
Rosc = 1e6
f = np.logspace(2, 6)
G0_diff = Rc / (2*(re_stima + Re))
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Gdiff_teo = (G0_diff*Zosc)/(Rc + Zosc)

# modello teorico common mode
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Zs_teo = 1/(1/Rs_stima + 1j*w*Cs_stima)
G0_cm = -Rc/(2*Zs_teo)
Gcm_teo = (G0_cm*Zosc)/(Rc + Zosc)


# plots
b1 = bodeplot(freqs, Amp = Gcm, Phase = Gcm_fase, err=1, Amperr=dGcm, Phaseerr=dGcm_fase, deg=False, logyscale=True)
b2 = bodeplot(freqs, H=-Gdiff_complesso, err = 1, Amperr = dH_amp, Phaseerr=dH_fase, color="blue", figure=b1)
b3 = bodeplot(f, Amp=abs(Gcm_teo), Phase=np.angle(Gcm_teo), deg=False, figure = b2, asline=True, linestyle='--', color="red")
b4 = bodeplot(f, H=Gdiff_teo, figure = b3, asline=True, linestyle = '--', color = "blue")
# set legend
ax = b3.axes[0]
handles,_ = ax.get_legend_handles_labels()
b3.legend(handles, labels=[r"Modello $G_{CM}$", r"Modello $G_{DIFF}$", r"$G_{CM}$", r"$G_{DIFF}$"], loc='lower center', ncol=2)
# adjust layout
b3.tight_layout()
b3.subplots_adjust(hspace=0.4, left=0.1)
plt.show()

# output
outfile = open(file_gdelta, 'w+')
for f, G in zip(freqs, Gdiff_complesso):
    outfile.write(str(f) + "," + str(G[0]) + '\n')
