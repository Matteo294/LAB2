from libphysics import *
from estrazione_segnale_funzione import *
from incertezza_H_funzione import *
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *

plt.rc('font', family='serif', size=18)


# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_R.csv"
[Gcm, Gcm_fase, dGcm, dGcm_fase] = readCSV(file_gcm)

freqs = [100, 250, 520, 980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3]
t_schermo=numpify([2.6e-3, 1e-3, 4.5e-4, 2.5e-4, 6.1e-5, 1.9e-5, 5.7e-5, 1e-6, 1e-6, 1e-6, 1e-6])
Rc = ufloat(9.830, 0.001)*1e3
Re = ufloat(119.25, 0.03)
R1 = 9.924e3

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gdiff = []
dGdiff = []
Gdiff_fase = []
dGdiff_fase = []
H = []
dH_amp =[]
dH_fase = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_diff/" + str(i+1) + ".csv"
    [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
    [_, A_in, B_in] = segnale_in
    [_, A_out, B_out] = segnale_out
    [_, dA_in, dB_in] = dsegnale_in
    [_, dA_out, dB_out] = dsegnale_out  

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H.append(np.complex(-C_out/C_in))
    dH = incertezza_H(C_in, C_out, t_schermo=t_schermo[i], freqs=freq)
    dH_amp.append(dH["abs"]) 
    dH_fase.append(dH["arg"])

    # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    Gdiff_complesso = -H[i] + 1/2*Gcm[i]*np.exp(1j*Gcm_fase[i]/2)
    Gdiff.append(float(abs(Gdiff_complesso)))
    Gdiff_fase.append(normalize_angle(np.angle(Gdiff_complesso)))
    dGdiff.append(dH_fase[i])


Gdiff = numpify(Gdiff)
Gdiff_fase = numpify(Gdiff_fase)
Gcm = numpify(Gcm)
Gcm_fase = numpify(Gcm_fase)
freqs = numpify(freqs)
H = numpify(H)
dH_amp=numpify(dH_amp)
dH_fase = numpify(dH_fase)



#stima Cmrr e re
Gdiff_u = unumpy.uarray(Gdiff, dGdiff)
Gcm_u = unumpy.uarray(Gcm, dGcm)
Cmrr = np.average(unumpy.nominal_values(Gdiff_u/Gcm_u))
dCmrr = np.std(unumpy.nominal_values(Gdiff_u/Gcm_u))
re = Rc/(2*Gdiff_u) - Re
re_stima = np.mean(unumpy.nominal_values(re)[:6])
dre_stima = np.std(unumpy.nominal_values(re)[:6])

print("re = {}+-{}".format(re_stima, dre_stima))
print("cmrr = {}+-{}".format(Cmrr, dCmrr))

# modelli teorici
Cosc = 110e-12
Rosc = 1e6
G_0_diff = Rc.n / (2*(re_stima + Re.n))
G_0_cm = -Rc.n/(2*R1)
f = np.logspace(2, 6)
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Gdiff_teo = (G_0_diff*Zosc)/(Rc.n + Zosc)
Gcm_teo = (G_0_cm*Zosc)/(Rc.n + Zosc)

# plots
b1 = bodeplot(freqs, Amp = Gcm, Phase = Gcm_fase, err=1, Amperr=dGcm, Phaseerr=dGcm_fase, deg=False, logyscale=True)
b2 = bodeplot(freqs, H=H, err=1, Amperr=dH_amp, Phaseerr=dH_fase, color="blue", figure=b1)
b3 = bodeplot(f, H=Gcm_teo, figure = b2, asline=True, linestyle='--', color="red")
b4 = bodeplot(f, H=Gdiff_teo, figure = b3, asline=True, linestyle = '--', color = "blue")
# set legend
ax = b3.axes[0]
handles,_ = ax.get_legend_handles_labels()
b3.legend(handles, labels=[r"Modello $G_{CM}$", r"Modello $G_{DIFF}$", r"$G_{CM}$", r"$G_{DIFF}$"], loc='lower center', ncol=2)
# adjust layout
b3.tight_layout()
b3.subplots_adjust(hspace=0.4, left=0.1)
plt.show()

print("Stima re: {} +- {}".format(re_stima, dre_stima))