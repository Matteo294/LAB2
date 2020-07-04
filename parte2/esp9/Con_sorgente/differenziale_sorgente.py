from libphysics import *
from estrazione_segnale_funzione import *
from incertezza_H_funzione import *
import uncertainties
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *

plt.rc('font', family='serif', size=18)

freqs = [980, 3.6e3, 11.4e3, 38.1e3, 88.6e3, 141.3e3, 200e3, 159e3, 100, 250, 520]
Rc = 9.830e3
Re = 119.25

# prendo Gcm dall'output dell'altro script
file_gcm = "output_common_mode_sorgente.csv"
file_gdelta = "output_differential_mode_sorgente.csv"
[Gcm, Gcm_fase, dGcm, dGcm_fase] = readCSV(file_gcm)

offset = 3
Gcm = numpify(Gcm)
Gcm = np.concatenate((Gcm[offset:], Gcm[:offset]))
freqs = numpify(freqs)
Gcm_fase = numpify(Gcm_fase)
Gcm_fase = np.concatenate((Gcm_fase[offset:], Gcm_fase[:offset]))

# preparo gli array per il guadagno (complesso, contiene dentro fase e ampiezza non separate)
Gdiff = []
dGdiff = []
Gdiff_fase = []
dGdiff_fase = []
Gdiff_complesso = []
trans = []

for i, freq in enumerate(freqs):    # ciclo sulle frequenze
    input_file = "Data_diff/" + str(i+1) + ".csv"
    [segnale_in, segnale_out, dsegnale_in, dsegnale_out] = estrazione_segnale(input_file, freq)
    [_, A_in, B_in] = segnale_in
    [_, A_out, B_out] = segnale_out
    [_, dA_in, dB_in] = dsegnale_in
    [_, dA_out, dB_out] = dsegnale_out    

    C_in = A_in - 1j*B_in
    C_out = A_out - 1j*B_out
    H = C_out/C_in
    dH = incertezza_H(A_in, B_in, A_out, B_out, dA_in, dB_in, dA_out, dB_out)
    [dH_amp, dH_fase] = [dH["abs"], dH["arg"]]

    # Gdiff = Vout/Vin - Gcm/2, ma sono numeri complessi
    trans.append(np.complex(H))
    Gdiff_complesso.append(H - Gcm[i]/2*np.exp(1j*Gcm_fase[i]/2))
    Gdiff.append(float(abs(Gdiff_complesso[i])))
    Gdiff_fase.append(float(np.angle(Gdiff_complesso[i])))
    dGdiff.append(sqrt(dH_amp**2 + dGcm[i]**2))

Gdiff = numpify(Gdiff, column = False)
Gdiff_complesso = numpify(Gdiff_complesso, column=False)
Gdiff_fase = numpify(Gdiff_fase, column = False)
Gcm = numpify(Gcm, column = False)
Gcm_fase = numpify(Gcm_fase, column = False)
freqs = numpify(freqs, column = False)
trans = numpify(trans)

#stima Cmrr e re
Gdiff_u = unumpy.uarray(Gdiff, dGdiff)
Gcm_u = unumpy.uarray(Gcm, dGcm)
Cmrr = abs(Gdiff_u/Gcm_u)
re = Rc/(2*Gdiff_u) - Re
re_stima = np.mean(unumpy.nominal_values(re)[:6])
dre_stima = np.std(unumpy.nominal_values(re)[:6])
print(re_stima)

# modello teorico
Cosc = 110e-12
Rosc = 1e6
G0_diff = Rc / (2*(45 + Re))
f = np.logspace(2, 6)
w = 2*pi*f
Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
Gdiff_teo = (G0_diff*Zosc)/(Rc + Zosc)
Zs_stima = 1/(1/5e6 + 1j*w*11e-12)
G0_cm = -Rc/(2*Zs_stima)
Gcm_teo = (G0_cm*Zosc)/(Rc + Zosc)


# plots
b1 = bodeplot(freqs, Amp = Gcm, Phase = Gcm_fase, deg=False, logyscale=True)
b2 = bodeplot(freqs, H=-trans, color="blue", figure=b1)
b3 = bodeplot(f, Amp=abs(Gcm_teo), Phase=np.angle(Gcm_teo), deg=False, figure = b2, asline=True, linestyle='--', color="red")
b4 = bodeplot(f, H=Gdiff_teo, figure = b3, asline=True, linestyle = '--', color = "blue")
# set legend
ax = b3.axes[0]
handles,_ = ax.get_legend_handles_labels()
b3.legend(handles, labels=[r"$G_{CM}$", r"$G_{DIFF}$", r"Modello $G_{CM}$", r"Modello $G_{DIFF}$"], loc='lower center', ncol=2)
# adjust layout
b3.tight_layout()
b3.subplots_adjust(hspace=0.4, left=0.1)
plt.show()

# output
outfile = open(file_gdelta, 'w+')
for f, G in zip(freqs, Gdiff_complesso):
    outfile.write(str(f) + "," + str(G[0]) + '\n')
