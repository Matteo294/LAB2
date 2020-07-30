''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np
import sys
from uncertainties import ufloat
from uncertainties import unumpy
from uncertainties.umath import *
from incertezza_H_funzione import *

plot_flag = 0
print_flag = 1
if len(sys.argv) > 1:
    plot_flag = int(sys.argv[1])

fileimmagine = '../esp9/Immagini/Faraday.png'

def Gdiff(f):
    w = 2*pi*f
    # carico = oscilloscopio
    Cosc = 110e-12
    Rosc = 1e6
    Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
    # parametri che servono
    Rc = 9.830e3
    Re = 119.25
    re = 45
    # Gdiff
    G0_diff = Rc / (2*(re + Re))
    return (G0_diff*Zosc)/(Rc + Zosc)

def dGdiff(f):
    return Gdiff(f) / 100

# Misure e costanti
R_lim = 10
dR_lim = 0.001
n_samples = 5 # numero misure ripetute

# Parametri bobine
mu0 = 4*pi*1e-7
sigma1 = np.pi * (17.5e-3/2)**2
sigma2 = sigma1
n_S = 30
n_R = 28
M_dipolo = lambda d: 2 * mu0/(4*np.pi) * n_S * n_R * sigma1 * sigma2 / d**3

distanze = numpify([1.35e-2, 2.3e-2, 4.6e-2, 10.5e-2, 4.4e-2, 1.8e-2])
d_distanze = [1 for _ in range(len(distanze))]
frequenze = [1e3, 50e3, 150e3]
omegas = [2*np.pi*f for f in frequenze]

t_schermo = numpify([2e-4, 3.7e-6, 1.5e-6])
''' Analisi accoppiamento '''
Z_ctrl = []
dZ_ctrl_fase = []
dZ_ctrl_abs = []
for i, f in enumerate(frequenze):

    Z = []
    dH_abs = []
    dH_fase = []
    for j in range(n_samples):

        filename = "Data/h/f" + str(i+1) + "/" + str(j+1) + ".csv"
        s_in, s_out, ds_in, ds_out = estrazione_segnale(filename, f, showplots=0)
        [_, A_in, B_in] = s_in 
        [_, A_out, B_out] = s_out 
        [_, dA_in, dB_int] = ds_in
        [_, dA_out, dB_out] = ds_out

        # Calcolo ampiezza complessa del segnale in ingresso
        C_in = A_in - 1j*B_in

        # Calcolo ampiezza complessa segnale in uscita
        C_out = A_out - 1j*B_out
        
        
        # Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)
        Z.append(C_out / C_in / Gdiff(f) * R_lim) 
        dH = incertezza_H(C_in, C_out, t_schermo[i], freqs=frequenze[i])
        dH_abs.append(dH["abs"])
        dH_fase.append(dH["arg"])

    dH_abs = np.mean(numpify(dH_abs))
    dH_fase = np.mean(numpify(dH_fase))
    Z = numpify(Z)
    Z_ctrl.append(np.mean(Z))
    dZ_ctrl_abs.append(np.std(abs(Z), ddof=1) + dH_abs)
    dZ_ctrl_fase.append(np.std(np.angle(Z), ddof=1) + dH_fase) 

Z_ctrl = numpify(Z_ctrl)
params = linreg(omegas, np.imag(Z_ctrl), dZ_ctrl_abs)
M_ctrl = params['m']
dM_ctrl = params['dm']
#print("Induttanza di controllo accoppiamento\n\tM_ctrl={}+-{}".format(M_ctrl, dM_ctrl))
'''-----------------------------------'''

# Array induttanza mutua a diverse distanze
Mrs = [] 
dMrs = []
Zeff_d = [] #array di array: per ogni distanza ho un array di Zeff in funzione della frequenza
dZeff_d_abs = []
dZeff_d_fase = []


for d in range(len(distanze)):

    base_input_file = "Data/d" + str(d+1)
    Z_eff = []
    dZ_eff = []
    Zeff_complessa = []
    dZeff_abs = []
    dZeff_fase = []

    for i, f in enumerate(frequenze):
        Z = []
        Z_complessa = []
        dH_abs = []
        dH_fase = []
        for j in range(n_samples):

            filename = base_input_file + "/f" + str(i+1) + "/" + str(j+1) + '.csv'
            plots=0
            if d==2:
                plots = 0
            s_in, s_out, ds_in, ds_out = estrazione_segnale(filename, f, showplots=plots)
            [_, A_in, B_in] = s_in 
            [_, A_out, B_out] = s_out 
            [_, dA_in, dB_int] = ds_in
            [_, dA_out, dB_out] = ds_out

            # Calcolo ampiezza complessa del segnale in ingresso
            C_in = A_in - 1j*B_in

            # Calcolo ampiezza complessa segnale in uscita
            C_out = A_out - 1j*B_out

            # Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)
            Z_complessa.append(C_out / C_in / Gdiff(f) * R_lim)
            Z.append(np.imag(C_out / C_in / Gdiff(f) * R_lim)) # (Z dovrebbe essere puramente immaginaria, c'è solo la mutua induzione)
            dH = incertezza_H(C_in, C_out, t_schermo[i], freqs=frequenze[i])
            dH_abs.append(dH["abs"])
            dH_fase.append(dH["arg"])
        # complesse per il bodeplot
        dH_abs = np.mean(numpify(dH_abs))
        dH_fase = np.mean(numpify(dH_fase))
        Zeff_complessa.append(np.mean(numpify(Z_complessa)))
        dZeff_abs.append(np.std(abs(numpify(Z_complessa))) + dH_abs)
        dZeff_fase.append(np.std(abs(numpify(Z_complessa))) + dH_fase)
        # parte immaginaria per il fit
        Z_eff.append(np.mean(Z))
        dZ_eff.append(np.std(Z, ddof=1) + dH_abs) # Non solo questa, anche 1% fondoscala

    
    Zeff_complessa = numpify(Zeff_complessa) - Z_ctrl
    dZeff_abs = numpify(dZeff_abs)
    dZeff_fase = numpify(dZeff_fase)
    Zeff_d.append(Zeff_complessa)
    dZeff_d_abs.append(dZeff_abs)
    dZeff_d_fase.append(dZeff_fase)

    Ze = linreg(omegas, Z_eff, dZ_eff)


    Mrs.append(Ze['m'])
    dMrs.append(np.sqrt(Ze['dm']**2 + dM_ctrl**2))
    

# modelli teorici Zeff = j*f/2pi*Mrs
f = np.logspace(3, 6)
Zeff_teo = [f*1j*2*pi *2 * mu0/(4*np.pi) * n_S * n_R * sigma1 * sigma2 / d**3 for d in distanze]

# bodeplot Zeff in funzione della frequenza
b00 = bodeplot(frequenze, H = Z_ctrl, err=1, Amperr=numpify(dZ_ctrl_abs), Phaseerr=numpify(dZ_ctrl_fase), logyscale=1, color="black")
b0 = bodeplot(frequenze, H = Zeff_d[0], err=1, Amperr=dZeff_d_abs[0], Phaseerr=dZeff_d_fase[0], logyscale=1, color="firebrick", figure=b00)
b1 = bodeplot(frequenze, H = Zeff_d[1], err=1, Amperr=dZeff_d_abs[1], Phaseerr=dZeff_d_fase[1], figure=b0, color="royalblue")
b2 = bodeplot(frequenze, H = Zeff_d[2], err=1, Amperr=dZeff_d_abs[2], Phaseerr=dZeff_d_fase[2], figure=b1, color = "darkorange")
b3 = bodeplot(frequenze, H = Zeff_d[3], err=1, Amperr=dZeff_d_abs[3], Phaseerr=dZeff_d_fase[3], figure=b2, color="mediumseagreen")
b4 = bodeplot(frequenze, H = Zeff_d[4], err=1, Amperr=dZeff_d_abs[4], Phaseerr=dZeff_d_fase[4],figure=b3, color="yellowgreen")
b5 = bodeplot(frequenze, H = Zeff_d[5], err=1, Amperr=dZeff_d_abs[5], Phaseerr=dZeff_d_fase[5],figure=b4, color="red")
# modello teorico
b6 = bodeplot(f, H=Zeff_teo[0], color = "firebrick", asline=1, figure=b5)
b7 = bodeplot(f, H=Zeff_teo[1], color = "royalblue", asline=1, figure=b6)
b8 =  bodeplot(f, H=Zeff_teo[2], color = "darkorange", asline=1, figure=b7)
b9 = bodeplot(f, H=Zeff_teo[3], color = "mediumseagreen", asline=1, figure=b8)
b10 = bodeplot(f, H=Zeff_teo[4], color = "yellowgreen", asline=1, figure=b9)
b11 = bodeplot(f, H=Zeff_teo[5], color = "red", asline=1, figure=b10)
ax = b5.axes[0]
handles,_ = ax.get_legend_handles_labels()
labels = [r"$d=${} mm".format(distanze[i]*1000) for i in range(len(distanze))]
labels.append("Controllo")
b3.legend(handles, labels=labels, loc='lower center', ncol=4)
#plt.show()


# Leggo i valori dei modelli teorici
d, val, approx = readCSV('Mrs/induzione.csv')


# Cambio udm ai dati
d = [dist/1e-3 for dist in d]
distanze = [dist/1e-3 for dist in distanze]
approx = [a/1e-6 for a in approx]
val = [v/1e-6 for v in val]
Mrs = [m/1e-6 for m in Mrs]
dMrs = [np.real(d)/1e-6 for d in dMrs]

# Plot dati, salvo nella cartella immagini
fig = plt.figure()
plt.semilogy(d, approx, label="Approssimazione dipolo", linewidth=1.8, color = 'red', ls='--') #color=[0, 1, 0])
plt.semilogy(d, val, label="Formula di Neumann", linewidth=1.8, color = 'red')#color='blue')
plt.errorbar(distanze, Mrs, yerr=dMrs, xerr=d_distanze, fmt='.', markersize=8, color='black', label="Dati sperimentali")#markerfacecolor='red'
plt.xlabel(r'Distanza  [mm]')
plt.ylabel(r'$M_{RS}$  $[\mu H]$')
plt.ylim((5e-4, 1e1))
plt.title("Mutua induzione tra bobine", fontsize=20)
#plt.legend()

fig.legend(loc='lower center', ncol=3)
# adjust layout
#fig.tight_layout()
#plt.savefig(fileimmagine, bbox_inches='tight', dpi=1000)
plt.show()
