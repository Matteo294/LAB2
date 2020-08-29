''' Esperienza di Faraday '''

from libphysics import *
from estrazione_segnale_funzione import *
import numpy as np
import sys
from uncertainties import ufloat
from uncertainties import unumpy, wrap
from uncertainties.umath import *
from incertezza_H_funzione import *

wangle = wrap(np.angle)
wreal = wrap(np.real)
wimag = wrap(np.imag)

def mod(Re, Im):
    return np.sqrt(Re**2 + Im**2)
wabs = wrap(mod)

def ang(Re, Im):
    return np.angle(Re+1j*Im)
wfase = wrap(ang)

plot_flag = 0
print_flag = 1

if len(sys.argv) > 1:
    plot_flag = int(sys.argv[1])

fileimmagine = '../esp9/Immagini/Faraday.png'

def Gdiff(f, flag='val'):
    w = 2*pi*f
    # carico = oscilloscopio
    Cosc = 110e-12
    Rosc = 1e6
    Zosc = Rosc/(1 + 1j*w*Rosc*Cosc)
    # parametri che servono
    Rc = ufloat(9.830e3, 0.002e3)
    Re = ufloat(119.25, 0.002e3)
    re = ufloat(45, 1)
    # Gdiff
    G0_diff = Rc / (2*(re + Re))
    G_abs = (G0_diff*np.absolute(Zosc))/ufloat(np.absolute(Rc.n + Zosc), Rc.s)
    G_fase = ufloat( np.angle((G0_diff.n*Zosc)/(Rc.n + Zosc)), 0)
    if flag == 'df':
        return G_abs.s
    else:
        return G_abs.n*np.exp(1j*G_fase.n) 

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

distanze = numpify([1.35e-2, 2.3e-2, 4.6e-2, 10.5e-2, 4.4e-2, 1.8e-2]) #m
d_distanze = [0.7 for _ in range(len(distanze))] #mm incertezza aggiunta dalla difficoltà di allineamento
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
        Z_complessa_abs = []
        Z_complessa_fase = []
        dH_abs = []
        dH_fase = []
        dGdiff_abs = []
        for j in range(n_samples):

            filename = base_input_file + "/f" + str(i+1) + "/" + str(j+1) + '.csv'
            plots=0
            if d==2:
                plots = 0
            s_in, s_out, ds_in, ds_out = estrazione_segnale(filename, f, showplots=plots)
            [_, A_in, B_in] = s_in 
            [_, A_out, B_out] = s_out 
            [_, dA_in, dB_in] = ds_in
            [_, dA_out, dB_out] = ds_out

            _A_in = ufloat(A_in, dA_in)
            _A_out = ufloat(A_out, dA_out)
            _B_in = ufloat(B_in, dB_in)
            _B_out = ufloat(B_out, dB_out)

            # Calcolo ampiezza complessa del segnale in ingresso
            C_in = A_in - 1j*B_in

            # Calcolo ampiezza complessa segnale in uscita
            C_out = A_out - 1j*B_out

            C_in_abs = ufloat(np.absolute(C_in), 1/np.absolute(C_in) * np.sqrt( (A_in*dA_in)**2 + (B_in*dB_in)**2 ) )
            C_out_abs = ufloat(np.absolute(C_out), 1/np.absolute(C_out) * np.sqrt( (A_out*dA_out)**2 + (B_out*dB_out)**2 ) )
            C_in_fase = wfase(_A_in, -1*_B_in)
            C_out_fase =  wfase(_A_out, -1*_B_out)

            # H e incertezze
            H_abs = C_out_abs/C_in_abs
            H_fase = C_out_fase - C_in_fase
            # G e incertezze
            dGdiff_abs = Gdiff(f, 'df')
            dGdiff_fase = 0
            Gdiff_abs = ufloat(np.abs(Gdiff(f)), dGdiff_abs)
            Gdiff_fase = ufloat(np.angle(Gdiff(f)), dGdiff_fase)

            ''' Impedenza efficace (vediamo lo spazio tra le due bobine come un induttore di induttanza Mrs)'''
            # Z e incertezze
            Z_complessa_abs.append(H_abs / Gdiff_abs * R_lim)
            Z_complessa_fase.append(H_fase - Gdiff_fase)

        # Preparo pesi per medie pesate
        w_abs = [1/v.s**2 for v in Z_complessa_abs]
        w_fase = [1/v.s**2 for v in Z_complessa_fase]
        # Medie pesate
        Z_eff_abs = np.average([v.n for v in Z_complessa_abs], weights=w_abs)
        d_abs1 = np.sqrt(1/np.sum(w_abs))
        Z_eff_fase = np.average([v.n for v in Z_complessa_fase], weights=w_fase)
        d_fase1 = np.sqrt(1/np.sum(w_fase))
        # Deviazioni standard
        d_abs2 = np.std([v.n for v in Z_complessa_abs], ddof=1)
        d_fase2 = np.std([v.n for v in Z_complessa_fase], ddof=1)
        # Combino gli errori
        d_abs = np.sqrt(d_abs1**2 + d_abs2**2)
        d_fase = np.sqrt(d_fase1**2 + d_fase2**2)
        # Carico negli array per il bodeplot
        Zeff_complessa.append(Z_eff_abs*np.exp(1j*Z_eff_fase))
        dZeff_abs.append(d_abs)
        dZeff_fase.append(d_fase)
        # parte immaginaria per il fit
        Z_eff.append(wimag(Z_eff_abs*np.exp(1j*Z_eff_fase)))
        dZ_eff.append(5*np.sqrt( (np.sin(Z_eff_fase)*d_abs)**2 + (Z_eff_abs*np.cos(Z_eff_fase)*d_fase)**2 ))
    
    Zeff_complessa = numpify(Zeff_complessa) - Z_ctrl
    dZeff_abs = numpify(dZeff_abs)
    dZeff_fase = numpify(dZeff_fase)
    Zeff_d.append(Zeff_complessa)
    dZeff_d_abs.append(dZeff_abs)
    dZeff_d_fase.append(dZeff_fase)

    print('distanza', distanze[d], 'frequenza', f,'fasi', np.angle(Zeff_complessa)*180/np.pi)

    Ze = linreg(omegas, Z_eff, dZ_eff)

    Mrs.append(Ze['m'])
    dMrs.append(np.sqrt(Ze['dm']**2 + dM_ctrl**2))   

# modelli teorici Zeff = j*f/2pi*Mrs
f = np.logspace(3, 6)
Zeff_teo = [f*1j*2*pi *2 * mu0/(4*np.pi) * n_S * n_R * sigma1 * sigma2 / d**3 for d in distanze]

# Bodeplot Zeff in funzione della frequenza
b00 = bodeplot(frequenze, H = Z_ctrl, err=1, Amperr=numpify(dZ_ctrl_abs), Phaseerr=numpify(dZ_ctrl_fase), logyscale=1, color="black")
b0 = bodeplot(frequenze, H = Zeff_d[0], err=1, Amperr=dZeff_d_abs[0], Phaseerr=dZeff_d_fase[0], logyscale=1, color="firebrick", figure=b00)
b1 = bodeplot(frequenze, H = Zeff_d[1], err=1, Amperr=dZeff_d_abs[1], Phaseerr=dZeff_d_fase[1], figure=b0, color="royalblue")
b2 = bodeplot(frequenze, H = Zeff_d[2], err=1, Amperr=dZeff_d_abs[2], Phaseerr=dZeff_d_fase[2], figure=b1, color = "darkorange")
b3 = bodeplot(frequenze, H = Zeff_d[3], err=1, Amperr=dZeff_d_abs[3], Phaseerr=dZeff_d_fase[3], figure=b2, color="mediumseagreen")
b4 = bodeplot(frequenze, H = Zeff_d[4], err=1, Amperr=dZeff_d_abs[4], Phaseerr=dZeff_d_fase[4],figure=b3, color="yellowgreen")
b5 = bodeplot(frequenze, H = Zeff_d[5], err=1, Amperr=dZeff_d_abs[5], Phaseerr=dZeff_d_fase[5],figure=b4, color="red")
# Modello teorico
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
plt.show()


# Leggo i valori dei modelli teorici
d, val, approx = readCSV('Mrs/induzione.csv')


# Cambio udm ai dati
d = [dist/1e-3 for dist in d]
distanze = [dist/1e-3 for dist in distanze]
approx = [a/1e-6 for a in approx]
val = [v/1e-6 for v in val]
Mrs = [m/1e-6 for m in Mrs]
dMrs = [np.real(dd)/1e-6 for dd in dMrs]

# Plot dati, salvo nella cartella immagini
fig = plt.figure()
plt.semilogy(d, approx, label="Approssimazione dipolo", linewidth=1.8, color = 'red', ls='--')
plt.semilogy(d, val, label="Formula di Neumann", linewidth=1.8, color = 'red')
plt.errorbar(distanze, Mrs, yerr=dMrs, xerr=d_distanze, fmt='.', markersize=8, color='black', label="Dati sperimentali")
plt.xlabel(r'Distanza  [mm]')
plt.ylabel(r'$M_{RS}$  $[\mu H]$')
plt.ylim((5e-4, 1e1))
plt.title("Mutua induzione tra bobine", fontsize=20)
plt.grid()

fig.legend(loc='lower center', ncol=3)
# adjust layout
fig.tight_layout()
plt.savefig(fileimmagine, bbox_inches='tight', dpi=1000)
plt.show()
