from libphysics import *
from scipy.optimize import curve_fit as fit

def fitfunc(t, C, A, B):
    return C + A*np.cos(w*t) + B*np.sin(w*t)

def estrazione_segnale(data_file, freq, showplots=False):
    """ 
    Prende la serie temporale di Vin e Vout, fa il fit a somma di seno e coseno, e ricava la funzione di trasferimento
    showplots e' False di default
    """
    global w 
    w = 2*pi*freq
	
    # presa dati
    # offset e dt dalla prima riga
    filetoread = os.path.join(data_file)
    with open(filetoread, 'r') as f:
        lines = f.readlines()
        line = lines[1]
        row = line.split(",")
        t_offset = float(row[3])
        dt = float(row[4])
    # dati veri e propri
    [t, Vin, Vout] = readCSV(data_file, skiprows=2, cols=[0,1,2])
    Vin = numpify(Vin)
    Vout = numpify(Vout)
    t = numpify(t)
    t = t_offset + t*dt
    dVin = float(max(Vin))/100
    dVout = float(max(Vout))/100
    dVin = numpify(dVin, dim=len(Vin))    
    dVout = numpify(dVout, dim=len(Vout))

    # funzioni per il fit
    func_const= numpify(np.ones(len(t)), column=True)   # costante
    func_sin = numpify(np.sin(w*t), column=True)        # seno
    func_cos = numpify(np.cos(w*t), column=True)        # coseno
    # le combino nella matrice da dare alla funzione
    matrix = np.hstack([func_const, func_cos, func_sin])

    # chiamo la funzione di fit
    fit_in = fit(fitfunc, t, Vin, sigma=dVin, absolute_sigma=True)
    fit_out = fit(fitfunc, t, Vout, sigma=dVout, absolute_sigma=True)
    [C_in, A_in, B_in] = fit_in[0]
    [dC_in, dA_in, dB_in] = np.sqrt([fit_in[1][0][0], fit_in[1][1][1], fit_in[1][2][2]])
    [C_out, A_out, B_out] = fit_out[0]
    [dC_out, dA_out, dB_out] = np.sqrt([fit_out[1][0][0], fit_out[1][1][1], fit_out[1][2][2]])

    # plots
    
    if(showplots):
        plt.plot(t, Vin, label="Vin")
        plt.plot(t, C_in + A_in*func_cos + B_in*func_sin, label="Vin fit")
    
        plt.plot(t, Vout, label="Vout")
        plt.plot(t, C_out + A_out*func_cos + B_out*func_sin, label="Vout fit")

        plt.legend()
        plt.show()

    return [[C_in, A_in, B_in], [C_out, A_out, B_out], [dC_in, dA_in, dB_in], [dC_out, dA_out, dB_out]]