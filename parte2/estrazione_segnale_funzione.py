from libphysics import *

def estrazione_segnale(data_file, freq, showplots=False):
    """ 
    Prende la serie temporale di Vin e Vout, fa il fit a somma di seno e coseno, e ricava la funzione di trasferimento
    showplots e' False di default
    """
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
    fit_Vin = lsq_fit(Vin, matrix, dVin)        # Vin = C + B*sin(wt) + A*cos(wt)
    fit_Vout = lsq_fit(Vout, matrix, dVout)      # Vout = C + B*sin(wt) + A*cos(wt)
    [C_in, A_in, B_in] = fit_Vin["fit_out"]
    [dC_in, dA_in, dB_in] = fit_Vin["dfit_out"]
    [C_out, A_out, B_out] = fit_Vout["fit_out"]
    [dC_out, dA_out, dB_out] = fit_Vout["dfit_out"]

    # plots
    
    if(showplots):
        plt.plot(t, Vin, label="Vin")
        plt.plot(t, C_in + A_in*func_cos + B_in*func_sin, label="Vin fit")
    
        plt.plot(t, Vout, label="Vout")
        plt.plot(t, C_out + A_out*func_cos + B_out*func_sin, label="Vout fit")

        plt.legend()
        plt.show()

    return [fit_Vin["fit_out"], fit_Vout["fit_out"], fit_Vin["dfit_out"], fit_Vout["dfit_out"]]
    

