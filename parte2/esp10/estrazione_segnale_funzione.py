from libphysics import *

def estrazione_segnale(data_file, freq, showplots=False):
    """ 
    Prende la serie temporale di Vin e Vout, fa il fit a somma di seno e coseno, e ricava la funzione di trasferimento
    showplots e' False di default
    """
    w = 2*pi*freq

    # presa dati
    [t, Vin, Vout] = readCSV(data_file, skiprows=10, untilrow=1900)
    dVin = numpify(1e-3, dim=len(Vin))    # da cambiare in lab!!
    dVout = numpify(1e-3, dim=len(Vout))

    # funzioni per il fit
    func_const= numpify(np.ones(len(t)))   # costante
    func_sin = numpify(np.sin(w*t))        # seno
    func_cos = numpify(np.cos(w*t))        # coseno
    # le combino nella matrice da dare alla funzione
    matrix = np.hstack([func_const, func_cos, func_sin])

    # chiamo la funzione di fit
    fit_Vin = lsq_fit(Vin, matrix, dVin)        # Vin = C + A*sin(wt) + B*cos(wt)
    fit_Vout = lsq_fit(Vin, matrix, dVout)      # Vout = C + A*sin(wt) + B*cos(wt)
    [C_in, A_in, B_in] = fit_Vin["fit_out"]
    [C_out, A_out, B_out] = fit_Vout["fit_out"]

    # plots
    
    if(showplots):
        plt.plot(t, Vin, label="Vin")
        plt.plot(t, C_in + A_in*func_cos + B_in*func_sin, label="Vin fit")
    
        plt.plot(t, Vout, label="Vout")
        plt.plot(t, C_out + A_out*func_cos + B_out*func_sin, label="Vout fit")

        plt.legend()
        plt.show()

    return {"A_in":A_in, "B_in":B_in, "C_in":C_in, "A_out":A_out, "B_out":B_out, "C_out":C_out}
    


    return H
