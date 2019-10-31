from labbclass import LinearFit
import os
import math

# !!!!! Non ho ancora traasferito le incertezze (non ho copiato la parte di funzione) !!!!!!

file1 = 'misure/monte_10V_5mA.csv'
file2 = 'misure/monte_50V_5mA.csv'
file3 = 'misure/valle_10V_500uA.csv'
file4 = 'misure/valle_50V_5mA.csv'

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Valore di fondo scala / n tacche / sqrt(12)
sigma_y_50V = 50 / 50 / math.sqrt(12)
sigma_y_10V = 10 / 50 / math.sqrt(12)
sigma_x_5ma = 5e-3 / 50 / math.sqrt(12)
sigma_x_500ua = 500e-6 / 50 / math.sqrt(12)



def main():

    # Creo gli oggetti dalla classe per il fit lineare
    monte_10v_5ma = LinearFit()
    monte_50v_5ma = LinearFit()
    valle_10v_500ua = LinearFit()
    valle_50v_5ma = LinearFit()

    # Nel file delle misure x e y sono invertiti (ho prima V e poi I), quindi attivo il flag di swap nella lettura

    # Analisi gruppo 1
    monte_10v_5ma.leggiDati(file1, scale_y=1, scale_x=10**(-3), swap_xy=1) # Le misure di x sono in mA quindi moltiplico tutto per 10^-3
    monte_10v_5ma.add_sigmas(sigmay=sigma_y_10V, sigmax=sigma_x_5ma) # 
    monte_10v_5ma.reg_lin(trasferisci=True)
    monte_10v_5ma.chi_quadro()

    # Analisi gruppo 2
    monte_50v_5ma.leggiDati(file2, scale_y=1, scale_x=10**(-3), swap_xy=1) # Le misure di x sono in mA quindi moltiplico tutto per 10^-3
    monte_50v_5ma.add_sigmas(sigmay=sigma_y_50V, sigmax=sigma_x_5ma)
    monte_50v_5ma.reg_lin(trasferisci=True)
    monte_50v_5ma.chi_quadro()

    # Analisi gruppo 3
    valle_10v_500ua.leggiDati(file3, scale_y=1, scale_x=10**(-6), swap_xy=1) # Le misure di x sono in uA quindi moltiplico tutto per 10^-6
    valle_10v_500ua.add_sigmas(sigmay=sigma_y_10V, sigmax=sigma_x_500ua)
    valle_10v_500ua.reg_lin(trasferisci=True)
    valle_10v_500ua.chi_quadro()

    # Analisi gruppo 4
    valle_50v_5ma.leggiDati(file4, scale_y=1, scale_x=10**(-3), swap_xy=1) # Le misure di x sono in mA quindi moltiplico tutto per 10^-3
    valle_50v_5ma.add_sigmas(sigmay=sigma_y_50V, sigmax=sigma_x_500ua)
    valle_50v_5ma.reg_lin(trasferisci=True)
    valle_50v_5ma.chi_quadro()

    # Quando printo gli oggetti della classe direttamente, viene invocata la funzione __str__ della classe
    print(monte_10v_5ma)
    print(monte_50v_5ma)
    print(valle_10v_500ua)
    print(valle_50v_5ma)

# Non so perchè vada fatta questa cosa ma consigliano di farla hahah
if __name__ == "__main__":
    main()
