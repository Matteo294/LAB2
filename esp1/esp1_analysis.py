from labbclass import LinearFit
import os
import math
from matplotlib import pyplot as plt

enable_plots = True # Mettere True per visualizzare i grafici, False per nasconderli

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

    # Creo gli oggetti dalla class e per il fit lineare
    monte_10v_5ma = LinearFit()
    monte_50v_5ma = LinearFit()
    valle_10v_500ua = LinearFit()
    valle_50v_5ma = LinearFit()

    # Nel file delle misure le grandezze x e y sono invertite (ho prima V e poi I), quindi attivo il flag di swap nella lettura

    # Analisi gruppo 1
    monte_10v_5ma.leggiDati(file1, scale_y=1, scale_x=10**(-3), swap_xy=1) # Le misure di x sono in mA quindi moltiplico tutto per 10^-3
    monte_10v_5ma.add_sigmas(sigmay=sigma_y_10V, sigmax=sigma_x_5ma) # 
    monte_10v_5ma.add_resistenza_tester_ICE(4, 3)
    monte_10v_5ma.reg_lin(trasferisci=True)
    monte_10v_5ma.chi_quadro()

    # Analisi gruppo 2
    monte_50v_5ma.leggiDati(file2, scale_y=1, scale_x=10**(-3), swap_xy=1)
    monte_50v_5ma.add_sigmas(sigmay=sigma_y_50V, sigmax=sigma_x_5ma)
    monte_50v_5ma.add_resistenza_tester_ICE(4, 4)
    monte_50v_5ma.reg_lin(trasferisci=True) # Le incertezze di x non sono trascurabili, quindi le trasferisco in y
    monte_50v_5ma.chi_quadro()

    # Analisi gruppo 3
    valle_10v_500ua.leggiDati(file3, scale_y=1, scale_x=10**(-6), swap_xy=1) # Le misure di x sono in uA quindi moltiplico tutto per 10^-6
    valle_10v_500ua.add_sigmas(sigmay=sigma_y_10V, sigmax=sigma_x_500ua)
    valle_10v_500ua.add_resistenza_tester_ICE(5, 3)
    valle_10v_500ua.reg_lin(trasferisci=True)
    valle_10v_500ua.chi_quadro()
    
    # Analisi gruppo 4
    valle_50v_5ma.leggiDati(file4, scale_y=1, scale_x=10**(-3), swap_xy=1)
    valle_50v_5ma.add_sigmas(sigmay=sigma_y_50V, sigmax=sigma_x_5ma)
    valle_50v_5ma.add_resistenza_tester_ICE(4, 4)
    valle_50v_5ma.reg_lin(trasferisci=True)
    valle_50v_5ma.chi_quadro()

    # Quando printo gli oggetti della classe direttamente, viene invocata la funzione __str__ della classe
    print(monte_10v_5ma)
    print(monte_50v_5ma)
    print(valle_10v_500ua)
    print(valle_50v_5ma)
    #print("Resistenza ICE \n5ma Amp = ", monte_10v_5ma.R_amp, "\n10V Vol = ", monte_10v_5ma.R_vol)


    # Correzione valori trovati tenendo conto di resistenza strumenti
    # Configurazione amperometro a monte: R_corretta = (1/R_misurata - 1/R_vol)^-1
    R_corr_monte_10v_5ma = 1/(1/monte_10v_5ma.B - 1/monte_10v_5ma.R_vol)
    R_corr_monte_50v_5ma = 1/(1/monte_50v_5ma.B - 1/monte_50v_5ma.R_vol)
    dR_corr_monte_10v_5ma = monte_10v_5ma.R_vol**2 / (monte_10v_5ma.R_vol - monte_10v_5ma.B)**2 * monte_10v_5ma.sigma_B
    dR_corr_monte_50v_5ma = monte_50v_5ma.R_vol**2 / (monte_50v_5ma.R_vol - monte_50v_5ma.B)**2 * monte_50v_5ma.sigma_B

    # Configurazione amperometro a valle: R_corretta = R_misurata - R_amp
    R_corr_valle_10v_500ua = valle_10v_500ua.B - valle_10v_500ua.R_amp
    R_corr_valle_50v_5ma = valle_50v_5ma.B - valle_50v_5ma.R_amp
    dR_corr_valle_10v_500ua = valle_10v_500ua.sigma_B
    dR_corr_valle_50v_5ma = valle_50v_5ma.sigma_B

    print("Valori della resistenza corretti con la resistenza dell'ICE\n{0:.4f} \t {1:.4f} \t {2:.4f} \t {3:.4f}".format(R_corr_monte_10v_5ma, R_corr_monte_50v_5ma, R_corr_valle_10v_500ua, R_corr_valle_50v_5ma))
    print("Incertezze: {0:.4f} \t {1:.4f} \t {2:.4f} \t {3:.4f}".format(dR_corr_monte_10v_5ma, dR_corr_monte_50v_5ma, dR_corr_valle_10v_500ua, dR_corr_valle_50v_5ma))
    if enable_plots:
        # Grafici  
        '''
            Per cambiare i parametri del grafico dei dati basta usare l'attributo data_plot (oggetto della classe matplotlib).
            Invece per cambiare quelli del grafico della regressione lineare basta usare l'attributo regression_plot
        '''
        monte_10v_5ma.plotData(xlabel='Corrente [mA]', ylabel='Tensione [V]', xscale=10**3,  title='Grafico I-V amperometro a monte con F.S. 50V - 5mA')
        monte_10v_5ma.data_plot.set_markerfacecolor('#ADFF2F') # data_plot è un oggetto matplotlib creato nella funzione plotData all'interno della classe
        monte_10v_5ma.regression_plot.set_color('gray')
        plt.legend(fontsize=18)
        plt.show()

        monte_50v_5ma.plotData(xlabel='Corrente [mA]', ylabel='Tensione [V]', xscale=10**3, title='Grafico I-V amperometro a monte con F.S. 50V - 5mA')
        monte_50v_5ma.data_plot.set_markerfacecolor('#ADFF2F')
        monte_50v_5ma.regression_plot.set_color('gray')
        plt.legend(fontsize=18)
        plt.show()

        valle_10v_500ua.plotData(xlabel='Corrente [mA]', ylabel='Tensione [V]', xscale=10**3, title='Grafico I-V amperometro a valle con F.S. 10V - 500uA')
        valle_10v_500ua.data_plot.set_markerfacecolor('#ADFF2F')
        valle_10v_500ua.regression_plot.set_color('gray')
        plt.legend(fontsize=18)
        plt.show()

        valle_50v_5ma.plotData(xlabel='Corrente [mA]', ylabel='Tensione [V]', xscale=10**3, title='Grafico I-V amperometro a valle con F.S. 50V - 5mA')
        valle_50v_5ma.data_plot.set_markerfacecolor('#ADFF2F')
        valle_50v_5ma.regression_plot.set_color('gray')
        plt.legend(fontsize=18)
        plt.show()

        #------------ Grafici dei residui -----------#
        monte_10v_5ma.residui(xlabel = 'Corrente [mA]', ylabel='Scarto [V]', xscale=10**3)
        monte_10v_5ma.residui_plot.set_markerfacecolor([0, 0.1, 0.6])
        plt.show()

        monte_50v_5ma.residui(xlabel = 'Corrente [mA]', ylabel='ScaScarto [V]', title='Grafico dei residui amperometro a monte con F.S. 50V - 5mA', xscale=10**3)
        monte_50v_5ma.residui_plot.set_markerfacecolor([0, 0.1, 0.6])
        plt.show()

        valle_10v_500ua.residui(xlabel = 'Corrente [mA]', ylabel='Scarto [V]', title='Grafico dei residui amperometro a valle con F.S. 10V - 500uA', xscale=10**3)
        valle_10v_500ua.residui_plot.set_markerfacecolor([0, 0.1, 0.6])
        plt.show()

        valle_50v_5ma.residui(xlabel = 'Corrente [mA]', ylabel='Scarto [V]', title='Grafico dei residui amperometro a valle con F.S. 50V - 5mA', xscale=10**3)
        valle_50v_5ma.residui_plot.set_markerfacecolor([0, 0.1, 0.6])
        plt.show()

if __name__ == "__main__":
    main()
