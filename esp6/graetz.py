from labbclass import Misura
from labbclass import LinearFit
import os
from matplotlib import pyplot as plt
import sys

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

C = Misura(220e-6, 1e-6) # Incertezza arbitraria
R3 = Misura(1.0015e-3, 0.1) # Incertezza sulla lettura
R = Misura(99.653, 0.01) # Incertezza sulla lettura
da_plottare = [] # Grafici da plottare

for grafico in sys.argv[1:]:
    da_plottare.append(int(grafico))




#-------------------------------- Analisi senza Zener -----------------------------------#
soloC = LinearFit()
soloC.leggiDati('Misure/soloC.csv')
soloC.add_sigmas(0, 0)
if 1 in da_plottare:
    plt.plot(soloC.ydata/soloC.xdata, soloC.ydata, '.')
    plt.grid()
    plt.show()
# ---------------------------------------------------------------------------------------#






#---------------------------------- Analisi con Zener -----------------------------------#
zener = LinearFit()
zener.leggiDati('Misure/zener.csv')
zener.add_sigmas(0, 0)
if 2 in da_plottare:
    plt.plot(zener.ydata/zener.xdata, zener.ydata, '.')
    plt.grid()
    plt.show()
# ---------------------------------------------------------------------------------------#