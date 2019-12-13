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
da_plottare = []


for grafico in sys,argv[1:]:





#-------------------------------- Analisi senza Zener -----------------------------------#
soloC = LinearFit()
soloC.leggiDati('Misure/soloC.csv')
soloC.add_sigmas(0, 0)
plt.plot(1/soloC.xdata, soloC.ydata, '.')
plt.show()
# ---------------------------------------------------------------------------------------#






#---------------------------------- Analisi con Zener -----------------------------------#
zener = LinearFit()
zener.leggiDati('Misure/zener.csv')
zener.add_sigmas(0, 0)
plt.plot(1/zener.xdata, zener.ydata, '.')
plt.show()
# ---------------------------------------------------------------------------------------#