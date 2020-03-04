from labbclass import LinearFit, leggiCSV
from matplotlib import pyplot as plt
import numpy as np

bjt = LinearFit()
bjt_lin = LinearFit()
file = "misure.csv"
bjt.Rb = leggiCSV(file, 0)
bjt.sigmaRb = .01*bjt.Rb
bjt.ydata = leggiCSV(file, 1)
bjt_lin.ydata = np.array([i for i in bjt.ydata if i<50e-3])
sigmay_lin = .01* bjt_lin.ydata
sigmay = .01 * bjt.ydata
bjt.xdata = leggiCSV(file, 2)
bjt_lin.xdata = np.array([bjt.xdata[i] for i in range(len(bjt.ydata)) if bjt.ydata[i]<50e-3])
sigmax = 0
sigmax_lin = 0
Rc = 0.1169e3  
bjt.add_sigmas(sigmay=sigmay, sigmax=sigmax) 
bjt_lin.add_sigmas(sigmay=sigmay_lin, sigmax = sigmax_lin)

bjt_lin.reg_lin(trasferisci = False)
bjt_lin.chi_quadro()
print(bjt_lin)

plt.plot(bjt.xdata, bjt.ydata, '.')
bjt_lin.plotData()
plt.show()