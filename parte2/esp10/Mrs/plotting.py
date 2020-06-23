from matplotlib import pyplot as plt 
from libphysics import *

plt.rc('text', usetex=True) 

d, val, approx = readCSV('induzione.csv')
plt.plot(d, val, label='Formula di Neumann')
plt.plot(d, approx, label='Approssimazione dipolo')
plt.xlabel('Distanza [m]')
plt.ylabel(r'Induzione mutua [$\mu$H]')
plt.legend()
plt.show()