from libphysics import *

data_file = "data_diff_sorgente.csv"
file_gcm = "output_common_mode_sorgente.csv"
Rc = 10e3
Re = 120
R1 = 10e3

[freq, Vin, Vout, fase] = readCSV(data_file, skiprows=1)

Vdiff = Vin
Vmedia = Vdiff/2
[Gcm] = readCSV(file_gcm)

Gdiff = Vout/Vdiff - Gcm/2
#Cmrr = abs(Gdiff/Gcm)

#re_stima = Rc/(2*Gdiff) - Re

b1 = bodeplot(freq, Amp=Gcm, Phase=fase)
b2=bodeplot(freq, Amp=[0,2], Phase=fase, figure=b1)
[ax1, ax2] = b2.axes
ax1.legend(["Gcm", "Gdiff"])
ax2.legend(["Gcm", "Gdiff"])
plt.show()
