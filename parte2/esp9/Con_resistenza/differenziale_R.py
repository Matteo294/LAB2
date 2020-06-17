from libphysics import *

data_file = "data_diff_R.csv"
file_gcm = "output_common_mode.csv"
Rc = 10e3
Re = 120
R1 = 10e3

[freq, Vin, Vout, fase] = readCSV(data_file, skiprows=1)

Vdiff = Vin
Vmedia = Vdiff/2
[Gcm] = readCSV(file_gcm)

Gdiff = Vout/Vdiff - Gcm/2
Cmrr = abs(Gdiff/Gcm)

re_stima = Rc/(2*Gdiff) - Re

bodeplot(freq, Amp=Gcm, Phase=fase)
#plt.show()
