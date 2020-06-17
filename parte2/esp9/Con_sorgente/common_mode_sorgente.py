from libphysics import *

data_file = "data_cm_sorgente.csv"
output_file = "output_common_mode_sorgente.csv"
Rc = 10e3
Re = 120

[freq, Vin, Vout, fase] = readCSV(data_file, skiprows=1)

Vcm = Vin
Gcm = Vout/Vcm
print(Gcm)

# write Gcm in file
with open(output_file, 'w') as f:
    writer = csv.writer(f)
    for element in Gcm:
        writer.writerow([element])
