import csv
import os
import numpy as np 

cartella = '../sistemati/' # cartella dove mettere i sistemati

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

rawDataDir = dname + "/grezzi"
os.chdir(rawDataDir)

onlyfiles = [f for f in os.listdir() if os.path.isfile(os.path.join(f))]
print(len(onlyfiles))

for f in onlyfiles:
    nomefile, estensione = f.split('.')
    if estensione == 'csv':

        file_lettura = f
        file_scrittura = cartella + f

        #myfile = os.path.join(file_lettura)
        if os.path.isfile(file_lettura):
            with open(file_lettura, 'r') as csvFile:
                reader = csv.reader(csvFile)
                i = 0
                xdata = np.asarray([])
                ydata = np.asarray([])
                for r in reader:
                    if i >= 2:
                        if r[0] and r[1]: # Controllo che non siano vuoti
                            xdata = np.append(xdata, r[0])
                            ydata = np.append(ydata, r[1])
                    i += 1

            
            with open(file_scrittura, 'w') as filew:
                writer = csv.writer(filew, delimiter=',')
                for x,y in zip(xdata, ydata):
                    writer.writerow([x, y])
