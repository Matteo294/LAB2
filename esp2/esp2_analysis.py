from labbclass import LinearFit
import numpy as np
from matplotlib import pyplot as plt
import os
import math

# cd nella directory di questo file (non sempre ci troviamo qui in automatico)
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# PROVA LETTURA FILE
scariche = np.ndarray((5,5), dtype=LinearFit)     # matrice: in una colonna hai i 5 scariche fatte con la stessa resistenza; colonne diverse per resistenze diverse
partial_path = 'Misure/Ctot/sistemati/r'
for i in range(5):
    suffisso_1 = str(i+1) + '_'
    min = 0
    epsilon = 1e-10
    for j in range(5):
        suffisso_2 = str(j+1) + '.csv'
        suffisso_2_txt = str(j+1) + '.txt'
        total_path = partial_path + suffisso_1 + suffisso_2
        total_path_txt = partial_path + suffisso_1 + suffisso_2_txt
        scariche[i,j] = LinearFit()
        scariche[i,j].leggiDati(total_path)
    #traslo in su se trovo valori negativi
        for value in scariche[i,j].ydata:
            if value < min:
                min = value
        scariche[i,j].ydata += -min + epsilon
    #le righe che seguono servono a beccare dal txt la risoluzione dell'oscilloscopio
        with open(total_path_txt, 'r') as filetxt:
            lines = filetxt.readlines()
            numero = lines[1][11:14]
            ris = numero[0] + numero[1] + numero[2]
            ris = int(ris)*1e-3
        scariche[i,j].add_sigmas(sigmay=ris)


# REGRESSIONE PER TROVARE scariche SULLE SINGOLE SCARICHE
for i in range(5):
    for j in range(5):
        scariche[i,j].reg_lin(cambiaVariabili=True, y=np.log(scariche[i,j].ydata), x=scariche[i,j].xdata, sigma=scariche[i,j].sigmay/scariche[i,j].ydata, trasferisci=True)

'''
tau_R = LinearFit()       # per fare la regressione tra 1/tau e 1/R
R = np.array([1,1,1,1,1]) # inserire valori resistenze usate davvero
sigmaR = 1
for i in range(5):      # per ogni colonna i, calcolo la media delle B (B=1/tau)
    B_colonna = np.empty(5)   # vettore delle B per ciascun elemento della colonna
    for j in range(5):         # sulla colonna, calcolo la media delle B
            B_colonna[j] = -scariche[j, i].B
    sigma_B = np.std(B_colonna)/(5)
    tau_R.ydata = np.append(tau_R.ydata, np.mean(B_colonna))
    tau_R.xdata = np.append(tau_R.xdata, 1/R[i])
    tau_R.sigmay = np.append(tau_R.sigmay, np.std(B_colonna)/(5))
    print(np.mean(B_colonna))
    tau_R.sigmax = np.append(tau_R.sigmax, 1/R[i]*sigmaR)        # non sono sicura del /n, o /(n-1), controllare

print(tau_R.sigmay)
tau_R.reg_lin(trasferisci=True)
print(1/tau_R.B)
'''
'''  
for tau in scariche:
    for r in tau:
        index_to_del = np.where(r.ydata <= 0)
        r.xdata = np.delete(r.xdata, index_to_del)
        r.ydata = np.delete(r.ydata, index_to_del)
        r.sigmay = np.delete(r.sigmax, index_to_del)
        r.sigmax = np.delete(r.sigmax, index_to_del)
'''

# Risolvo con il metodo di Cramer la soluzione alle equazioni che massimizzano la likelihood (o minimizzano il chi2)
# d(chi2)/da = 0, dove a è uno dei coefficienti da ricavare
for tau in scariche:
    print('\n\n')
    for r in tau:

        # y = A + Bt + Cx   con x = 1/Vm  e  y = log(Vm)    
        x = 1/r.ydata
        t = r.xdata
        y = np.log(r.ydata)

        w = 1/r.sigmay**2

        M = np.asarray([
            [np.sum(w), np.sum(w*t), np.sum(w*x)],
            [np.sum(w*t), np.sum(w*(t**2)), np.sum(w*x*t)],
            [np.sum(w*x), np.sum(w*x*t), np.sum(w*(x**2))]
        ])

        delta = np.linalg.det(M)

        MA = np.copy(M)
        MB = np.copy(M)
        MC = np.copy(M)
        MA[:,0] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*x*y)]
        MB[:,1] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*y*x)]
        MC[:,2] = [np.sum(w*y), np.sum(w*y*t), np.sum(w*y*x)]

        A = np.linalg.det(MA) / delta
        B = np.linalg.det(MB) / delta
        C = np.linalg.det(MC) / delta

        #print(M)
        print(MA == M)

        print(delta)
        print("Risultati: A = {0:.8f} \t B = {1:.8f} \t C = {2:.8f}".format(A, B, C))


 
# print(1/scariche[0,0].B, scariche[0,0].sigma_B)
'''n_ripetute = 5 # numero di set per ogni misura
scariche_tot = np.asarray([])

for scariche in scariche_tot:
    set_misure = [LinearFit() for i in range(n_ripetute)] # Array di set: per ogni resistenza (per ogni scariche) ho 5 set
    scariche = set_misure # scariche è un array di oggetti della LinearFit: sarà il risultato dell'analisi di questi set'''