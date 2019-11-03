import os

def latex_writer(f, nfile='noname'):

    rows = len(f.splitlines())
    columns = len(f.splitlines()[0].split(','))
    print(str(rows) + ' ' + str(columns))
    fname = 'results/' + nfile + '.txt'
    filename = os.path.join(path, fname)
    latexfile = open(filename, 'w')

    ############# Testo prima della tabella ####################
    textbefore = chr(92) + 'begin' + chr(123) + 'table' + chr(125) + '[h]\n\t' +  chr(92) + 'centering\n\t'
    textbefore = textbefore + chr(92) + 'begin' + chr(123) + 'tabular' + chr(125) + '{'
    
    for i in range(columns):
        textbefore = textbefore + '|c'
    textbefore = textbefore + '|}' + '\n'
    latexfile.write(textbefore)
    ###########################################################

    for row in f.splitlines():
        latexline = '\t\t' + row.replace(',', '\t' + chr(38) + '\t') + '\t' + chr(92) + chr(92) + '\n'
        latexfile.write(latexline)  

    ################ Testo dopo la tabella ####################
    textafter = '\t' + chr(92) + 'end' + chr(123) + 'tabular' + chr(125) + '\n'
    textafter = textafter + chr(92) + 'end' + chr(123) + 'table' + chr(125)
    latexfile.write(textafter)
    ###########################################################

    latexfile.close() 

# Mi metto nella cartella di questo file (non è sempre detto che lo faccia di default)
path = os.path.dirname(os.path.realpath(__file__))
print('Mi trovo in ' + path)
# Lista dei file in questa cartella
files = os.listdir(path)
# Rimuovo dalla lista questo file
files.remove('main.py')
files.remove('info.txt')
if 'results' in files:
    files.remove('results')
print(files)

print('\n')
table = [] # inizializzo l'array per contenere i dati
for f in files:
    # Combino questa directory con ciascun file
    filepath = os.path.join(path, f)
    # Provo a leggere il file, altrimenti stampo un errore
    try:
        ffile = open(filepath, 'r')
        data = ffile.read()
        print(f)
        print(data)
        print('\n')
        # chiamo la funzione definita a inizio programma
        print(f[:-4])
        latex_writer(data, f[:-4])
        ffile.close()
    except:
        print('Non sono riuscito a leggere il file ' + f + '\n')

input('Premere un tasto per terminare')

    
