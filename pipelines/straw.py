#DA ESEGUIRE DOPO JUICER



#fonti: https://github.com/aidenlab/straw/wiki/Python
#https://github.com/aidenlab/straw/tree/master/pybind11_python

#sembrano due versioni diverse, la seconda più recente. Ma la prima sembra spiegare leggermente più chiaramente i vari passaggi

#ESECUZIONE: python3 straw.1.py /home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic 1000000 /home/alessio/hic/complete_hic_pipe_exe/straw
#ESECUZIONE E22: python3 straw.1.py /mnt/externalhd2/giulio/HiC/hic10/E22/inter_30.hic 1000000 /home/alessio/hic/complete_hic_pipe_exe/straw


import hicstraw

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.pyplot as plt 

import sys
import os



hic_map = sys.argv[1]  #'/home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic'
resolution = sys.argv[2] #1000000
out_dir = sys.argv[3] #/home/alessio/hic/complete_hic_pipe_exe/straw

#making output dir
counts_plots = out_dir+'/'+resolution+'/counts_plots/'
chr_counts = out_dir+'/'+resolution+'/chr_counts/'

os.mkdir(counts_plots)
os.mkdir(chr_counts)



#file di prova fornito nell'usage:
#hic_map = hicstraw.HiCFile("https://www.encodeproject.org/files/ENCFF718AWL/@@download/ENCFF718AWL.hic")


#HiCFile class object methods. First, we need to create the object giving the .hic matrix to the method
hic = hicstraw.HiCFile(hic_map)
print(hic.getGenomeID())    #ottengo una stringa con il genome ID
print(hic.getResolutions()) #ottengo una lista di digits (tutte le risoluzioni contenute nel file .hic)


#individuo quali cromosomi sono presenti nella matrice hic:
chrom_names = []

#"getChromosomes" method generate a "chromosome" object for each chromosomes HiCFile class file (here, "hic")
#class "chromosomes" objects are characterized by "index", "name", "length" data descriptors 
for chrom in hic.getChromosomes():  
  #print(chrom.name, chrom.length, chrom.index);
  if chrom.name not in ['ALL', 'MT']:  #escludo la dicitura ALL e il chr mitocondriale
  	chrom_names.append(chrom.name)
print(chrom_names)


#normalizations utilized by juicer: VC, VC_SQRT, KR,SCALE, GW_KR, GW_SCALE, GW_VC, INTER_KR, INTER_SCALE, INTER_VC [default: VC,VC_SQRT,KR,SCALE]
normalization_types = ['NONE', 'VC', 'VC_SQRT', 'KR'] #juicer default normalizations (ho rimosso 'SCALE' perchè nella nostra matrice hic delle shallow, e nelle E22 non c'era (File did not contain SCALE normalization vectors for one or both chromosomes at 1000000 BP)). Idem per le altre normalizzazioni non incluse nella lista delle default



REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])

#helper function for plotting numpy.ndarray object generated with .strawAsMatrix
def plot_hic_map(dense_matrix, maxcolor):
    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
    plt.savefig(out_dir+'/'+resolution+'/counts_plots/'+'observed_'+chrom+'_'+norm)
    #plt.show()

#--
#a volte si blocca dopo aver elaborato il chr1 con il seguente errore:
#Traceback (most recent call last):
#  File "straw.1.py", line 74, in <module>
#    result_straw = hicstraw.straw('observed', str(norm), str(hic_map), str(chrom), str(chrom), 'BP', int(resolution))
#ValueError: cannot create std::vector larger than max_size()

#altre volte il suddetto errore appare all' observed NONE 5, e la Mem raggiunge i 252gb e lo swap a 5gb, CPU
#arriva ad observed SCALE 6 e lo swap raggiunge il limite di 8gb e il processo viene killato.
#tutti gli output fino al chr6 incluso vengono prodotti, tranne SCALE del 6, sia .tsv che i plot .png

#--




#.tsv con i counts e plot della matrice con i counts dei soli singoli cromosomi interi (1, 2, 3 ecc)
# !! ipotesi introduzione delle coordinate da input in quanto i campi chr1 e chr2 possono essere compilati anche cosi: '4:1000000:2000000', '4:1000000:2000000' !!
for chrom in chrom_names:
	for norm in normalization_types:
		print("observed", norm, chrom)
		result_straw = hicstraw.straw('observed', str(norm), str(hic_map), str(chrom), str(chrom), 'BP', int(resolution))
		result_strawAsMatrix = hicstraw.strawAsMatrix('observed', str(norm), str(hic_map), str(chrom), str(chrom), 'BP', int(resolution))
		with open(out_dir+'/'+resolution+'/chr_counts/'+'chr'+chrom+'_'+norm+'.tsv', 'w', newline='\n') as tsvfile: 
			for i in range(len(result_straw)):
				tsvfile.write("%s %s %s\n" % (result_straw[i].binX, result_straw[i].binY, result_straw[i].counts))
		plot_hic_map(result_strawAsMatrix, 30) #30 è il valore utilizzato nell'usage della funzione
	

	
		
	
#-- EXTRA --

#getMatrixZoomData + getRecordsAsMatrix dovrebbe essere la stessa cosa che eseguire .strawAsMatrix...



#This HiCFile object retains information for fast future queries. Here's an example that pick the counts from the intrachromosomal region for chr4 with KR normalization at 5kB resolution.
#We can do it using an HiCFile class specific method ("getMatrixZoomData") that generate a MatrixZoomData class object through an HiCFile object.
#mzd = hic.getMatrixZoomData(chrom1, chrom2, data_type, normalization, "BP", resolution)
#mzd = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000) #in questo caso vogliamo i counts esclusivamente del chr4, e quindi per i campi chrom1 e chrom2 indichiamo sempre '4'. Se avessimo voluto i counts per una regione più grande, ad es da chr4 a chr6, avremmo dovuto indicare '4' e '6'. Credo appunto che questa sintassi mi faccia ottenere i counts dal chr4 al chr6 (quindi includendo anche il ch5), e non solo del chr4 e del chr6 
#mzd2 = hic2.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000)

### GETRECORDASMATRIX ### - works

#possiamo interrogare l'oggetto "mzd" tramite gli specifici metodi degli oggetti di classe "MatrixZoomData":
#numpy_matrix_chr4 = mzd.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000) #crea un oggetto di classe "numpy.ndarray", plottabile
#numpy_matrix_chr4_2 = mzd2.getRecordsAsMatrix(10000000, 12000000, 10000000, 12000000)

#print("plottiamo numpy_matrix_chr4 e salviamo il plot")

#plottiamo numpy_matrix_chr4 e salviamo il plot (nostra matrice hic)
#REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
#def plot_hic_map(dense_matrix, maxcolor):
#    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
#    plt.savefig("numpy_matrix_chr4_getRecordAsMatrix")
    #plt.show()

#plot_hic_map(numpy_matrix_chr4, 30)
#print("plot OK")

#plottiamo numpy_matrix_chr4_2 e salviamo il plot (matrice hic dell'usage) (ridefinisco la funzione perchè devo ancora trovare un modo per inserire il nome della varibile in "plt.savefig()")
#REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
#def plot_hic_map(dense_matrix, maxcolor):
#    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
#    plt.savefig("numpy_matrix_chr4_getRecordAsMatrix_2")
    #plt.show()
#plot_hic_map(numpy_matrix_chr4_2, 30)


### GETEXPECTEDVALUES ### - not working


#print("plottiamo numpy_matrix_chr4_expected_values e salviamo il plot")
#interroghiamo l'oggetto mzd tramite il metodo "getExpectedValues"
#numpy_matrix_chr4_expected_values = mzd.getExpectedValues()

#not working:
#plottiamo numpy_matrix_chr4_expected_values e salviamo il plot
#REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
#def plot_hic_map(dense_matrix, maxcolor):
#    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
#    plt.savefig("numpy_matrix_chr4_getExpectedValues")
    #plt.show()

#plot_hic_map(numpy_matrix_chr4_expected_values, 30)
#print("plot OK")

### GETNORMVECTOR ### - not working

#l'argomento 0 da fornire a getNormVector è un integer. dev'essere uno dei due index di cromosomi dati al momento della creazione dell'oggetto di classe "MatrixZoomData" (mzd = hic.getMatrixZoomData('4', '4', "observed", "KR", "BP", 5000) quindi nel nostro caso o 4 o 4. Se ci fossero stati '4', '5', avrei potuto mettere o 4 o 5. Fonte: se fornisco un integer errato al metodo getNormVector (ad esempio 10) mi dice "Invalid index provided: 10. Should be either 4 or 4"
#numpy_matrix_chr4_norm_vector = mzd.getNormVector(4) #numpy_matrix_chr4_norm_vector è un oggetto di classe "numpy.ndarray"

#not working
#print("plottiamo numpy_matrix_chr4_norm_vector")
#REDMAP = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
# helper function for plotting
#def plot_hic_map(dense_matrix, maxcolor):
#    plt.matshow(dense_matrix, cmap=REDMAP, vmin=0, vmax=maxcolor)
#    plt.savefig("numpy_matrix_chr4_getNormVector")
    #plt.show()

#plot_hic_map(numpy_matrix_chr4_norm_vector, 30)
#print("plot OK")






### straw ### - documentation
#hicstraw.straw(data_type, normalization, file, region_x, region_y, 'BP', resolution)

#Extract all reads on chromosome X at 1MB resolution with no normalization ("NONE") in file "hic".
#detto in un'altra maniera: to fetch a list of all the raw (definiti come "observed") contacts on chrX at 1MB resolution:
#result = hicstraw.straw('observed', 'NONE', '/home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic', 'X', 'X', 'BP', 1000000)
#print("observed, no normalization")
#for i in range(len(result)):
#    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))




#To fetch a list of KR normalized contacts (reads)  for the same region:
#per utilizzare il modulo ".straw", il file hic non dev'essere importato con hicstraw.HiCFile, basta il suo path assoluto o una variabile in cui esso è salvato
#result = hicstraw.straw('observed', 'KR', '/home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic', 'X', 'X', 'BP', 1000000)
#print("observed, KR normalization")
#for i in range(len(result)):
#    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))


#To query observed/expected KR normalized data (quindi non più i raw, gli observed, ma i count sottoforma di ratio observed/expected, calcolati non si sa come)
#result = hicstraw.straw('oe', 'KR', '/home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic', 'X', 'X', 'BP', 1000000)
#print("observed/expected, KR normalization")
#for i in range(len(result)):
#    print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))



#Stessa cosa, ma con un'altra regione e risoluzione:
#extracts all reads from chromosome 4 at 50KB resolution with KR (Knight-Ruiz balancing algorithm) normalization from the combined MAPQ 30 map from Rao and Huntley et al. 2014
#result = hicstraw.straw("observed", 'KR', '/home/alessio/hic/hic_pipe_results/2_1_S1_L001/juicer_results/aligned/inter_30.hic', '4:1000000:2000000', '4:1000000:2000000', 'BP', 50000)
#print("observed, KR normalization, chr4")
#for i in range(10):
#  print("{0}\t{1}\t{2}".format(result[i].binX, result[i].binY, result[i].counts))




### strawC ### - mi da lo stesso output di .straw


#possibilità di eseguire la funzione .straw, e quindi di ottenere i counts, su più cromosomi (es: da chr2 a chr4)

