#TODO - fornire in input anche l'enzima di restrizione (variabile "enzymes")

#argomenti da fornire in input:
#1) sample name
#2) path assoluto r1
#3) path assoluto r2
#4) path directory in cui finiranno gli output
#5) risoluzione per il plot della matrice hic in "hic_map" function
#6) genoma index gem (hg19.gem)
#7) genoma normale (hg19.fa)

#prova d'uso: 
#python3 tadbit_automatiz2.py 2_1_S1_L001 /mnt/externalhd2/giulio/HiC/hic13/tadbit_quality_prova/fastq/2_1_S1_L001_R1_001.fastq.gz /mnt/externalhd2/giulio/HiC/hic13/tadbit_quality_prova/fastq/2_1_S1_L001_R2_001.fastq.gz /mnt/externalhd2/giulio/HiC/hic13/tadbit_quality_prova/results 1000000 /opt/genome/human/hg19/index/gem/hg19.gem /opt/genome/human/hg19/index/hg19.fa


import pytadbit
from pytadbit.mapping.mapper import full_mapping
import os
import sys
#import re
#from re import search

from pytadbit.utils.fastq_utils import quality_plot
from pytadbit.parsers.genome_parser import parse_fasta
from pytadbit.parsers.map_parser import parse_map
from pytadbit.mapping import get_intersection

from pytadbit.mapping.analyze import plot_distance_vs_interactions
from pytadbit.mapping.analyze import plot_genomic_distribution
from pytadbit.mapping.analyze import hic_map
from pytadbit.mapping.analyze import insert_sizes
from pytadbit.mapping.analyze import plot_strand_bias_by_distance

# Get the current working directory -> os.getcwd()

# Print the current working directory
#print("Current working directory: {0}".format(os.getcwd()))

#salvo nelle variabili gli argomenti dati in input
sample = sys.argv[1]   #dato che i nomi dei file fastq sono spesso diversi, è difficile capire come individuare il nome del sample. Per il momento propongo di darglielo direttamente in input
r1_path = sys.argv[2] #abs path di r1 di un sample
r2_path = sys.argv[3] #abs path di r2 di un sample
out_dir = sys.argv[4]  #directory dei results (al cui interno avrò le 3 subdirectories citate nei prerequisiti)
resolution = int(sys.argv[5])  #1000000   #ciò che do in input come argomento è sempre una stringa. quindi la trasformo in integer
res = str(resolution)
genome_gem = sys.argv[6]
genome = sys.argv[7]
threads = sys.argv[8]

#salvo nelle variabili i path delle directory dove verranno salvati gli output di tadbit
mapped_reads_path = os.path.join(out_dir, 'mapped_reads')					#dir in cui salverà i file full e frag prodotti dal mapping con gem
uniquely_mapped_reads_path = os.path.join(out_dir, 'uniquely_mapped_reads') #dir in cui salverà le uniquely mapped reads estratte con "parse_map" e poi il file .tsv con le suddette mergiate da get_intersection
plots_path = os.path.join(out_dir, 'qc_hic_experiment_plots')				#dir in cui salverà tutti i plot specifici di tadbit sul qc di un hic experiment
fastq_quality_plots_path = os.path.join(out_dir, 'fastq_quality_plots')		#dir in cui salverò i plot sulla qualità dei fastq				

#Crea le suddette cartelle se non sono già state create.
new_dirs = [mapped_reads_path, uniquely_mapped_reads_path, plots_path, fastq_quality_plots_path]

for d in new_dirs:
	isExist = os.path.exists(d) # Check whether the specified path exists or not
	if not isExist:
		os.makedirs(d) # Create a new directory because it does not exist
		print("The directory", d, "is created!")
	else:
		print("The directory", d, "already exist", "\n")


enzymes = ['DpnII']


#Tadbit reads fastq quality plot
print("--- 	fastq quality plots", "\n")
quality_plot(r1_path, r_enz=enzymes, nreads=1000000, savefig=fastq_quality_plots_path+'/'+sample+'_R1.png')
quality_plot(r2_path, r_enz=enzymes, nreads=1000000, savefig=fastq_quality_plots_path+'/'+sample+'_R2.png')


#### 1) MAP of fastq reads on the reference genome #### -- PROVARE con frag_map=False (ovvero effettuando un iterative mapping e non un fragmented mapping, per vedere cosa cambia)

print("--- 	Mapping reads of sample:", sample, "on the reference genome", "\n")
#le variabili "mapped_r1" e "mapped_r2" sono delle liste contenenti i path assoluti dei file frag e full, rispettivamente per R1 ed R2
mapped_r1 = full_mapping(mapper_index_path=genome_gem, fastq_path=r1_path,  out_map_dir=mapped_reads_path, windows=(1, 151), r_enz=enzymes, frag_map=True, nthreads=threads, clean=True, temp_dir=out_dir+'/temp_r1')
mapped_r2 = full_mapping(mapper_index_path=genome_gem, fastq_path=r2_path,  out_map_dir=mapped_reads_path, windows=(1, 151), r_enz=enzymes, frag_map=True, nthreads=threads, clean=True, temp_dir=out_dir+'/temp_r2')


#### 2) PARSE - extraction of uniquely mapped reads + 3) Merge files R1 and R2 with uniquely mapped reads in one unique file - reads found in both files are merged ####

genome_seq = parse_fasta(genome)


reads1 = uniquely_mapped_reads_path+"/"+sample+"_R1.tsv"	       #path del file output contenente le uniquely mapped R1 del sample in questione
reads2 = uniquely_mapped_reads_path+"/"+sample+"_R2.tsv"		   #path del file output contenente le uniquely mapped R1 del sample in questione
parse_map(mapped_r1, mapped_r2, out_file1=reads1, out_file2=reads2, genome_seq=genome_seq, re_name='DpnII', verbose=True)
	
merged_reads = uniquely_mapped_reads_path+"/"+sample+"_merged_reads.tsv"   #path del file che conterrà le reads mergiate (ex reads12.tsv)
get_intersection(reads1, reads2, merged_reads, verbose=True)

#### 4) Quality check of the Hi-C experiment ####
# - plot_distance_vs_interactions - andamento delle interazioni tra due regioni genomiche in base alla loro distanza (non escono i numeri che si vedono sul tutorial)
print("--- Plotting distance vs interactions", "\n")
ddi = plot_distance_vs_interactions(merged_reads, resolution=100, max_diff=1000, show=True, savefig=plots_path+'/'+sample+'_plot_distance_vs_interactions.png')  #parametri di default

    
# - plot genomic distribution - n° of reads per bin (coverage)
print("--- Plotting genomic distribution - n° of reads per bin", "\n")
plot_genomic_distribution(merged_reads, resolution=100, ylim=None, show=True, savefig=plots_path+'/'+sample+'_plot_genomic_distribution.png') #nel tutorial era resolution=500000, ylim=(0, 100000)
#plt.tight_layout() #non funziona se non inserisco anche la parte del "decay by cromosome"
    
# - plot hic map
print("--- Plotting hic map", "\n")
hic_map(merged_reads, resolution=resolution, show=True, cmap='viridis', savefig=plots_path+'/'+sample+'_hic_map_res_'+res+'.png')  #resolution = 1000000
	
#TROVARE IL N° DI READS DA MAPPED_READS, PER SAMPLE
#nreads=$((`zcat ${R1_FASTQ} | wc -l`/4)) ;
	
#From the reads that are mapped in a single RE fragment (dangling-end reads) we can infer the average insert size:
print("--- Plotting distirbution of dangling ends length", "\n")
sizes = insert_sizes(merged_reads, show=True, nreads=None, savefig=plots_path+'/'+sample+'_distribution_of_dangling_ends_lengths.png')  #il n° di reads dovrebbe essere contenuto nel file "reads12.tsv" (ora da me chiamato "merged_reads") (prima usato 417838) il plot viene tuttavia viene identico sia con none che con quel valore. come fa la funzione ad identificare solo le dangling?)

#This function separates each read-end pair into 4 categories depending of the orientation of the strand in which each maps.
print("--- Plotting strand bias by distance", "\n")
plot_strand_bias_by_distance(merged_reads, valid_pairs=False, full_step=None, nreads=None, savefig=plots_path+'/'+sample+'_plot_strand_bias_by_distance.png')

#Capire quale sia la distanza minima tra le read-ends alla quale la loro distribuzione nelle 4 categorie è uniforme 
plot_strand_bias_by_distance(merged_reads, nreads=None, savefig=plots_path+'/'+sample+'_plot_strand_bias_by_distance2.png')


#remove the "uniquely_mapped_reads" directory to save space (we don't need it anymore)
os.rmdir(uniquely_mapped_reads_path)




