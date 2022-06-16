#!/bin/bash

#HIC PIPELINE (script che eseguirà juicer, tadbit, hicexplorer ecc)


#scriptsDir (-D): la directory contenente gli SCRIPT di juicer (/opt/applications/scripts/juicer)
#topDir (-d): la directory che conterrà gli output prodotti con juicer (lo script crea "juicer_out" e la rendo cwd. cwd è topDir di default)

#site_file (option -y): il file con i restriction enzyimes sites. Devo averne uno nella apposita dir "restriction_sites" nella juice_

assoc_file=$1  #file di associazione (.tsv)
scriptsDir=$2   #path alla directory contenente gli script di juicer, tadbit ed hicexplorer (script.sh e dir common DEVONO essere nella stessa directory) (ex scriptsDir)
genomeID=$3  #es: hg19 (option -g juicer)
genome_fa=$4 #path al file .fa del genoma di riferimento (option -z juicer) (/opt/genome/human/hg19/index/hg19.fa)
site_file=$5  #path al restriction site file 

threads=$6 #num of threads (option -t juicer)
tadbit_resolution=$7  #tadbit resolution at the moment (17.03.2022) is 1000000
genome_gem=$8  #abs_path ref genome .gem (/opt/genome/human/hg19/index/gem/hg19.gem)
hicexplorer_resolution=$9   #sottoforma di lista di valori divisi da virgola, o di singolo valore. Es: 5000,10000,25000

### NEW ### Introducing OPTARG instead of constant input from command line

# Activate Anaconda Environment (HiC)
#conda activate hic

##Read arguments
usageHelp="Usage: ${0##*/} [-a assoc_file] [-j scriptsDir] [-i genomeID] \n [-f genome_fa]  [-g genome_gem] [-s site_file] \n [-r tadbit_resolution] [-R hicexplorer_resolution] [-t threads] [-h help]"
scriptsDirHelp="* [scriptsDir] is the absolute path of the directory containing the main script, which must be present in the same directory of "scripts" directory, in which must be hicexplorer_hicrep_TIGET.sh, tadbit_TIGET.sh, juicer_TIGET.sh and the "common" directory containing all the scripts called by juicer"
genomeHelp="* [genomeID] must be defined in the script, e.g. \"hg19\" or \"mm10\" \n [genome_fa] is the absolute path of the .fa file of the reference genome. This file should be present in the same directory with the .fa index files. \n [genome_gem] is the absolute path to the .gem file of the reference genome"
siteFileHelp="* [restriction site file]: enter path for restriction site file (locations of\n  restriction sites in genome; can be generated with the script\n  misc/generate_site_positions.py)"
tadbitResolutionHelp="* [tadbit_resolution] is the resolution that will be utilized by tadbit to make plots"
hicexplorerResolutionHelp="* [hicexplorer_resolution] could be a single value or a list of value divided by a comma. hicexplorer will be executed for each sample, for each indicated resolution"
threadsHelp="* [threads]: number of threads when running tadbit full_mapping, juicer BWA alignment, hicexplorer hicFindTADs, hicDetectLoops, hicDifferentialTAD "
helpHelp="* -h: print this help and exit"


printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$scriptsDirHelp"
    echo -e "$genomeHelp"  #(description of genomeID, genome_fa e genome_gem)
    echo -e "$siteFileHelp"
	echo -e "$tadbitResolutionHelp"
	echo -e "$hicexplorerResolutionHelp"
    echo -e "$threadsHelp"
    echo "$helpHelp"
    exit "$1"
     #echo -e "$stageHelp"
   #echo -e "$pathHelp"       (al momento inutile, perchè diamo il genomeID tramite cui retriva i chrom.sizes)
}


while getopts "a:j:i:f:g:s:r:R:h:t:" opt; do
    case $opt in
	a) assoc_file=$OPTARG ;;   #file di associazione (.tsv)
	j) scriptsDir=$OPTARG ;;	   #path alla directory contenente gli script di juicer (script.sh e dir common DEVONO essere nella stessa directory)	(ex juiceDir)
	i) genomeID=$OPTARG ;;		#es: hg19 (option -g juicer)	
	f) genome_fa=$OPTARG ;;	#path al file .fa del genoma di riferimento (option -z juicer) (/opt/genome/human/hg19/index/hg19.fa)	 
	g) genome_gem=$OPTARG ;;     #abs_path ref genome .gem (/opt/genome/human/hg19/index/gem/hg19.gem)
	s) site_file=$OPTARG ;;     #path al restriction site file
	r) tadbit_resolution=$OPTARG ;;  #tadbit resolution at the moment (17.03.2022) is 1000000    
	R) hicexplorer_resolution=$OPTARG ;;
	t) threads=$OPTARG ;;       #num of threads (option -t juicer)
	h) printHelpAndExit 0;;
	
	[?]) printHelpAndExit 1;;
    esac
done


#----------




while IFS=$'\t' read -r sample ID path_dir T_UT T_Type fastq1 fastq2  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le righe
do
  arr=($sample $ID $path_dir $T_UT $T_Type $fastq1 $fastq2) #salvo le variabili in un array per comodità
  
  echo -e "--- Analyzing sample ${arr[0]}\n"
  mkdir ${arr[2]}/${arr[0]} #creo la dir del sample 
  mkdir ${arr[2]}/${arr[0]}/tadbit_results #creo la dir con gli output di tadbit per il sample in analisi
  mkdir ${arr[2]}/${arr[0]}/fastqc_results
  
  #FASTQC analysis on fastq reads
  fastqc ${fastq1} ${fastq2} -o ${arr[2]}/${arr[0]}/fastqc_results
  
  #TADBIT 
  echo -e "--- TADBIT\n"
  #spiegazione: python3 tadbit_finale1.py sample_name abs_path_R1 abs_path_R2 abs_path_tadbit_results tadbit_resolution abs_path_ref_genome.gem abs_path_ref_genome.fa 
  python3 ${scriptsDir}/scripts/tadbit.py ${arr[0]} ${fastq1} ${fastq2} ${arr[2]}/${arr[0]}/tadbit_results ${tadbit_resolution} ${genome_gem} ${genome_fa} ${threads}
   
																												    
  #JUICER
  echo -e "--- JUICER ${arr[0]}\n"
  mkdir ${arr[2]}/${arr[0]}/juicer_results #creo la dir dove salverò gli output di juicer, specifica di un sample
  mkdir ${arr[2]}/${arr[0]}/straw_results
  cd ${arr[2]}/${arr[0]}/juicer_results #cd alla main dir di un sample (ovvero alla directory dove salverò gli output). importante perchè "topDir" in juicer.sh di default è la cwd. Quindi topDir sarà questa directory (e varierà ad ogni iterazione, per ogni sample) 

  # -d è topDir, la directory in cui finiranno gli output. Non la specifico perchè di default è la cwd, specificata prima. -n = sample name
  bash ${scriptsDir}/scripts/juicer.sh -D ${scriptsDir} -p ${genomeID} -z ${genome_fa} -y ${site_file} -n ${arr[0]} -u ${fastq1} -v ${fastq2} -t 8 

  

done < <(tail -n +2 ${assoc_file}) #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)
  
#HICEXPLORER & HICREP - devo runnarlo al di fuori del loop precedente perchè questo script è fatto per funzionare da solo, al suo interno si looppa a sua volta, quindi
#non ha senso loopare uno script che loopa tra samples. Creerei due loop uno dentro l'altro, che in questo caso è un errore.
echo -e "--- STRAW, HICEXPLORER & HICREP\n"

#eseguo hicexplorer (inclusa analisi tad differenziali) ed hicrep una volta per ogni risoluzione desiderata. Ad es, risoluzione 5000, 10000 ecc.
res_arr=() 
IFS=',' read -r -a res_arr <<< "${hicexplorer_resolution}" 
echo "Resolutions choosed are: ${res_arr[@]}" 
 
arr_len="${#res_arr[@]}"  #lunghezza dell array 1 (che è uguale a quella dell'array 2 di solito) 
for (( i=0; i<=${arr_len}-1; i++ ))    #faccio -1 perchè arr_len è il numero di elementi dell'array, ma gli indici iniziano da 0, quindi devo sottrarre 1 al num totale di elementi 
do
echo -e "--- Executing hicexplorer with resolution ${res_arr[i]}"
bash ${scriptsDir}/scripts/hicexplorer_hicrep.sh ${assoc_file} ${res_arr[i]} ${threads}

#STRAW - eseguito per ogni risoluzione indicata nel vettore "hicexplorer_resolution"
echo -e "--- STRAW\n"
python3 ${scriptsDir}/straw.py ${arr[2]}/${arr[0]}/juicer_results/aligned/inter_30.hic ${res_arr[i]} ${arr[2]}/${arr[0]}/straw_results

done




### V4 ###
#prova di esecuzione - introduzione delle optarg, input dei files utilizzando una lettera
#./HIC_PIPE_V4.sh -a /home/alessio/hic/complete_hic_pipe_exe/association_file2.tsv -j /home/alessio/hic/complete_hic_pipe_exe -i hg19 -f /opt/genome/human/hg19/index/hg19.fa -s /home/alessio/hic/complete_hic_pipe_exe/restriction_sites/hg19_DpnII.txt -t 12 -r 1000000 -g /opt/genome/human/hg19/index/gem/hg19.gem -R 1000000,500000



### V3 ###
#prova di esecuzione
#./HIC_PIPE_V3.sh /home/alessio/hic/complete_hic_pipe_exe/association_file2.tsv /home/alessio/hic/complete_hic_pipe_exe hg19 /opt/genome/human/hg19/index/hg19.fa /home/alessio/hic/complete_hic_pipe_exe/restriction_sites/hg19_DpnII.txt 12 1000000 /opt/genome/human/hg19/index/gem/hg19.gem 1000000,500000

#--------------

### V2 ###

#prova di esecuzione
#./HIC_PIPE_V2.sh /home/alessio/hic/complete_hic_pipe_exe/association_file2.tsv /home/alessio/hic/complete_hic_pipe_exe hg19 /opt/genome/human/hg19/index/hg19.fa /home/alessio/hic/complete_hic_pipe_exe/restriction_sites/hg19_DpnII.txt 12 1000000 /opt/genome/human/hg19/index/gem/hg19.gem


#--------------
#V1

#comando tipo per juicer originale
#./juicer.sh -z /opt/genome/human/hg19/index/hg19.fa -p hg19 -D /home/alessio/hic/juicer_prova


#comando runnando juicer_CPU_MOD_OK senza lo script bash superiore
#./juicer_CPU_MOD.sh -D /opt/applications/scripts/juicer -z /opt/genome/human/hg19/index/hg19.fa -p hg19 -t 8 -n 2_1_S1_L001 -u /home/alessio/hic/juicer_prova/fastq/2_1_S1_L001_R1_001.fastq.gz -v /home/alessio/hic/juicer_prova/fastq/2_1_S1_L001_R2_001.fastq.gz


#comando di esempio con cui runnare questo script.
#./HIC_PIPE_V1.sh /home/alessio/hic/juicer_prova/juicer_AF.tsv /home/alessio/hic/juicer_prova hg19 /opt/genome/human/hg19/index/hg19.fa /home/alessio/hic/juicer_prova/restriction_sites/hg19_DpnII.txt 12




