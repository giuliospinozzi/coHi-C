#!/bin/bash

#HIC PIPELINE (script che eseguirà juicer, tadbit, hicexplorer ecc)


#scriptsDir (-D): la directory contenente gli SCRIPT di juicer (/opt/applications/scripts/juicer)
#topDir (-d): la directory che conterrà gli output prodotti con juicer (lo script crea "juicer_out" e la rendo cwd. cwd è topDir di default)




# Activate Anaconda Environment (HiC)
#conda activate hic

##Read arguments
usageHelp="Usage: ${0##*/} [-a assoc_file] [-D scriptsDir] [-g genomeID] \n [-z genome_fa]  [-m genome_gem] [-s site] \n [-r tadbit_resolution] [-R hicexplorer_resolution] [-t threads] [-h help]"
scriptsDirHelp="* [scriptsDir] is the absolute path of the directory containing the main script, which must be present in the same directory of "scripts" directory, in which must be hicexplorer_hicrep_TIGET.sh, tadbit_TIGET.sh, juicer_TIGET.sh and the "common" directory containing all the scripts called by juicer"
genomeHelp="* [genomeID] e.g. \"hg19\" or \"mm10\" \n [genome_fa] is the absolute path of the .fa file of the reference genome. This file should be present in the same directory with the .fa index files. \n [genome_gem] is the absolute path to the .gem file of the reference genome"
siteHelp="* [restriction enzyme]: enter the name of the restriction enzyme" 
tadbitResolutionHelp="* [tadbit_resolution] is the resolution that will be utilized by tadbit to make plots"
hicexplorerResolutionHelp="* [hicexplorer_resolution] could be a single value or a list of value divided by a comma. hicexplorer will be executed for each sample, for each indicated resolution"
threadsHelp="* [threads]: number of threads when running tadbit full_mapping, juicer BWA alignment, hicexplorer hicFindTADs, hicDetectLoops, hicDifferentialTAD "
helpHelp="* -h: print this help and exit"


printHelpAndExit() {
    echo -e "$usageHelp"
    echo -e "$scriptsDirHelp"
    echo -e "$genomeHelp"  #(description of genomeID, genome_fa e genome_gem)
    echo -e "$siteHelp"
	echo -e "$tadbitResolutionHelp"
	echo -e "$hicexplorerResolutionHelp"
    echo -e "$threadsHelp"
    echo "$helpHelp"
    exit "$1"
     #echo -e "$stageHelp"
   #echo -e "$pathHelp"       (al momento inutile, perchè diamo il genomeID tramite cui retriva i chrom.sizes)
}


while getopts "a:D:g:z:m:s:p:r:R:t:l:y:h:" opt; do
    case $opt in
	a) assoc_file=$OPTARG ;;   #file di associazione (.tsv)
	D) scriptsDir=$OPTARG ;;	   #path alla directory contenente gli script di juicer (script.sh e dir common DEVONO essere nella stessa directory)	(ex juiceDir)
	g) genomeID=$OPTARG ;;		#es: hg19 (option -g juicer)	
	z) genome_fa=$OPTARG ;;	#path al file .fa del genoma di riferimento (option -z juicer) (/opt/genome/human/hg19/index/hg19.fa)	 
	m) genome_gem=$OPTARG ;;     #abs_path ref genome .gem (/opt/genome/human/hg19/index/gem/hg19.gem)
	s) site=$OPTARG ;;     #Restriction enzyme (es: DpnII)
	p) genomePath=$OPTARG ;; #path to the chrom.size file
	r) tadbit_resolution=$OPTARG ;;  #tadbit resolution at the moment (17.03.2022) is 1000000    
	R) hicexplorer_resolution=$OPTARG ;;
	t) threads=$OPTARG ;;       #num of threads (option -t juicer)
	l) is_shallow=$OPTARG ;;    #true or false. se true, allora tadbit verrà eseguito per intero. Se false, verrà eseguito solo il quality plot. Questo perchè è molto oneroso sui fastq enormi. Se non true e non false, raise error
        y) site_file=$OPTARG ;; #restriction site file absolute path, for juicer execution
	h) printHelpAndExit 0;;
	
	[?]) printHelpAndExit 1;;
    esac
done


#----------

### Assembly-stats, tadbit and juicer for each sample included in association file ###
while IFS=$'\t' read -r sample ID path_dir T_UT T_Type fastq1 fastq2  #saving each values of row i in some corresponding variables (row i = sample i)
do
  arr=($sample $ID $path_dir $T_UT $T_Type $fastq1 $fastq2) #saving those variables in a single array 
  
  echo -e "--- Analyzing sample ${arr[0]}\n"
  mkdir ${arr[2]}/${arr[0]} #creating sample dir 
  mkdir ${arr[2]}/${arr[0]}/tadbit_results #creating dir that will contain tadbit outputs for sample i
  
  #ASSEMBLY-STATS
  echo -e "--- Assembly-stats on R1 & R2 of sample ${arr[0]}\n"
  assembly-stats <(zcat ${fastq1}) &
  assembly-stats <(zcat ${fastq2}) &
  wait

  #TADBIT 
  echo -e "--- TADBIT\n"
  #spiegazione: python3 tadbit.py sample_name abs_path_R1 abs_path_R2 abs_path_tadbit_results tadbit_resolution abs_path_ref_genome.gem abs_path_ref_genome.fa threads true/false
  python3 ${scriptsDir}/tadbit.py ${arr[0]} ${fastq1} ${fastq2} ${arr[2]}/${arr[0]}/tadbit_results ${tadbit_resolution} ${genome_gem} ${genome_fa} ${site} ${threads} ${is_shallow}
   																											    
  #JUICER
  echo -e "--- JUICER ${arr[0]}\n"
  mkdir ${arr[2]}/${arr[0]}/juicer_results #creating dir that will contain juicer outputs for sample i 
  mkdir ${arr[2]}/${arr[0]}/straw_results
  cd ${arr[2]}/${arr[0]}/juicer_results #setting juicer_results dir of sample i as working directory. Important because default "topDir" in juicer.sh is the cwd. So topDir will be juicer_results (and it will change at each interaction, for each sample)

  # -d juicer arguments is "topDir". We don't specity it since, by default it is cwd. -n = sample name
  bash ${scriptsDir}/juicer.sh -D ${scriptsDir} -g ${genomeID} -p ${genomePath} -z ${genome_fa} -n ${arr[0]} -s ${site} -y ${site_file} -u ${fastq1} -v ${fastq2} -t ${threads} 

done < <(tail -n +2 ${assoc_file}) #providing association file and skipping header


  
### HICEXPLORER & HICREP - executing outside of the previous while cycle since this script has been implemented to be executed also as stand alone. It loops on association file by itself. ###
echo -e "--- STRAW, HICEXPLORER & HICREP\n"

# executing hicexplorer, differential TADs analysis ed hicrep for each provided resolution #
res_arr=() 
IFS=',' read -r -a res_arr <<< "${hicexplorer_resolution}" 
echo "Resolutions choosed are: ${res_arr[@]}" 
 
arr_len="${#res_arr[@]}"  #lunghezza dell array 1 (che è uguale a quella dell'array 2 di solito) 
for (( i=0; i<=${arr_len}-1; i++ ))    #faccio -1 perchè arr_len è il numero di elementi dell'array, ma gli indici iniziano da 0, quindi devo sottrarre 1 al num totale di elementi 
do
echo -e "--- Executing hicexplorer with resolution ${res_arr[i]}"
bash ${scriptsDir}/hicexplorer_hicrep.sh ${assoc_file} ${res_arr[i]} ${scriptsDir} ${threads}

# STRAW - executed for each the same array of resolution provided for hicexplorer # - DA SCORPORARE
#echo -e "--- STRAW\n"
#mkdir -p ${arr[2]}/${arr[0]}/straw_results/${res_arr[i]}
#python3 ${scriptsDir}/straw.py ${arr[2]}/${arr[0]}/juicer_results/aligned/inter_30.hic ${res_arr[i]} ${arr[2]}/${arr[0]}/straw_results

done



