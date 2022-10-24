#!/bin/bash


FILE=$1
RESOLUTION=$2
R_PLOTS_SCRIPT_DIR=$3
MAXTHREADS=$4


#Arrays for Normalization of muliple samples
types_t=()     #initializing empty array that will contain the type of treatment of treated samples (es: T1, T2 ecc)


#Arrays for HICREP steps
all_samples_names=()   #initializing an empty array that will contain all samples names, both treated and untreated ones
project_path=()        #initializing an empty array that will contain the path of the project directory, the one with all the samples files 


treated=()
untreated=()


### 1° CYCLE OF ASSOCIATION FILE READING ###

### Leggo tutte le righe del file .tsv (una per iterazione). Se il tipo di Treatment non c'è nell array "types_t", lo aggiungo ad esso
### Se il tipo di trattamento è indicato come "None", allora essi sono sample "Untreated", quindi salvo  i loro nomi nell'array "untreated"
### Nello stesso ciclo, converto le matrici dal formato .hic a quello .cool, .h5 ed .mcool, che serviranno successivamente


while IFS=$'\t' read -r sample ID path_dir T_UT T_Type Fastq1_Dir Fastq2_Dir  #saving each values of row i in some corresponding variables
do
  arr=($sample $ID $path_dir $T_UT $T_Type $Fastq1_Dir $Fastq2_Dir) #saving those variables in a single array 
  
  #checking if those directory are already present. If not, create them. 
  if [ ! -e "${arr[2]}/${arr[0]}/hicexplorer_results/${RESOLUTION}_resolution" ]; then
  	mkdir -p ${arr[2]}/${arr[0]}/hicexplorer_results/${RESOLUTION}_resolution
  fi
  
  if [ ! -e "${arr[2]}/hicrep_results/${RESOLUTION}_resolution" ]; then
  	mkdir -p ${arr[2]}/hicrep_results/${RESOLUTION}_resolution
  fi

  #salvo i path piu utilizzati in variabili, cosi da snellire il codice
  juicer_results_dir=${arr[2]}/${arr[0]}/juicer_results/aligned
  hicexp_outdir=${arr[2]}/${arr[0]}/hicexplorer_results/${RESOLUTION}_resolution
  hicrep_outdir=${arr[2]}/hicrep_results/${RESOLUTION}_resolution
  

  printf "\n Sample ${arr[0]} is TREATED with treatment ${arr[4]} \n"  #colonna "T_Type"
  if [[ ! " ${types_t[*]} " =~ " ${arr[4]} " ]]; then
     types_t+=(${arr[4]})                                    #saving treatment type name in the "types_t" array
  fi


  printf "\n >>>>>>>>> ${arr[0]} --> hicConvertFormat \n"
  hicConvertFormat -m ${juicer_results_dir}/inter_30.hic  --inputFormat hic --outputFormat cool -o ${hicexp_outdir}/inter_30.cool --resolution ${RESOLUTION}
  hicConvertFormat -m ${hicexp_outdir}/inter_30_${RESOLUTION}.cool --inputFormat cool --outputFormat h5 -o ${hicexp_outdir}/inter_30_${RESOLUTION}.h5
  hicConvertFormat -m ${hicexp_outdir}/inter_30_${RESOLUTION}.cool --inputFormat cool --outputFormat mcool -o ${hicexp_outdir}/inter_30_${RESOLUTION}.mcool --resolutions ${RESOLUTION}

done < <(tail -n +2 ${FILE})                                     #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)

echo "Samples treatments types are: ${types_t[@]}"


### 2° CYCLE OF ASSOCIATION FILE READING ###
## For each type of treatment, for each row of association file, save the name of samples having that type of treatment in "sample_names" variable ##
## Normalizing together .h5 matrices of samples having the same treatment type, with hicNormalize ##


for t in ${types_t[@]}; #iterating for each treatment type
do
   printf "\n * Considering treatment $t to indicate the cycle * \n"
   sample_names=()                                       #will contain sample names of samples having the t treatment type 
   norm_input_trt=()						             #will contain path of .h5 matrices of sample having t treatment type 
   norm_output_trt=()			                       	 #will contain path of normalized.h5 matrices, output of hicNormalize
   while IFS=$'\t' read -r sample ID path_dir T_UT T_Type Fastq1_Dir Fastq2_Dir
   do
      arr=($sample $ID $path_dir $T_UT $T_Type $Fastq1_Dir $Fastq2_Dir)
      
      hicexp_outdir=${arr[2]}/${arr[0]}/hicexplorer_results/${RESOLUTION}_resolution  #lo devo rimettere perchè sennò hicexp_outdir sarà uguale al path dell'ultimo sample analizzato nel primo ciclo while
      
  	  if [[ ${t} = ${arr[4]} ]]; then
    	  sample_names+=(${arr[0]})
          norm_input_trt+=(${hicexp_outdir}/inter_30_${RESOLUTION}.h5)
          norm_output_trt+=(${hicexp_outdir}/inter_30_${RESOLUTION}_normalized.h5)
      fi

      done < <(tail -n +2 ${FILE})

      echo "Treated samples with treatment ${t} are: ${sample_names[@]}"
      echo "Comand to normalize is: ${norm_input_trt[@]}"
      echo "Final command is:  hicNormalize -m ${norm_input_trt[@]} -n smallest -o ${norm_output_trt[@]}"

      printf "\n"
      echo ">>>>>>>>> hicNormalize - normalizing together samples treated with treatment ${t}: ${sample_names[@]}"
      printf "\n"
      hicNormalize -m $(echo "${norm_input_trt[@]}") -n smallest -o $(echo "${norm_output_trt[@]}")                                

done


### 3° CYCLE OF ASSOCIATION FILE READING ###
## Post NORMALIZE -> CorrectMatrix, FindTADs, DetectLoops ##

while IFS=$'\t' read -r sample ID path_dir T_UT T_Type Fastq1_Dir Fastq2_Dir
do
  arr=($sample $ID $path_dir $T_UT $T_Type $Fastq1_Dir $Fastq2_Dir)

  hicexp_outdir=${arr[2]}/${arr[0]}/hicexplorer_results/${RESOLUTION}_resolution

  printf "\n >>>>>>>>> ${arr[0]} --> hicCorrectMatrix \n"
  hicCorrectMatrix diagnostic_plot --matrix ${hicexp_outdir}/inter_30_${RESOLUTION}_normalized.h5 -o   ${hicexp_outdir}/inter_30_${RESOLUTION}_normalized_diagnostic.png
  hicCorrectMatrix correct -m ${hicexp_outdir}/inter_30_${RESOLUTION}_normalized.h5 -o ${hicexp_outdir}/inter_30_${RESOLUTION}_corrected.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicFindTADs \n"
  hicFindTADs -m ${hicexp_outdir}/inter_30_${RESOLUTION}_corrected.h5 --outPrefix ${hicexp_outdir}/tads_hic_corrected --numberOfProcessors ${MAXTHREADS} --correctForMultipleTesting fdr --maxDepth $((${RESOLUTION}*10)) --thresholdComparison 0.05 --delta 0.01
  printf "\n >>>>>>>>> ${arr[0]} --> hicDetectLoops \n"
  hicDetectLoops -m ${hicexp_outdir}/inter_30_${RESOLUTION}.cool -o ${hicexp_outdir}/loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 --threads ${MAXTHREADS}

  printf "\n >>>>>>>>>> ${arr[0]} --> make_tracks_file \n"
  make_tracks_file --trackFiles ${hicexp_outdir}/inter_30_${RESOLUTION}_corrected.h5 ${hicexp_outdir}/tads_hic_corrected_boundaries.bed -o ${hicexp_outdir}/tracks.ini 

  #save all sample's name in "all_samples_names" array for hicrep subsequent step
  all_samples_names+=(${arr[0]})

  #save project folder path in "project_path" array for hicrep and hicDifferentialTADs subsequent step
  if [ "${#project_path[@]}" -lt 1 ]; then
     project_path+=(${arr[2]})
  fi

# saving path of T samples in "treated" array and path of UT samples in "untreated" array for differential TADs analysis #
   if [[ ${arr[3]} = "T" ]]; then
      echo "Adding sample name "${arr[0]}" to treated array" 
      treated+=(${arr[0]})
   else
      echo "Adding sample name "${arr[0]}" to untreated array"
      untreated+=(${arr[0]})
   fi


done < <(tail -n +2 ${FILE})


echo "Treated samples are: ${treated[@]}"
echo "Untreated samples are: ${untreated[@]}"



### If samples are all treated or all untreated, no differential TADs analysis is performed. Passing "diffTADs" variable as arguments of TADs_loops_plots.R script. If analysis is performed, it creates differential TADs barplot and dataframe. ###
if [[ ${#treated[@]} -eq 0 ]] || [[ ${#untreated[@]} -eq 0 ]]; then
	diffTADs="false"
	printf "\n !! Differential TADs analysis will not be executed: hicDifferentialTAD needs both treated and untreated samples !! \n"
else
	diffTADs="true"
	
	if [ ! -e "${arr[2]}/diff_TADs_analysis/${RESOLUTION}_resolution" ]; then
  		mkdir -p ${arr[2]}/diff_TADs_analysis/${RESOLUTION}_resolution
	fi

### Differential TADs analysis ###
	printf "\n >>>>>>>>>> hicDifferentialTAD \n"
	for t_sample in ${treated[@]}; do
		for ut_sample in ${untreated[@]}; do

    hicDifferentialTAD -cm ${project_path[0]}/${ut_sample}/hicexplorer_results/${RESOLUTION}_resolution/inter_30_${RESOLUTION}_corrected.h5 -tm ${project_path[0]}/${t_sample}/hicexplorer_results/${RESOLUTION}_resolution/inter_30_${RESOLUTION}_corrected.h5 -td ${project_path[0]}/${t_sample}/hicexplorer_results/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o ${project_path[0]}/diff_TADs_analysis/${RESOLUTION}_resolution/differential_tads_${ut_sample}-${t_sample} -p 0.01 -t 4 -mr all --threads ${MAXTHREADS}

		done
	done

fi


### R script that creates dataframes and barplots of number of TADs and loops per sample and of number of differential TADs per samples comparison if differential TADs analysis is performed###
printf "\n >>>>>>>>>> Creating dataframes and barplots of number of TADs and loops per sample and number of differential TADs per sample comparison if differential TADs analysis is performed: \n"

  if [ ! -e "${project_path[0]}/stats_plots/${RESOLUTION}_resolution" ]; then
  	mkdir -p ${project_path[0]}/stats_plots/${RESOLUTION}_resolution
  fi

#need association file, resolution used in hicexplorer, project path and diffTADs parameter
Rscript --vanilla ${R_PLOTS_SCRIPT_DIR}/TADs_loops_plots.R ${FILE} ${RESOLUTION} ${project_path[@]} ${diffTADs}
																																																																																	



### HICREP ###  - usa la matrice .mcool, quindi quest' analisi non differisce in base alla normalizzazione utilizzata, ma differisce in base alla risoluzione
#migliorare questa parte dei parametri in base alla risoluzione: e se ho risoluzioni diverse da 5000 e 10000? E se ne ho più di due? Risolvere


#setup of "h" parameter according to resolution (binSize)
if [[ ${RESOLUTION} = 5000 ]]; then
    echo "With resolution=5000, h=10"
    h=10
elif [[ ${RESOLUTION} = 10000 ]]; then
    echo "With resolution=10000, h=20"
    h=20
fi

echo "Samples to compare are: ${all_samples_names[@]}"

for sample1 in ${all_samples_names[@]}; do
    for sample2 in ${all_samples_names[@]}; do

    if [[ "${sample1}" != "${sample2}" ]]; then   #avoid comparison between the same sample 
    printf "\n ------> Comparison between samples: ${sample1}-${sample2} \n"
    hicrep ${project_path[0]}/${sample1}/hicexplorer_results/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/${sample2}/hicexplorer_results/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/hicrep_results/${RESOLUTION}_resolution/hicrep_${sample1}-${sample2}_SCC1.txt  --binSize ${RESOLUTION}  --h ${h} --dBPMax 500000
    fi

   done
        samples=(${all_samples_names[@]:1:100})  #non so specificare diversamente che deve fare un subset dall'elemento 1 all ultimo elemento dell array
done

																																																					    












