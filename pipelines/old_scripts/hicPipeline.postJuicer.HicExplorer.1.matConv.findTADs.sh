#!/bin/bash

#le colonne del file tabulato sono fisse
#non è importante in quale cartella lo si fa partire

FILE=$1
RESOLUTION=$2
MAXTHREADS=$3

#Arrays for hicrep step
samples=()
project_path=()

while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT plotMatrix_coor plotTADs_coor CorrectMatrix_Chr  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le righe
do
  arr=($sample $ID $path_dir $matrix_name $T_UT $plotMatrix_coor $plotTADs_coor $CorrectMatrix_Chr) #salvo le variabili in un array per comodità
  cd ${arr[2]}
  mkdir hicrep_analysis_${RESOLUTION}
  echo "Processing sample ${arr[0]}"
  cd ${arr[2]}/${arr[0]} #3° colonna: path project directory + sample name (indici partono da 0)
  echo "la working directory è: ${PWD}"
  mkdir ${RESOLUTION}_resolution


  printf "\n >>>>>>>>> ${arr[0]} --> hicConverFormat \n"
  hicConvertFormat -m ${arr[3]}.hic  --inputFormat hic --outputFormat cool -o ${RESOLUTION}_resolution/${arr[3]}.cool --resolution ${RESOLUTION}
  hicConvertFormat -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool --inputFormat cool --outputFormat h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5
  hicConvertFormat -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool --inputFormat cool --outputFormat mcool -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.mcool --resolutions ${RESOLUTION}

  printf "\n >>>>>>>>> ${arr[0]} --> hicNormalize \n"
  hicNormalize -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5 -n smallest -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicCorrectMatrix \n"
  hicCorrectMatrix diagnostic_plot --matrix ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.diagnostic.png
  hicCorrectMatrix correct -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicFindTADs \n"
  hicFindTADs -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 --outPrefix ${RESOLUTION}_resolution/tads_hic_corrected --numberOfProcessors ${MAXTHREADS} --correctForMultipleTesting fdr --maxDepth $((${RESOLUTION}*10)) --thresholdComparison 0.05 --delta 0.01
  printf "\n >>>>>>>>> ${arr[0]} --> hicDetectLoops \n"
  hicDetectLoops -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 --threads ${MAXTHREADS}

  printf "\n >>>>>>>>>> ${arr[0]} --> make_tracks_file \n"
# make_tracks_file --trackFiles ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 ${RESOLUTION}_resolution/tads_hic_corrected_boundaries.bed -o tracks.ini 


#save sample's name in "samples" array for hicrep subsequent step
samples+=(${arr[0]})

#save project folder path in "project_path" array for hicrep subsequent step
   if [ "${#project_path[@]}" -lt 1 ]; then
    project_path+=(${arr[2]})
   fi

done < <(tail -n +2 ${FILE}) #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)


echo "Samples to compare are: ${samples[@]}"

for sample1 in ${samples[@]}; do
    for sample2 in ${samples[@]}; do

    if [[ "${sample1}" != "${sample2}" ]]; then   #evita confronti tra lo stesso campione
    printf "\n ------> Comparison between samples: ${sample1}-${sample2} \n"
    hicrep ${project_path[0]}/${sample1}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/${sample2}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/hicrep_analysis_${RESOLUTION}/hicrep_${sample1}-${sample2}_SCC1.txt  --binSize 10000  --h 20 --dBPMax 500000
    fi

   done
        samples=(${samples[@]:1:100})  #non so specificare diversamente che deve fare un subset dall'elemento 1 all ultimo elemento dell array
done
