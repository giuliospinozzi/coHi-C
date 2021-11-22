#!/bin/bash

#le colonne del file tabulato sono fisse 

FILE=$1
RESOLUTION=$2
MAXTHREADS=$3

while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT plotMatrix_coor plotTADs_coor CorrectMatrix_Chr  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le righe
do
  arr=($sample $ID $path_dir $matrix_name $T_UT $plotMatrix_coor $plotTADs_coor $CorrectMatrix_Chr) #salvo le variabili in un array per comodità
  echo "Processing sample ${arr[0]}"

  echo ${arr[2]} #3° colonna: path (indici partono da 0)
  cd ${arr[2]}
  echo "la working directory è: ${arr[2]}"
  mkdir ${RESOLUTION}_resolution


  printf "\n >>>>>>>>> ${arr[0]} --> hicConverFormat \n"
  hicConvertFormat -m ${arr[3]}.hic  --inputFormat hic --outputFormat cool -o ${RESOLUTION}_resolution/${arr[3]}.cool --resolution ${RESOLUTION}
  hicConvertFormat -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool --inputFormat cool --outputFormat h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicNormalize \n"
  hicNormalize -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5 -n smallest -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicCorrectMatrix \n"
  hicCorrectMatrix diagnostic_plot --matrix ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.diagnostic.png
  hicCorrectMatrix correct -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicFindTADs \n"
  hicFindTADs -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 --outPrefix ${RESOLUTION}_resolution/tads_hic_corrected --numberOfProcessors ${MAXTHREADS} --correctForMultipleTesting fdr --maxDepth $((${RESOLUTION}*10)) --thresholdComparison 0.05 --delta 0.01
  printf "\n >>>>>>>>> ${arr[0]} --> hicDetectLoops \n"
  hicDetectLoops -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05

  printf "\n >>>>>>>>>> ${arr[0]} --> make_tracks_file \n"
  make_tracks_file --trackFiles ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 ${RESOLUTION}_resolution/tads_hic_corrected_boundaries.bed -o tracks.ini




done < <(tail -n +2 ${FILE}) #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)

