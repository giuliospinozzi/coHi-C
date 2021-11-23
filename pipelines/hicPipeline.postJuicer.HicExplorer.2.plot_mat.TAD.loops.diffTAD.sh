#!/bin/bash

#il file deve avere 2 colonne di regioni per plotMatrix e 2 colonne di regioni per plotTADs. analizziamo tutte e 4 le regioni per ogni campione

FILE=$1
RESOLUTION=$2
MAXTHREADS=$3

treated=()
untreated=()
project_path=()

while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT plotMatrix_coor1 plotMatrix_coor2 plotTADs_coor1 plotTADs_coor2 CorrectMatrix_Chr  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le ri>
do
  arr=(${sample} ${ID} ${path_dir} ${matrix_name} ${T_UT} ${plotMatrix_coor1} ${plotMatrix_coor2} ${plotTADs_coor1} ${plotTADs_coor2} ${CorrectMatrix_Chr}) #salvo le variabili in un array per comodità
  echo "Processing sample ${arr[0]}"
  cd ${arr[2]}
  mkdir analysis  #creo una cartella in cui andranno i file rejected, accepted prodotti da hicDifferentialTADs
  cd ${arr[2]}/${arr[0]}
  echo "la working directory è: ${arr[2]}/${arr[0]}"

   ##ciclo tra le 2 diverse regioni da plottare per lo stesso plot (colonna 5 e 6) e salvo la coordinata per intero in "coor", poi in due while separo "chr" da "start" ed "end" per poterli inserire nel nome del file output. sfrutto i cicli per scrivere una sola volta il comando hicPlotMatrix
   printf "\n>>>>>>>>>> ${arr[0]} --> hicPlotMatrix \n"
   if [ ${arr[5]} !=  "None" ]; then   #not equal to
     for coor_matrix in ${arr[@]:5:2};
       do
          while IFS=":" read -r chr start_end; do
            while IFS="-" read -r start end; do

          echo "region-specific command for ${arr[0]}"  #region-specific command
          hicPlotMatrix -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 --clearMaskedBins --region $coor_matrix -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_log2_${chr}_${start}-${end}_matrix_plot.png --log1p --dpi 300

            done <<< ${start_end}
           done <<< ${coor_matrix}
          done
   else
	echo "genome-wide analysis"
   fi

   printf "\n>>>>>>>>>> ${arr[0]} --> hicPlotTADs and plotMatrix with Loops\n" #nello script di giulio, le stesse regioni plottate con plotTADs vengono plottate in plotMatrix con parametro --loops
   if [ ${arr[7]} !=  "None" ]; then
       echo "region-specific command for ${arr[0]}" #region-specific commands
       for coor_TADs in ${arr[@]:7:2};
       do
         echo "le coordinate intere sono: ${coor_TADs}"
           while IFS=":" read -r chr start_end; do
             while IFS="-" read -r start end; do

           hicPlotTADs --tracks tracks.ini -o ${RESOLUTION}_resolution/${arr[3]}_TADs_${chr}_${start}-${end}_track.png --region ${coor_TADs} --dpi 300
           hicPlotMatrix -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 -o ${RESOLUTION}_resolution/${arr[0]}_${RESOLUTION}_${chr}_${start}-${end}_matrix_loop.png --log1p --region ${coor_TADs} --loops ${RESOLUTION}_resolution/loops.bedgraph --dpi 300

             done <<< ${start_end}
           done <<< ${coor_TADs}
        done
    else
        echo "genome-wide analysis"
    fi

#salviamo i path dei campioni T nell'array "treated" e quelli dei campioni UT nell array "untreated" per l analisi differenziale successiva

   if [ ${arr[4]} = "T" ]; then
      echo "aggiungo il nome del campione all array treated"
      treated+=(${arr[0]})
   else
      echo "aggiungo il nome del campione all array untreated"
      untreated+=(${arr[0]})
   fi

#salvo il path del progetto nella variabile "project_path" inizializzata all'esterno del while all inizio
   if [ "${#project_path[@]}" -lt 1 ]; then
    project_path+=(${arr[2]})
   fi



done < <(tail -n +2 ${FILE})  #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)


echo " i campioni treated sono: ${treated[@]}"
echo " i campioni untreated sono: ${untreated[@]}"

# ANALISI DIFFERENZIALE #
#capire se va bene farla dopo aver chiuso il while che legge il file .tsv e che considera un sample alla volta (che quindi può essere T o UT) e dire "se è T, allora fallo confrontare con tutti gli UT". Poi passa all'altro sample (altra riga). Se è UT skippo, se è un altro T lo faccio confrontare con tutti gli UT. Tuttavia mi serve fare un altro ciclo che mi controlli se sono UT o T, quindi boh, forse meglio fare il tutto a parte.
#però se faccio a parte ho il problema del nome del sample e del nome della matrice (inter_30)

  printf "\n >>>>>>>>>> hicDifferentialTAD \n"
  for t_sample in $treated; do
     for ut_sample in $untreated; do

      hicDifferentialTAD -cm ${project_path[0]}/${ut_sample}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -tm ${project_path[0]}/${t_sample}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td ${project_path[0]}/${t_sample}/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o ${project_path[0]}/analysis/differential_tads_${ut_sample}-${t_sample} -p 0.01 -t 4 -mr all --threads ${MAXTHREADS}

      done
  done
