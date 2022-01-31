#!/bin/bash

#le colonne del file tabulato sono fisse
#non è importante in quale cartella si runna lo script


FILE=$1
RESOLUTION=$2
MAXTHREADS=$3

#salvo il path in cui eseguo lo script perchè in esso è presente anche il file .tsv nei cicli while successivi sarà necessario darglielo
#altrimenti, dato che nel primo while cambiamo cd, non riuscirà più a leggerlo.
execution_dir=$(pwd)
echo "$execution_dir"

#Arrays for Normalization of muliple samples
types_t=()     #array che conterrà i tipi di trattamento dei sample treated (T1, T2, ecc)
num_rows_no_none=0  #inizializzo un contatore che conterà quante sono le righe del file .tsv riguardanti campioni trattamento. e quindi quanti sono i campioni trattamento, qualsiasi sia il trt.


#Arrays for HICREP steps
all_samples_names=()   #conterrà i nomi di tutti i samples, sia untreated che treated (ogni tipo di treatment)
project_path=()        #conterrà il path della cartella progetto, quella che contiene tutti samples

## 1° CICLO DI LETTURA DEL FILE .TSV ##

### Leggo tutte le righe del file .tsv (una per iterazione). Se il tipo di Treatment non c'è nell array "types_t", lo aggiungo ad esso
### Se il tipo di trattamento è indicato come "None", allora essi sono sample "Untreated", quindi salvo  i loro nomi nell'array "untreated"
### Nello stesso ciclo, converto le matrici dal formato .hic a quello .cool, .h5 ed .mcool, che serviranno successivamente

while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT T_Type  plotMatrix_coor plotTADs_coor CorrectMatrix_Chr  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le righe
do
  arr=($sample $ID $path_dir $matrix_name $T_UT $T_Type $plotMatrix_coor $plotTADs_coor $CorrectMatrix_Chr) #salvo le variabili in un array per comodità
  cd ${arr[2]}
  mkdir hicrep_analysis_${RESOLUTION}    #spostare questa linea fuori dal while, sennò prova a crearla ogni volta e da l errore -  #l'analisi hicrep NON varia al variare della normalizzazione in quanto necessita di matrici .mcool e non normalized.h5
  printf "\n Processing sample ${arr[0]} \n"
  cd ${arr[2]}/${arr[0]}                                        #3° colonna: path project directory + sample name (indici partono da 0)
  echo "la working directory è: ${PWD}"
  mkdir ${RESOLUTION}_resolution

  printf "\n Sample ${arr[0]} is TREATED with treatment ${arr[5]} \n"
  if [[ ! " ${types_t[*]} " =~ " ${arr[5]} " ]]; then
     types_t+=(${arr[5]})                                    #salvo il nome del tipo ti trattamento nell'apposito array
  fi
     num_rows_no_none=$((num_rows_no_none + 1))              #incrementa il contatore ogni volta che legge una riga di un campione treated


  printf "\n >>>>>>>>> ${arr[0]} --> hicConverFormat \n"
  hicConvertFormat -m ${arr[3]}.hic  --inputFormat hic --outputFormat cool -o ${RESOLUTION}_resolution/${arr[3]}.cool --resolution ${RESOLUTION}
  hicConvertFormat -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool --inputFormat cool --outputFormat h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5
  hicConvertFormat -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool --inputFormat cool --outputFormat mcool -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.mcool --resolutions ${RESOLUTION}

done < <(tail -n +2 ${FILE})                                     #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)

echo "il num di righe non None (e quindi il numero di sample trattati) è: ${num_rows_no_none}"
echo "i tipi di trattamento sono:  ${types_t[@]}"
echo "i campioni UNTREATET sono: ${untreated[@]}"

cd ${execution_dir}
echo "la current directory è:"
pwd

## 2° CICLO DI LETTURA DEL FILE .TSV ##

####Per ogni TIPO DI TRATTAMENTO precedentemente identificato nel .tsv, per ogni riga del .tsv raggruppo nella variabile "sample_names" i nomi dei sample aventi quel tipo di trattamento
####E utilizzo i sample names per generare i comandi di hicNormalize specifici per normalizzare assieme le matrici dei sample treated  aventi lo stesso tipo di trattamento
for t in ${types_t[@]};
do
   echo "considero il trattamento $t per indicare il ciclo"
   sample_names=()                                                           #conterrà i sample names di sample trattati con uno stesso tipo di trattamento ad ogni iterazione "for t in ${types_t[@]}"
   norm_input_trt=()						             #conterrà i path delle matrici .h5 degli specifici treated samples da normalizzare
   norm_output_trt=()			                       	             #conterrà i path delle matrici output normalizzate " normalized.h5" degli specifici treated samples normalizzati
   while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT T_Type plotMatrix_coor plotTADs_coor CorrectMatrix_Chr
   do
      arr=($sample $ID $path_dir $matrix_name $T_UT $T_Type $plotMatrix_coor $plotTADs_coor $CorrectMatrix_Chr)
        cd ${arr[2]}
  	if [[ ${t} = ${arr[5]} ]]; then
    	sample_names+=(${arr[0]})
        norm_input_trt+=(${arr[0]}/${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.h5)
        norm_output_trt+=(${arr[0]}/${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_normalized.h5)
        fi

      done < <(tail -n +2 ${FILE})

      echo "i sample di trattamento ${t} sono: ${sample_names[@]}"
      echo "il comando per normalizzare è: ${norm_input_trt[@]}"
      echo "il comando finale è:  hicNormalize -m ${norm_input_trt[@]} -n smallest -o ${norm_output_trt[@]}"

      printf "\n"
      echo ">>>>>>>>> hicNormalize - normalizzo assieme i campioni treated con treatment ${t}: ${sample_names[@]} "
      hicNormalize -m $(echo "${norm_input_trt[@]}") -n smallest -o $(echo "${norm_output_trt[@]}")                                

      cd ${execution_dir}               #torno nella directory da dove ho runnato lo script ed in cui è presente il file .tsv per rileggerlo nuovamente, cercando i samples trattati con il trattamento t successivo (es, al primo ciclo ha cercato i T1, ora i T2 ecc fino a Tn)
      echo "la current directory è:"
      pwd

done



## 3° CICLO DI LETTURA DEL FILE .TSV ##

###DA DOPO NORMALIZE - CorrectMatrix, FindTADs, DetectLoops

while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT T_Type  plotMatrix_coor plotTADs_coor CorrectMatrix_Chr
do
  arr=($sample $ID $path_dir $matrix_name $T_UT $T_Type  $plotMatrix_coor $plotTADs_coor $CorrectMatrix_Chr)
  cd ${arr[2]}/${arr[0]}

  printf "\n >>>>>>>>> ${arr[0]} --> hicCorrectMatrix \n"
  hicCorrectMatrix diagnostic_plot --matrix ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_normalized_diagnostic.png
  hicCorrectMatrix correct -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_normalized.h5 -o ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_corrected.h5

  printf "\n >>>>>>>>> ${arr[0]} --> hicFindTADs \n"
  hicFindTADs -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}_corrected.h5 --outPrefix ${RESOLUTION}_resolution/tads_hic_corrected --numberOfProcessors ${MAXTHREADS} --correctForMultipleTesting fdr --maxDepth $((${RESOLUTION}*10)) --thresholdComparison 0.05 --delta 0.01
  printf "\n >>>>>>>>> ${arr[0]} --> hicDetectLoops \n"
  hicDetectLoops -m ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/loops.bedgraph --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05 --threads ${MAXTHREADS}

  #printf "\n >>>>>>>>>> ${arr[0]} --> make_tracks_file \n"
  #make_tracks_file --trackFiles ${RESOLUTION}_resolution/${arr[3]}_${RESOLUTION}.corrected.h5 ${RESOLUTION}_resolution/tads_hic_corrected_boundaries.bed -o tracks.ini 

  #save all sample's name in "all_samples_names" array for hicrep subsequent step
  all_samples_names+=(${arr[0]})

  #save project folder path in "project_path" array for hicrep subsequent step
  if [ "${#project_path[@]}" -lt 1 ]; then
     project_path+=(${arr[2]})
  fi

done < <(tail -n +2 ${FILE})



### HICREP ###

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

    if [[ "${sample1}" != "${sample2}" ]]; then   #evita confronti tra lo stesso campione
    printf "\n ------> Comparison between samples: ${sample1}-${sample2} \n"
    hicrep ${project_path[0]}/${sample1}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/${sample2}/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.mcool ${project_path[0]}/hicrep_analysis_${RESOLUTION}/hicrep_${sample1}-${sample2}_SCC1.txt  --binSize ${RESOLUTION}  --h ${h} --dBPMax 500000
    fi

   done
        samples=(${all_samples_names[@]:1:100})  #non so specificare diversamente che deve fare un subset dall'elemento 1 all ultimo elemento dell array
done
