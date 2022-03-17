#!/bin/bash

#HIC PIPELINE (script che eseguirà juicer, tadbit, hicexplorer ecc)


#juiceDir (-D): la directory contenente gli SCRIPT di juicer (/opt/applications/scripts/juicer)
#topDir (-d): la directory che conterrà gli output prodotti con juicer (lo script crea "juicer_out" e la rendo cwd. cwd è topDir di default)

#site_file (option -y): il file con i restriction enzyimes sites. Devo averne uno nella apposita dir "restriction_sites" nella juice_

assoc_file=$1  #file di associazione (.tsv)
juiceDir=$2   #path alla directory contenente gli script di juicer (script/common)
genomeID=$3  #es: hg19 (option -g)
genome_path=$4 #path al file .fa del genoma di riferimento (option -z)
site_file=$5  #path al restriction site file 

threads=$6 #num of threads (option -t)



while IFS=$'\t' read -r sample ID path_dir matrix_name T_UT T_Type fastq1_dir fastq2_dir  #salvo ogni valore di ogni colonna in una variabile mentre leggo tutte le righe
do
  arr=($sample $ID $path_dir $matrix_name $T_UT $T_Type $fastq1_dir $fastq2_dir) #salvo le variabili in un array per comodità
  echo -e "--- Analyzing sample ${arr[0]}\n"
  mkdir ${arr[2]}/${arr[0]} #creo la dir del sample 
  mkdir ${arr[2]}/${arr[0]}/juicer_out #creo la dir dove salverò gli output di juicer, specifica di un sample
  cd ${arr[2]}/${arr[0]}/juicer_out #cd alla main dir di un sample (ovvero alla directory dove salverò gli output). importante perchè "topDir" in juicer.sh di default è la cwd. Quindi topDir sarà questa directory (e varierà ad ogni iterazione, per ogni sample) 

  # -d è topDir, la directory in cui finiranno gli output. Non la specifico perchè di default è la cwd, specificata prima. -n = sample name
 bash ${juiceDir}/scripts/juicer_CPU_MOD.sh -D ${juiceDir} -p ${genomeID} -z ${genome_path} -y ${site_file} -n ${arr[0]} -u ${fastq1_dir} -v ${fastq2_dir} -t 8 
 echo ${arr}  
  
done < <(tail -n +2 ${assoc_file}) #fornisco il file da leggere e dico che voglio leggere dalla linea 2 (skippo l'header)


#comando tipo per juicer originale
#./juicer.sh -z /opt/genome/human/hg19/index/hg19.fa -p hg19 -D /home/alessio/hic/juicer_prova


#comando runnando juicer_CPU_MOD_OK senza lo script bash superiore
#./juicer_CPU_MOD.sh -D /opt/applications/scripts/juicer -z /opt/genome/human/hg19/index/hg19.fa -p hg19 -t 8 -n 2_1_S1_L001 -u /home/alessio/hic/juicer_prova/fastq/2_1_S1_L001_R1_001.fastq.gz -v /home/alessio/hic/juicer_prova/fastq/2_1_S1_L001_R2_001.fastq.gz


#comando di esempio con cui runnare questo script.
#./HIC_PIPE_V1.sh /home/alessio/hic/juicer_prova/juicer_AF.tsv /home/alessio/hic/juicer_prova hg19 /opt/genome/human/hg19/index/hg19.fa /home/alessio/hic/juicer_prova/restriction_sites/hg19_DpnII.txt 12




