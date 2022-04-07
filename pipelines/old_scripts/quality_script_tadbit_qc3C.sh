#!/bin/bash

r1=$1
r2=$2

arr1=()
arr2=()

IFS=',' read -r -a arr1 <<< "${r1}"
IFS=',' read -r -a arr2 <<< "${r2}"

echo "l'array1 contenente gli R1 è:  ${arr1[@]}"
echo "l'array2 contenente gli R2 è:  ${arr2[@]}"

arr_len="${#arr1[@]}"  #lunghezza dell array 1 (che è uguale a quella dell'array 2 di solito)


for (( i=0; i<=${arr_len}-1; i++ ))    #faccio -1 perchè arr_len è il numero di elemento dell'array, ma gli indici iniziano da 0, quindi devo sottrarre 1 al num totale di elementi
do
   echo "Tadbit for couple $i of fastq"
   tadbit map -w tadbit2 --fastq fastq/${arr1[i]} \
                         --fastq2 fastq/${arr2[i]} \
                         --index /opt/genome/human/hg19/index/gem/hg19.gem \
                         --renz BfuCI Bsp143I BssMI BstMBI DpnII Kzo9I MboI NdeII Sau3AI \
                         --read 1


   echo "qc3C for couple $i of fastq"
   qc3C kmer --mean-insert 300 \
		  --lib kmer_db.jf  \
		  --reads fastq/${arr1[i]} \
		  --reads fastq/${arr2[i]} \
		  --library-kit arima \
		  --output-path qc3C_output_couple_${i}  #devo indicare una cartella nuova per ogni esecuzione di qc3C in cui salverà gli output. La cartella la crea lui.

done
