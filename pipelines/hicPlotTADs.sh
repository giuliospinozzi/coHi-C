#!/bin/bash

### HIC PLOT TADS ### - for multiple resolution, for multiple regions, for multiple samples.

#Script in cui plotto soltanto i TADs. 
#Input necessari: tracks.ini (e quindi la matrice .corrected.h5 e tads_hic_corrected_domains.bed) - lo trova in automatico

#EXECUTION EXAMPLE:
#./hicPlotTADs.sh /home/alessio/hic/complete_hic_pipe_exe/momentaneo_ADA_da_eliminare/association_file_shallow_ADA.tsv 5000,10000 8 chr11:33875000-34135000,chr11:33864303-34012900


FILE=$1	#association file
RESOLUTION=$2 #one or more resolution, comma separated (e. g. 5000,10000)
MAXTHREADS=$3  #credo che qui non serva e debba quindi essere rimosso, da confermare
COORDINATES=$4  #coordinate delle regioni da plottare (comma separated)


arr_coor=()
IFS=',' read -r -a arr_coor <<< "${COORDINATES}" 
echo "l'array contenente le coordinate è:  ${arr_coor[@]}" 


printf "\n>>>>>>>>>> ${arr[0]} --> hicPlotTADs \n" #plotto ogni regione, per ogni campione, una volta per ogni risoluzione desiderata. Ad es, risoluzione 5000, 10000 ecc.

res_arr=() 
IFS=',' read -r -a res_arr <<< "${RESOLUTION}" 
echo "Resolutions choosed are: ${res_arr[@]}" 

arr_len="${#res_arr[@]}"  #lunghezza dell array 1 (che è uguale a quella dell'array 2 di solito) 

#ciclo per ogni risoluzione
for (( i=0; i<=${arr_len}-1; i++ ))    #faccio -1 perchè arr_len è il numero di elementi dell'array, ma gli indici iniziano da 0, quindi devo sottrarre 1 al num totale di elementi 
do
echo -e "--- Executing hicPlotTADs with resolution ${res_arr[i]}"
	#ciclo per ogni coordinata
	for coor in ${arr_coor[@]}
	do
		echo "la coordinata da plottare è:" ${coor}
		
		#ciclo per ogni sample
    	while IFS=$'\t' read -r sample ID path_dir T_UT T_Type Fastq1_Dir Fastq2_Dir
	  	do
	  	arr=(${sample} ${ID} ${path_dir} ${T_UT} ${T_Type} ${Fastq1_Dir}	${Fastq2_Dir} )
	  	printf "\n #### PROCESSING SAMPLE \"${arr[0]}\" #### \n"
	  	cd ${arr[2]}
	  	
	  	#Creating "fastq" dir in "topDir" if its not already present 
		if [ ! -d "plotTADs_${res_arr[i]}_${coor}" ]; then
    		mkdir plotTADs_${res_arr[i]}_${coor}  
		else 
    		echo -e "---  Using already created plotTADs_${res_arr[i]}_${coor}\n"		
		fi
	  	
	  		while IFS=":" read -r chr start_end; do
				while IFS="-" read -r start end; do

		       		hicPlotTADs --tracks ${arr[0]}/hicexplorer_results/${res_arr[i]}_resolution/tracks.ini -o ${arr[2]}/plotTADs_${res_arr[i]}_${coor}/${arr[0]}_${res_arr[i]}_TADs_${chr}_${start}-${end}.png --region ${coor} --dpi 300


				done <<< ${start_end}
		  	done <<< ${coor}


        done < <(tail -n +2 ${FILE})



	done  #chiudo ciclo per coordinate

done	#chiudo ciclo per risoluzioni





