#!/bin/bash
source /etc/environment
source /etc/profile
source /opt/applications/scripts/miniconda3/etc/profile.d/conda.sh


echo "
  +--------------------------------------------------------+
  |                                                        |
  |                 HiC Pipeline (SR-Tiget)                |
  |                                                        |
  +--------------------------------------------------------+
  |  Author:   Giulio Spinozzi                             |
  |  Date:     July 2021                                   |
  |  Contact:  spinozzi.giulio@hsr.it                      |
  |  Version:  0.1 - from Juicer Data                      |
  +--------------------------------------------------------+

  REQUIRED VARS and relative ORDER POSITION -> REMEMBER NO SPACES!!!!!!!!!
	1. WORKING_DIR [/opt/ngs/results] NO / at the end !!
	2. MAXTHREADS [12]

  NOTES:
  	- Input are post Juicer run
  	- Based on HiCexplorer 
"



echo "<`date +'%Y-%m-%d %H:%M:%S'`> [SR-Tiget] Preprocessing input variables (delimiters:<>)"
## print input variables (check for log utils)
INPUTVARNUM=0
for INPUTVAR in "$@"; do
	let INPUTVARNUM++; 
	printf -v INPUTNUM '%02d' $INPUTVARNUM;
    echo "  => Input Variable: Order number = <${INPUTNUM}> ; Var Content = <${INPUTVAR}>
    ";
done



##### =================== Start SETTINGS ======================= #####
WORKING_DIR="${1}";
MAXTHREADS="${2}";

#printf "Folder creation --> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}\n"
#mkdir ${RESULTS_DIR}/${PROJECT_NAME}

##### ==================== End SETTINGS ======================== #####




# Activate Anaconda Environment (HiC)
conda activate hic

for k in $(ls -d $WORKING_DIR/*); do
    printf "\n###############################################"
    printf "\n ----- Processins Sample --> ${k#"$WORKING_DIR/"}"
    printf "\n############################################### \n"
    SAMPLE=${k#"$WORKING_DIR/"}
    cd $WORKING_DIR/$SAMPLE

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicConvertFormat \n"
    hicConvertFormat -m inter_30.hic --inputFormat hic --outputFormat cool -o inter_30.cool --resolutions 10000
    hicConvertFormat -m inter_30_10000.cool --inputFormat cool --outputFormat h5 -o inter_30_10000.h5

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicCorrectMatrix \n"
    hicCorrectMatrix correct -m inter_30_10000.h5 --filterThreshold -1.5 5 -o inter_30_10000.corrected.h5

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicMergeMatrixBins \n"
    hicMergeMatrixBins -m inter_30_10000.corrected.h5 -o inter_30_10000.corrected.nb50.h5 -nb 50

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPCA \n" ##### Too slow!!!!
    #hicPCA --matrix inter_30_10000.corrected.h5 -o pca1.bedgraph pca2.bedgraph -f bedgraph

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotMatrix \n"
    # Edit for specific project
    hicPlotMatrix -m inter_30_10000.corrected.h5 --clearMaskedBins --region chr2:60221755-61279177 -o ${SAMPLE}_10000_log2_chr2_60221755-61279177_matrix_plot.png --log1p
    hicPlotMatrix -m inter_30_10000.corrected.h5 --clearMaskedBins --region chr11:4769502-5825416 -o ${SAMPLE}_10000_log2_chr11_4769502-5825416_matrix_plot.png --log1p

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicFindTADs \n"
    hicFindTADs -m inter_30_10000.corrected.h5 --outPrefix tads_hic_corrected --numberOfProcessors $MAXTHREADS --correctForMultipleTesting fdr --maxDepth 100000 --thresholdComparisons 0.05 --delta 0.01
    
    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotTADs \n"
    # Edit for specific project
    hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr11_5019502-5575416_track.png --region chr11:5019502-5575416
    hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr2_60471755-61029177_track.png --region chr2:60471755-61029177
done



# Edit for specific project (All Steps)
printf "Starting Comparative Analysis --> $WORKING_DIR/analysis \n"
mkdir $WORKING_DIR/analysis
cd $WORKING_DIR/

#All Samples
#hicCompareMatrices -m E22/inter_30_10000.corrected.nb50.h5 E23/inter_30_10000.corrected.nb50.h5 --operation log2ratio -o analysis/E22-E23_log2_m50.h5

#hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr11_5019502-5575416_track.png --region chr11:5019502-5575416
#hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr2_60471755-61029177_track.png --region chr2:60471755-61029177

#hicDifferentialTAD -tm UT1/inter_30_10000.corrected.h5 -cm E22/inter_30_10000.corrected.h5 -td UT1/UT1_domains.bed -o anayslis/differential_tads -p 0.01 -t 4 -mr all
