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
  3. RESOLUTION [10000]

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
RESOLUTION="${3}";

#printf "Folder creation --> ${RESULTS_DIR}/${PROJECT_NAME}/${POOL_NAME}\n"
#mkdir ${RESULTS_DIR}/${PROJECT_NAME}

##### ==================== End SETTINGS ======================== #####




# Activate Anaconda Environment (HiC)
conda activate hic


for k in $(ls -d $WORKING_DIR/*/); do
    printf "\n###############################################"
    printf "\n ----- Processins Sample --> ${k#"$WORKING_DIR/"}"
    printf "\n############################################### \n"
    SAMPLE=${k#"$WORKING_DIR/"}
    SAMPLE=${SAMPLE::len-1}
    cd $WORKING_DIR/$SAMPLE
    mkdir $WORKING_DIR/$SAMPLE/${RESOLUTION}_resolution

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicConvertFormat \n"
    hicConvertFormat -m inter_30.hic --inputFormat hic --outputFormat cool -o ${RESOLUTION}_resolution/inter_30.cool --resolutions ${RESOLUTION}
    hicConvertFormat -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --inputFormat cool --outputFormat h5 -o ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.h5

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicCorrectMatrix \n"
    hicCorrectMatrix correct -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.h5 --filterThreshold -1.5 5 -o ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicMergeMatrixBins \n"
    hicMergeMatrixBins -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -o ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.nb50.h5 -nb 50

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPCA \n" ##### Too slow!!!!
    hicAdjustMatrix --matrix ${RESOLUTION}_resolution/inter_30_5000.corrected.h5 --outFileName ${RESOLUTION}_resolution/inter_30_5000_chr2-11.corrected.h5 --chromosomes 2 11 --action keep
    hicPCA -m ${RESOLUTION}_resolution/inter_30_5000_chr2-11.corrected.h5 --outputFileName ${RESOLUTION}_resolution/pca1.bw ${RESOLUTION}_resolution/pca2.bw --format bigwig

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicTransform \n"
    hicTransform -m ${RESOLUTION}_resolution/inter_30_5000_chr2-11.corrected.h5 --outFileName ${RESOLUTION}_resolution/pearson_chr2-11.h5 --method pearson --perChromosome

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotMatrix: A/B compartments can be plotted \n"
    hicPlotMatrix -m ${RESOLUTION}_resolution/pearson_chr2-11.h5 --outFileName ${RESOLUTION}_resolution/pca1.png --perChromosome --bigwig ${RESOLUTION}_resolution/pca1.bw --dpi 300

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotMatrix \n"
    # Edit for specific project
    hicPlotMatrix -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 --clearMaskedBins --region chr2:60221755-61279177 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_log2_chr2_60221755-61279177_matrix_plot.png --log1p --dpi 300
    hicPlotMatrix -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 --clearMaskedBins --region chr11:4769502-5825416 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_log2_chr11_4769502-5825416_matrix_plot.png --log1p --dpi 300

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicFindTADs \n"
    hicFindTADs -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 --outPrefix ${RESOLUTION}_resolution/tads_hic_corrected --numberOfProcessors $MAXTHREADS --chromosomes 2 11 --correctForMultipleTesting fdr --maxDepth $((${RESOLUTION}*10)) --thresholdComparisons 0.05 --delta 0.01
    
    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotTADs \n"
    # Edit for specific project
    hicPlotTADs --tracks ${RESOLUTION}_resolution/track.ini -o ${RESOLUTION}_resolution/${SAMPLE}_TADs_chr11_5019502-5575416_track.png --region chr11:5019502-5575416 --dpi 300
    hicPlotTADs --tracks ${RESOLUTION}_resolution/track.ini -o ${RESOLUTION}_resolution/${SAMPLE}_TADs_chr2_60471755-61029177_track.png --region chr2:60471755-61029177 --dpi 300

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicDetectLoops \n"
    hicDetectLoops -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/loops.bedgraph --chromosomes 2 11 --maxLoopDistance 2000000 --windowSize 10 --peakWidth 6 --pValuePreselection 0.05 --pValue 0.05
    hicPlotMatrix -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr11_5019502-5575416_matrix_loop.png --log1p --region 11:5019502-5575416 --loops ${RESOLUTION}_resolution/loops.bedgraph --dpi 300
    hicPlotMatrix -m ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr2_60471755-61029177_matrix_loop.png --log1p --region 2:60471755-61029177 --loops ${RESOLUTION}_resolution/loops.bedgraph --dpi 300

    printf "\n>>>>>>>>>> ${SAMPLE} --> hicPlotViewpoint \n"
    hicPlotViewpoint --matrix ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --region 2:60471755-61029177 --referencePoint 2:60721755-60722677 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr2_60471755-61029177_matrix_ViewPoint.png --dpi 300
    hicPlotViewpoint --matrix ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --region 2:60471755-61029177 --referencePoint 2:60778398-60779177 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr2_60778398-60779177_matrix_ViewPoint.png --dpi 300
    hicPlotViewpoint --matrix ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --region 11:5019502-5575416 --referencePoint 11:5269502-5271087 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr11_5269502-5271087_matrix_ViewPoint.png --dpi 300
    hicPlotViewpoint --matrix ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --region 11:5019502-5575416 --referencePoint 11:5274418-5276011 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr11_5274418-5276011_matrix_ViewPoint.png --dpi 300
    hicPlotViewpoint --matrix ${RESOLUTION}_resolution/inter_30_${RESOLUTION}.cool --region 11:5019502-5575416 --referencePoint 11:5291155-5325416 -o ${RESOLUTION}_resolution/${SAMPLE}_${RESOLUTION}_chr11_5291155-5325416_matrix_ViewPoint.png --dpi 300
done



# Edit for specific project (All Steps)
printf "Starting Comparative Analysis --> $WORKING_DIR/analysis \n"
mkdir $WORKING_DIR/analysis
mkdir $WORKING_DIR/analysis/${RESOLUTION}_resolution/
cd $WORKING_DIR/

#All Samples
#hicCompareMatrices -m E22/inter_30_${RESOLUTION}.corrected.nb50.h5 E23/inter_30_${RESOLUTION}.corrected.nb50.h5 --operation log2ratio -o analysis/E22-E23_log2_m50.h5

#hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr11_5019502-5575416_track.png --region chr11:5019502-5575416
#hicPlotTADs --tracks track.ini -o ${SAMPLE}_TADs_chr2_60471755-61029177_track.png --region chr2:60471755-61029177

#hicAdjustMatrix -m inter_30_${RESOLUTION}.corrected.h5 --action keep --chromosomes 2 -o matrix_chr1.h5

hicDifferentialTAD -tm UT1/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm E22/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT1/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT1-E22 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT1/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm E23/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT1/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT1-E23 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT1/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm P10/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT1/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT1-P10 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT1/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm P11/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT1/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT1-P11 -p 0.01 -t 4 -mr all

hicDifferentialTAD -tm UT2/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm E22/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT2/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT2-E22 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT2/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm E23/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT2/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT2-E23 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT2/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm P10/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT2/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT2-P10 -p 0.01 -t 4 -mr all
hicDifferentialTAD -tm UT2/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -cm P11/${RESOLUTION}_resolution/inter_30_${RESOLUTION}.corrected.h5 -td UT2/${RESOLUTION}_resolution/tads_hic_corrected_domains.bed -o analysis/${RESOLUTION}_resolution/differential_tads_UT2-P11 -p 0.01 -t 4 -mr all

