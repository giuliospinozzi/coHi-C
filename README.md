# coHi-C - a Comprehensive Pipeline for HiC Data Processing developed at SR-Tiget

## Introduction

**coHi-C** is a comprehensive pipeline for Hi-C data analysis that covers all stages of analysis integrating both existing tools and novel algorithms to ensure versatility, accuracy, and efficiency. It can be executed through **Bash** command line, is fully automated and multithreaded. To ensure its proper functioning, a series of **arguments** must be provided, including the absolute path of the related **association file**, a tab-separated text file that must be manually prefilled with all the information related to samples and case studies. **coHi-C** is made by eight scripts: each one executes a set of functions provided by the various adopted tools with the only exception of Juicer.sh, which is a modified version of the original one to fit in the flow, to perform specific tasks not foreseen by its developers and with some shrewdness from the point of view of memory space optimization.

***coHiC.sh*** is the principal script, the one that the user needs to run. It is written in bash and what it does is automatically execute all the other scripts sequentially, except for *hicPlotTADs.sh* and *hicPlotLoops.sh* which are intended to be executed only after the end of the complete run, in a separate way.

First, **Assembly-stats** provided information about the input FASTQs like the number of reads and their average length.

Then TADbit assesses the quality of the performed HiC experiment through the generation of plots and statistics as the percentage of dangling-end and undigested reads and the trend of genomic interactions based on the distance between two loci.

Juicer is next in charge to build up the *.hic* contact matrices at multiple resolutions (bin size) from raw sequencing data (pairs of genomic positions that were adjacent to each other in 3D space during the experiment). It includes a collection of algorithms designed to annotate contact matrices and hence recognize features of genome folding including loops, loop anchor motifs, and contact domain. 

Post-processing of the interaction matrices is taken into account by HiCExplorer’s tools: it converts and normalizes *.hic* matrices into its own native formats and into *.mcool* format for the subsequent reproducibility assessment through HiCRep. Then, after detecting TADs and loops, differential TADs analysis is computed across different sample matrices whose results are summarized through some descriptive statistics plots.

The entire workflow is executed once for each sample.

## Prerequisites

...

### Applications

**For Graphical User Interface and scripts running:**

- Python 3.8 (modules: pytadbit, os, sys) ([https://www.python.org/downloads/](https://www.python.org/downloads/))
- R version 4.1.3 ([https://www.r-project.org/](https://www.r-project.org/)) (Packages: ggplot2, ggpubr, stringr, RColorBrewer)

**For quality control on HiC experiment, on fastq files and pre-processing:**

- TADbit 1.0.1 ([https://github.com/3DGenomes/TADbit](https://github.com/3DGenomes/TADbit))
- assembly-stats 1.0.1 ([https://github.com/sanger-pathogens/assembly-stats](https://github.com/sanger-pathogens/assembly-stats))

**For alignment, BAM files manipulation, building .hic interaction matrix:** 

- Juicer 1.6 ([https://github.com/aidenlab/juicer](https://github.com/aidenlab/juicer))
- samtools 1.15.1 ([http://www.htslib.org/](http://www.htslib.org/))

**For TADs & Loops plots, and differential TADs analysis:**

- HiCExplorer 3.7.2 ([https://github.com/deeptools/HiCExplorer](https://github.com/deeptools/HiCExplorer))

**For HiC experiment reproducibility assessment:**

- hicrep 0.2.6 ([https://github.com/TaoYang-dev/hicrep](https://github.com/TaoYang-dev/hicrep))


### Input files
The necessary input files are as follows:

**Association file (AF):** A tab-separated **file** in which each row contains information about a single sample.  All samples must belong to the same experimental study. Columns contain:

- **Sample_name:** name of the sequenced sample.
- **ID:** ID of the study
- **Outdir_Abs_Path_Dir:** Absolute path of the pipeline output’s directory.
- **T/UT:** Label indicating if a sample comes from a treated (T) or an untreated (UT) sample. It is necessary to perform differential TADs analysis.
- **T_Type:** Treated samples could be treated in different ways in the same experiment. For this reason, this variable contains the type of treatment to which a sample was subjected. Became relevant during the matrices normalization step.
- **Fastq_R1_Abs_Path:** Absolute path to sample forward (R1) FASTQ file.
- **Fastq_R2_Abs_Path:** Absolute path to sample reverse (R2) FASTQ file.

Other input files that must be provided as arguments at the execution are:

- Reference genome necessary for alignment with GEM2. Its absolute path must be provided as -m argument at the execution and its relative indexes must be in its same directory.
- Reference genome, necessary for alignment with BWA. Its absolute path must be provided as -z argument at the execution and its relative indexes must be in its same directory.
- Locations of restriction sites in reference genome necessary for alignment with BWA file. Its absolute path must be provided as -y argument at the execution.
- File containing all chromosome sizes of the reference genome necessary for alignment with BWA file. Its absolute path must be provided as -p argument at the execution.

**Creating indexes:**

The indexes can be created starting from the Fasta file and require different commands based on the alignment software for which they will be used.

In the case of Gem2, the necessary command is: (not sure)

`gem-indexer -i genome.fasta -o genome_index`

In the case of BWA, the necessary command is:

`bwa index file.fa`

## Running Application

### Command line

After installing the necessary tools and getting the required input files, you can launch the application using command-line. The following is the use for analysis with the arguments that must necessarily be inserted, as they do not provide default options.

### Usage

`coHIC.sh -a <association_file.tsv> -D <abs_path_script_dir> -g <reference_genome_ID> -z <abs_path_reference_genome_fasta_file> -m <abs_path_reference_genome_gem_file> -s <restriction_enzyme_name> -y <abs_path_restriction_site_file> -p <abs_path_reference_genome_chrom.sizes_file> -r <tadbit_resolution> -R <hicexplorer_resolutions> -t <max_threads> -l <is_shallow>` 

### Arguments

| Command | Description |
| --- | --- |
| `-a` | Association file absolute path. No default options. |
| `-D` | Absolute path to the directory containing all the scripts (both with Juicer "common" directory). No default options. |
| `-g` | Reference genome ID. Default option: "hg19" |
| `-z` | Reference genome fasta file absolute path. No default options. |
| `-m` | Reference genome .gem file absolute path. No default options. |
| `-s` | Restriction enzyme name. No default options |
| `-y` | Absolute path to restriction site file, containing genomic coordinates the of restriction sites. |
| `-p` | Path to the chrom.sizes file of the corresponding reference genome.  |
| `-r` | Resolution used by TADbit |
| `-R` | HiCExplorer resolutions. Can be one or more commas divided (es: 5000,10000). No default options. |
| `-t` | Max thread number. No default options. |
| `-l` | TADbit execution mode: complete execution for shallow sequences (-l true), only quality plots for full sequences (-l false). No default options |
| `-h` | Print usage page |
