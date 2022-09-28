#!/usr/bin/env Rscript

#This script create barplots of number of total TADs and loops detected in each sample, both with number of differential TADs in each sample comparison performed through TADs differential analysis



#### DATAFRAME of number of TADs, number of Loops and TADs mean length of each sample ####

args = commandArgs(trailingOnly=TRUE)

#saving arguments in variables
path_tsv <- args[1] #path of association file
resolution <- args[2] #resolution used in hicexplorer
project_path <- args[3] #path of project directory that contains all the samples output


association_file <- read.table(file = path_tsv, sep = '\t', header = TRUE) #open association file
samples <- association_file[, 1] #saving samples names

TADs_Loops_df <- data.frame(matrix(ncol=3, nrow=0)) #initializing dataframe that will contains TADs and loops stats
col_names <- c("N_TADs", "Mean_Length_TADs", "N_Loops") 

#cycle to get stats of each sample
for (sample in samples) {
  row_sample <- c()  #initializing an empty row of TADs_Loops_df
  
  sample_tads_domains_bed <- read.table(paste(project_path, "/", sample, "/hicexplorer_results/", resolution, "_resolution/tads_hic_corrected_domains.bed", sep = "")) #opening sample tads_hic_corrected_domains.bed file
  cat(sample, "sample stats:", "\n")
  
  n_TADs <- length(sample_tads_domains_bed$V2) #extracting number of TADs of the sample
  cat("Number of TADs:", n_TADs, "\n")
  
  #mean length of sample TADs computation
  length_TADs <- sample_tads_domains_bed$V3-sample_tads_domains_bed$V2
  mean_length_TADs <- mean(length_TADs)
  cat("TADs mean length:", mean_length_TADs, "\n")
  
  #Extracting number of Loops of the sample
  sample_loops_bedgraph <- read.table(paste(project_path, "/", sample, "/hicexplorer_results/", resolution, "_resolution/loops.bedgraph", sep = ""))
  n_LOOPS <- length(sample_loops_bedgraph$V1)
  cat("Number of Loops:", n_LOOPS, "\n", "\n")
  
  #Appending sample stats to TADs_Loops_df
  row_sample <- append(row_sample, c(n_TADs, mean_length_TADs, n_LOOPS))
  TADs_Loops_df <- rbind(TADs_Loops_df, row_sample)
  
}

colnames(TADs_Loops_df) <- col_names  
rownames(TADs_Loops_df) <- samples

name <- paste(project_path, "/stats_plots/", resolution, "_resolution/", "TADs_Loops_df.csv", sep = "") #creating dataframe file name

cat("Saving TADs and Loops dataframe as .csv file", "\n")
write.csv(TADs_Loops_df, name, row.names = TRUE)



#### BARPLOTS of number of TADs, number of Loops and TADs mean length of each sample #####

library(ggplot2)
library(ggpubr) #for "ggarange" function

cols <- hcl.colors(length(samples), "Zissou 1") #using this palette for barplot

#barplot of number of TADs
plot_N_TADs <- ggplot(TADs_Loops_df, aes(x = samples, y = N_TADs, fill = samples)) +
  geom_bar(stat = "identity" ) +
  theme(axis.title = element_text(size = 22 )) + #axis labels 
  theme(axis.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18 )) +
  theme(legend.text = element_text(size = 16 )) + 
  scale_fill_manual(values=cols) # or + scale_fill_brewer(palette = "Dark2")

#barplot of TADs mean length
plot_mean_length_TADs <- ggplot(TADs_Loops_df, aes(x = samples, y = Mean_Length_TADs, fill = samples)) +
  geom_bar(stat = "identity" ) +
  theme(axis.title = element_text(size = 22 )) + #axis labels
  theme(axis.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18 )) +
  theme(legend.text = element_text(size = 16 )) +
  scale_fill_manual(values=cols) # or + scale_fill_brewer(palette = "Dark2")

#barplot of number of loops
plot_N_Loops <- ggplot(TADs_Loops_df, aes(x = samples, y = N_Loops, fill = samples)) +
  geom_bar(stat = "identity" ) +
  theme(axis.title = element_text(size = 22 )) + #axis labels
  theme(axis.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18 )) +
  theme(legend.text = element_text(size = 16 )) +
  scale_fill_manual(values=cols) # or + scale_fill_brewer(palette = "Dark2")

# PNG device

#getting the name of the project (es "hicADA" to write it as title on the plot) using strsplit base R function 
project_name <- unlist(strsplit(project_path, "/"))
project_name <- project_name[length(project_name)]

png(paste(project_path, "/stats_plots/", resolution, "_resolution/TADs_loops_plots.png",  sep = ""), width = 1400, height = 1000, units = "px")

plot_title <- paste(project_name, "-", "Res:", resolution,  sep = " ")


#arranging the three separated plots in one single figure 
grouped_plots <- ggarrange(plot_N_TADs, plot_mean_length_TADs, plot_N_Loops, 
                           labels = c("A", "B", "C"),
                           ncol = 2, nrow = 2)

#annotating figure with title
annotate_figure(grouped_plots, top = text_grob(plot_title, 
                                               color = "firebrick2", face = "bold", size = 30))

cat("Saving TADs and Loops barplot as .png file", "\n")

# Close device
dev.off()



#### DIFFERENTIAL TADS DATAFRAME AND BARPLOT ####
library(stringr)

diff_TADs_path <- paste(project_path, "/diff_TADs_analysis/", resolution, "_resolution", sep = "") #path containing differential TADs files
diff_TADs_df <- data.frame(matrix(ncol=6, nrow=0)) #initializing differential TADs dataframe

rejected_files <- list.files(path=diff_TADs_path, pattern="*rejected.diff_tad", full.names=TRUE, recursive=FALSE)  

comparisons <- c() #initializing samples comparison vector
N_diff_TADs <- c() #initializing number of differential TADs vector

#cycle to get number of differential TADs for each samples comparison
for(i in 1:length(rejected_files)) {
  
  rej_diff_TAD_table <- read.table(rejected_files[i])
  
  #extracting sample comparison from a differential tads file (es: "SAMPLE1-SAMPLE2")
  conf <- rejected_files[i]
  conf <- unlist(strsplit(conf, "/"))
  conf <- conf[length(conf)]
  conf <- str_remove(conf, "_rejected.diff_tad")
  conf <- str_remove(conf, "differential_tads_")
  
  comparisons <- append(comparisons, conf) #append sample comparison to comparison vector
  N_diff_TADs <- append(N_diff_TADs, length(rej_diff_TAD_table$V1)) #appending number of differential TADs for this samples comparison (it is equal to the number of rows in that file. each row is a differential TAD)
  
}

cat("Comparisons:", comparisons, "\n")
cat("Number of differential TADs per comparison:", N_diff_TADs, "\n", "\n")

#creating dataframe with number of differential TADs
diff_TADs_df <- data.frame(comparisons, N_diff_TADs)
rownames(diff_TADs_df) <- diff_TADs_df[, 1]
diff_TADs_df[, 1] <- NULL
print(diff_TADs_df)
cat("\n")

name_csv <- paste(project_path, "/stats_plots/", resolution, "_resolution/", "diff_TADs_df.csv", sep = "") #name of dataframe containing number of differential TADS

cat("Saving number of differential TADs dataframe as .csv file", "\n")
write.csv(diff_TADs_df, name_csv, row.names = TRUE)



#### BARPLOT of number of differential TADs for each samples comparison ####

library(RColorBrewer)
display.brewer.all()

# PNG device
png(paste(project_path, "/stats_plots/", resolution, "_resolution/diffTADs_plot.png",  sep = ""), width = 1400, height = 1000, units = "px")


cols <- hcl.colors(length(comparisons), "Zissou 1")

diff_TADs_plot <- ggplot(diff_TADs_df, aes(x = comparisons, y = N_diff_TADs, fill = comparisons, )) +
  geom_bar(stat = "identity" ) +
  theme(axis.title = element_text(size = 22 )) + #i labes degli assi 
  theme(axis.text = element_text(size = 16)) +
  theme(legend.title = element_text(size = 18 )) +
  theme(legend.text = element_text(size = 16 )) + 
  scale_fill_manual(values=cols) # or + scale_fill_brewer(palette = "Dark2")

plot_title <- paste(project_name, "-", "Res:", resolution,  sep = " ")

#Adding title to figure with plot
annotate_figure(diff_TADs_plot, top = text_grob(plot_title, 
                                               color = "firebrick2", face = "bold", size = 30))

cat("Saving number of differential TADs barplot as .png file", "\n")

# Close device
dev.off()



















