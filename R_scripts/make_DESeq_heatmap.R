# make_DESeq_heatmap.R
# Created 6/27/16
# Arguments that need to be specified: working directory (ARGV1), save filename (ARGV2)

args <- commandArgs(TRUE)

library(DESeq2, quietly=TRUE)
library("pheatmap", quietly = TRUE)
library("RColorBrewer", quietly = TRUE)
library("ggplot2", quietly = TRUE)

# check for necessary specs
if (!is.na(args[1])) {
  cat ("Working directory is ", args[1], "\n")
  wd_location <- args[1]
} else {
  print ("WARNING: No working directory specified as ARGV1") 
  stop() }

if (!is.na(args[2])) {
  save_filename <- args[2]
} else {
  print ("WARNING: No name of saved PCA plot specified as ARGV2")
  stop() }

setwd(wd_location)

# GET FILE NAMES
control_files <- list.files(
  pattern = "control_*", full.names = T, recursive = FALSE)
control_names = ""
for (name in control_files) {
  control_names <- c(control_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
control_names <- control_names[-1]
control_names_trimmed = ""
for (name in control_names) {
  control_names_trimmed <- c(control_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
control_names_trimmed <- control_names_trimmed[-1]

exp_files <- list.files(
  pattern = "experimental_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
exp_names <- exp_names[-1]
exp_names_trimmed = ""
for (name in exp_names) {
  exp_names_trimmed <- c(exp_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
exp_names_trimmed <- exp_names_trimmed[-1]

# READ IN FILES
# loading the control table
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header=F, quote = "", sep = "\t", fill = TRUE)
    if (ncol(control_table) == 4) {
      colnames(control_table) = c("DELETE", x, "V3", "md5")
      control_table <- control_table[,c(4,2,3)] 
    } else {
      colnames(control_table) = c("DELETE", x, "V3")
      control_table <- control_table[,c(3,2)] }
    if (nrow(control_table) > 1000) {
      control_table <- control_table[c(1:1000),] }}     # can be deleted, restricts to top 1k hits
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    if (nrow(temp_table) > 1000) {
      temp_table <- temp_table[c(1:1000),] }      # can be deleted, restricts to top 1k hits
    print (x)
    if (ncol(temp_table) == 4) {
      colnames(temp_table) = c("DELETE", x, "V3", "md5") 
    } else {
      colnames(temp_table) = c("DELETE", x, "V3") }
    control_table <- merge(control_table, temp_table[,c(2,3)], by = "V3", all.x = T)  }
}
control_table[is.na(control_table)] <- 0
rownames(control_table) = control_table$V3
control_table_trimmed <- control_table[,-1, drop = FALSE]

# loading the exp table
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header=F, quote = "", sep = "\t", fill = TRUE)
    if (ncol(exp_table) == 4) {
      colnames(exp_table) = c("DELETE", x, "V3", "md5")
      exp_table <- exp_table[,c(4,2,3)] 
    } else {
      colnames(exp_table) = c("DELETE", x, "V3")
      exp_table <- exp_table[,c(3,2)] }
    if (nrow(exp_table) > 1000) { exp_table <- exp_table[c(1:1000),] }} # can be deleted, restricts to top 1k hits
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    if (nrow(temp_table) > 1000) { temp_table <- temp_table[c(1:1000),] }    # can be deleted, restricts to top 1k hits
    print (x)
    if (ncol(temp_table) == 4) {
      colnames(temp_table) = c("DELETE", x, "V3", "md5") 
    } else {
      colnames(temp_table) = c("DELETE", x, "V3") }
    exp_table <- merge(exp_table, temp_table[,c(2,3)], by = "V3", all.x = T)  }
}
exp_table[is.na(exp_table)] <- 0
rownames(exp_table) = exp_table$V3
exp_table_trimmed <- exp_table[,-1]

# getting the column names simplified
colnames(control_table_trimmed) = control_names_trimmed
colnames(exp_table_trimmed) = exp_names_trimmed

complete_table <- merge(control_table_trimmed, exp_table_trimmed, by=0, all = TRUE)
complete_table[is.na(complete_table)] <- 1
rownames(complete_table) <- complete_table$Row.names
complete_table <- complete_table[,-1]
completeCondition <- data.frame(condition=factor(c(
  rep(paste("control", 1:length(control_files), sep=".")), 
  rep(paste("experimental", 1:length(exp_files), sep=".")))))
completeCondition1 <- t(completeCondition)
colnames(complete_table) <- completeCondition1
completeCondition2 <- data.frame(condition=factor(c(
  rep("control", length(control_files)), 
  rep("experimental", length(exp_files)))))

dds <- DESeqDataSetFromMatrix(complete_table, completeCondition2, ~condition)

dds <- DESeq(dds)
transformed_data <- rlog(dds, blind=FALSE)
complete_array <- data.matrix(complete_table)

# making the PCA plot

# calculate euclidean distances from the variance-stabilized data
dists <- dist(t(assay(transformed_data)))

# Create heatmap of distances
sampleDistMatrix <- as.matrix( dists )
rownames(sampleDistMatrix) <- colnames(complete_array)
colnames(sampleDistMatrix) <- transformed_data$condition
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=dists,
         clustering_distance_cols=dists,
         show_rownames = TRUE,
         show_colnames = TRUE,
         col=colors)

# saving and finishing up
cat ("Saving PCA plot as ", save_filename, " now.\n")
pdf(file = paste("macaque_function_heatmap_3-28", ".pdf", sep = ""), width=10, height=7)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=dists,
         clustering_distance_cols=dists,
         show_rownames = TRUE,
         show_colnames = TRUE,
         col=colors)
dev.off()
