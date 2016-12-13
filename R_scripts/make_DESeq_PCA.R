# make_DESeq_PCA.R
# Created 6/28/16
# Arguments that need to be specified: working directory (ARGV1), save filename (ARGV2)

args <- commandArgs(TRUE)

suppressPackageStartupMessages({
  library(DESeq2) })
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
  control_names <- c(control_names, unlist(strsplit(name, split='_', fixed=TRUE))[4])}
control_names <- control_names[-1]
control_names_trimmed = ""
for (name in control_names) {
  control_names_trimmed <- c(control_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
control_names_trimmed <- control_names_trimmed[-1]

exp_files <- list.files(
  pattern = "experimental_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[4])}
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
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(control_table) = c("DELETE", x, "V3")
    control_table <- control_table[,c(2,3)]      }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    temp_table <- temp_table[,c(2,3)]
    control_table <- merge(control_table, temp_table, by = "V3", all = T)  }
}
control_table[is.na(control_table)] <- 0
rownames(control_table) = control_table$V3
control_table_trimmed <- control_table[,-1]

# loading the experimental table
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header=F, quote = "", sep = "\t")
    colnames(exp_table) = c("DELETE", x, "V3")
    exp_table <- exp_table[,c(2,3)]  }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t")
    colnames(temp_table) = c("DELETE", x, "V3")
    exp_table <- merge(exp_table, temp_table[,c(2,3)], by = "V3", all = T)  }
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
complete_array <- data.matrix(data_table_filtered)

# making the PCA plot

# calculate euclidean distances from the variance-stabilized data
dists <- dist(t(assay(transformed_data)))
PCAplot <- plotPCA(transformed_data, intgroup = "Dose", returnData = TRUE)
percentVar <- round(100 * attr(PCAplot, "percentVar"))

# saving and finishing up
cat ("Saving PCA plot as ", save_filename, " now.\n")
pdf(file = paste(save_filename, ".pdf", sep = ""), width=10, height=7)
ggplot(func_PCAplot, aes(PC1, PC2, color=Individual, shape=Dose)) +
    geom_point(size=3) +
#    geom_text(aes(label=name), hjust=1, vjust=-1) +
    ggtitle("PCA Plot of control vs. experimental FUNCTIONAL data") +
    theme(legend.position = "bottom") +
#    xlim(-35, 25) + ylim(-25, 40) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()