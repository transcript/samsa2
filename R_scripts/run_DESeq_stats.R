# run_DESeq_stats.R
# Created 6/22/16
# Arguments that need to be specified: working directory (ARGV1), save filename (ARGV2)

args <- commandArgs(TRUE)

library(DESeq2, quietly=TRUE)

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
  print ("WARNING: No name of saved results file specified as ARGV2")
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

baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )

# This step creates the summary results output
res <- results(dds)
org_results <- data.frame(res)
org_results <- merge(org_results, baseMeanPerLvl, by="row.names")
org_results <- org_results[,c(1,2,8,9,3,4,5,6,7)]
colnames(org_results)[c(3,4)] <- c("controlMean", "experimentalMean")

sorted_org_results <- org_results[order(-org_results$baseMean),]
colnames(sorted_org_results)[1] <- "Organism Name"

# saving and finishing up
cat ("Saving results file as ", save_filename, "\n")
write.table(sorted_org_results, file = save_filename, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
