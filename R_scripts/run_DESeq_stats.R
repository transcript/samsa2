# Created 6/22/16, updated 6/16/2017
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="./",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="DESeq_results.tab", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-R", "--raw_counts"), type="character", default=NULL,
              help="raw (total) read counts for this starting file", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ run_DESeq_stats.R -I working_directory/ -O save.filename")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

if (is.null(opt$out)) {
  print ("WARNING: No save name for DESeq results specified; defaulting to 'DESeq_results.tab'.") 
  save_filename <- opt$out
} else { save_filename <- opt$out }

if (is.null(opt$raw_counts)) {
  print ("WARNING: no raw counts file specified, skipping this info for DESeq analysis.")
} else {
  counts_file <- opt$raw_counts
}

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
})

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

# merging the two tables
complete_table <- merge(control_table_trimmed, exp_table_trimmed, by=0, all = TRUE)
complete_table[is.na(complete_table)] <- 0
rownames(complete_table) <- complete_table$Row.names
complete_table <- complete_table[!(complete_table$Row.names == ""), ]
# reducing stuff down to avoid duplicates
complete_table <- complete_table[, lapply(.SD, sum), by=complete_table$Row.names]
complete_table <- aggregate(. ~  Row.names, data = complete_table, sum)
# removing extra Row.names column
complete_table <- complete_table[,-1]

# OPTIONAL: importing the raw counts
if (is.null(opt$raw_counts) == FALSE) {
  raw_counts_table <- read.table(counts_file, header=FALSE, sep = "\t", quote = "")
  raw_counts_table <- data.frame(raw_counts_table, 
        do.call(rbind, strsplit(as.character(raw_counts_table$V1),'_')))
  raw_counts_table$X2 <- as.numeric(as.character(raw_counts_table$X2))
  raw_counts_table <- t(raw_counts_table[,c("X2", "V2")])
  row.names(raw_counts_table) <- c("SAMPLE","RAW TOTAL")
  colnames(raw_counts_table) <- raw_counts_table[1,]
  raw_counts_table <- as.data.frame(raw_counts_table)
  raw_counts_table <- raw_counts_table[-1,]
  
  # Need to subtract off the total number of annotations
  raw_counts_table["ANNOTATION COUNT",] <- colSums(complete_table)
  raw_counts_table["OTHER",] <- raw_counts_table[1,] - raw_counts_table[2,]

  complete_table <- rbind(complete_table, raw_counts_table["OTHER",])
}

# DESeq statistical calculations
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
cat ("\nSuccess!\nSaving results file as ", save_filename, "\n")
write.table(sorted_org_results, file = save_filename, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
