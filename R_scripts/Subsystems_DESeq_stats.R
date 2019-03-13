# Subsystems_DESeq_stats.R
# Created 3/06/2017, by Sam Westreich
# Last updated 6/16/2017
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
              help="raw (total) read counts for this starting file", metavar="character"),
  make_option(c("-L", "--level"), type="integer", default=1,
              help="level of Subsystems hierarchy for DESeq stats [default=%default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ Subsystems_DESeq_stats.R -I input_directory/ -O save.filename -L level (1,2,3,4) [-R raw_counts_file]")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

cat ("Saving results as ", opt$out, "\n")
save_filename <- opt$out

if (is.null(opt$raw_counts)) {
  print ("WARNING: no raw counts file specified, skipping this info for DESeq analysis.")
} else {
  counts_file <- opt$raw_counts
}

cat ("Calculating DESeq results for hierarchy level ", opt$level, "\n")
levelname=paste("Level", opt$level, sep="")

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
  library("data.table")
  library(plyr)
})

# get list of files
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
  pattern = "experiment_*", full.names = T, recursive = FALSE)
exp_names = ""
for (name in exp_files) {
  exp_names <- c(exp_names, unlist(strsplit(name, split='_', fixed=TRUE))[2])}
exp_names <- exp_names[-1]
exp_names_trimmed = ""
for (name in exp_names) {
  exp_names_trimmed <- c(exp_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
exp_names_trimmed <- exp_names_trimmed[-1]

# sanity check
if (length(exp_files) == 0 || length(control_files) == 0) {
  print ("\nWARNING: No files found.  Is the directory correct?  Are the files named with 'control_' and 'experimental_' as prefixes?")
  stop()
}

# loading the control files in
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(control_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      control_table = control_table[,c(2, 5)]
    } else if (opt$level == 2) {
      control_table = control_table[,c(2, 6)]
    } else if (opt$level == 3) {
      control_table = control_table[,c(2, 4)]
    } else {
      control_table = control_table[,c(2, 3)]
    }
    control_table <- ddply(control_table, colnames(control_table)[2], numcolwise(sum))
    rownames(control_table) <- control_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    control_table <- merge(temp_table, control_table, by = colnames(temp_table)[1], all=T)
  }
}
control_table <- control_table[!is.na(names(control_table))]
control_table[is.na(control_table)] <- ""

# Need to convert NAs to 0s
control_data <- control_table[,c(2:(length(control_files)+1))]
control_data <- lapply(control_data, function(x) as.numeric(as.character(x)))
control_data <- as.data.frame(control_data)
control_data[is.na(control_data)] <- 0
control_table[,c(2:(length(control_files)+1))] <- control_data

# loading the experimental files in
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(exp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      exp_table = exp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      exp_table = exp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      exp_table = exp_table[,c(2, 4)]
    } else {
      exp_table = exp_table[,c(2, 3)]
    }
    exp_table <- ddply(exp_table, colnames(exp_table)[2], numcolwise(sum))
    rownames(exp_table) <- exp_table[,1] }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    if (opt$level == 1) {
      temp_table = temp_table[,c(2, 5)]
    } else if (opt$level == 2) {
      temp_table = temp_table[,c(2, 6)]
    } else if (opt$level == 3) {
      temp_table = temp_table[,c(2, 4)]
    } else {
      temp_table = temp_table[,c(2, 3)]
    }
    temp_table <- ddply(temp_table, colnames(temp_table)[2], numcolwise(sum))
    rownames(temp_table) <- temp_table[,1]
    exp_table <- merge(temp_table, exp_table, by = colnames(temp_table)[1], all=T)
  }
}
exp_table <- exp_table[!is.na(names(exp_table))]
exp_table[is.na(exp_table)] <- ""

# converting NAs to 0s
exp_data <- exp_table[,c(2:(length(exp_files)+1))]
exp_data <- lapply(exp_data, function(x) as.numeric(as.character(x)))
exp_data <- as.data.frame(exp_data)
exp_data[is.na(exp_data)] <- 0
exp_table[,c(2:(length(exp_files)+1))] <- exp_data

# reducing stuff down to avoid duplicates
colnames(control_table) <- c(colnames(control_table)[1], control_names_trimmed)
colnames(exp_table) <- c(colnames(exp_table)[1], exp_names_trimmed)
l1_table <- merge(control_table, exp_table, by=colnames(exp_table)[1], all.x = T)
levelname=colnames(l1_table)[1]
l1_table[,levelname][is.na(l1_table[,levelname])] <- ""
l1_table[,levelname] <- sub("^$", "NO HIERARCHY", l1_table[,levelname])
l1_table <- ddply(l1_table, colnames(l1_table)[1], numcolwise(sum))
rownames(l1_table) <- l1_table[,levelname]
l1_names <- l1_table[,levelname]
l1_table[,levelname] <- NULL
l1_table[is.na(l1_table)] <- 0

# OPTIONAL: importing the raw counts
cat ("Now importing raw counts, if provided.\n")
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
  raw_counts_table["ANNOTATION COUNT",] <- colSums(l1_table)
  raw_counts_table["OTHER",] <- raw_counts_table[1,] - raw_counts_table[2,]
  
  l1_table <- rbind(l1_table, raw_counts_table["OTHER",])
  l1_names <- c(l1_names, "OTHER")
  rownames(l1_table) <- l1_names
}

# now the DESeq stuff
cat ("Now running DESeq.\n")
completeCondition <- data.frame(condition=factor(c(rep("control", length(control_files)), 
  rep("experimental", length(exp_files)))))
dds <- DESeqDataSetFromMatrix(l1_table, completeCondition, ~ condition)
dds <- DESeq(dds)
baseMeanPerLvl <- sapply( levels(dds$condition), function(lvl) rowMeans( 
  counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )

res <- results(dds, contrast = c("condition", "experimental", "control"))
l1_results <- data.frame(res)
rownames(l1_results) <- rownames(l1_table)
rownames(baseMeanPerLvl) <- rownames(l1_table)
l1_results <- merge(l1_results, baseMeanPerLvl, by="row.names")
l1_results <- l1_results[,c(1,2,8,9,3,4,5,6,7)]
colnames(l1_results)[c(3,4)] <- c("controlMean", "experimentalMean")
l1_results <- l1_results[order(-l1_results$baseMean),]

# saving and finishing up
cat ("\nSuccess!\nSaving results file as ", save_filename, "\n")
write.table(l1_results, file = save_filename, 
  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
