# diversity_graphs.R
# Created 8/11/16
# Last updated 6/16/2017
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-I", "--input"), type="character", default="./",
              help="Input directory", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="diversity_graph.pdf", 
              help="output file name [default= %default]", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ diversity_graphs.R -I working_directory/ -O save.filename")

# check for necessary specs
if (is.null(opt$input)) {
  print ("WARNING: No working input directory specified with '-I' flag.")
  stop()
} else {  cat ("Working directory is ", opt$input, "\n")
  wd_location <- opt$input  
  setwd(wd_location)  }

cat ("Saving diversity graphs as", opt$out, "\n")
save_filename <- opt$out

# import other necessary packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(scales)
  library(reshape2)
  library(knitr)
  library(vegan)
  library(gridExtra)
  library(ggplot2)
})

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

# merging the two tables together
complete_table <- merge(control_table_trimmed, exp_table_trimmed, by=0, all = TRUE)
complete_table[is.na(complete_table)] <- 1
rownames(complete_table) <- complete_table$Row.names
complete_table <- complete_table[,-1]

# getting diversity statistics
flipped_complete_table <- data.frame(t(complete_table))

graphing_table <- data.frame(condition=factor(c(rep("control", length(control_files)), 
  rep("experimental", length(exp_files)))))
graphing_table[,"order"] <- c(1:nrow(graphing_table))
graphing_table[,"Shannon"] <- diversity(flipped_complete_table, index = "shannon")
graphing_table[,"Simpson"] <- diversity(flipped_complete_table, index = "simpson")

shannon_plot <- ggplot(data = graphing_table, aes(x=order, y=Shannon, 
  color = condition, fill = condition)) + 
  geom_bar(stat="identity", width = 0.8) +
  ggtitle("Shannon diversity of control vs experimental samples") +
  theme(legend.position = "bottom")

simpson_plot <- ggplot(data = graphing_table, aes(x=order, y=Simpson, 
  color = condition, fill = condition)) + 
  geom_bar(stat="identity", width = 0.8) +
  ggtitle("Simpson diversity of control vs experimental samples") +
  theme(legend.position = "bottom")

cat ("\nSuccess!\nSaving diversity graphs as ", save_filename, " now.\n")
pdf(file = save_filename, width=10, height=7)
grid.arrange(shannon_plot, simpson_plot, ncol=1)
dev.off()
