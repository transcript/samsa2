# make_combined_graphs.R
# Created 6/28/16
# Arguments that need to be specified: working directory (ARGV1), save filename (ARGV2)

args <- commandArgs(TRUE)

library(optparse)
option_list = list(
  make_option(c("-g", "--gtitle"), type="character", default=NULL, 
              help="title to be displayed on graph", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="combined_graph.pdf", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--directory"), type="character", default=NULL,
              help="working directory location", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# check for necessary specs
if (is.null(opt$directory)) {
  print ("WARNING: No working directory specified with '-d' flag.")
  stop()
} else {
  cat ("Working directory is ", opt$directory, "\n")
  wd_location <- opt$directory
}

if (is.null(opt$out)) {
  print ("WARNING: No save name for plot specified; defaulting to 'combined_graph.pdf'.") }

suppressPackageStartupMessages({
  library(DESeq2)
  library("RColorBrewer")
  library("ggplot2")
  library(gridExtra)
  library(scales)
  library(reshape2)
  library(knitr)
})

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

# melting tables
control_table_trimmed_m <- melt(cbind(control_table_trimmed, 
      Genus = rownames(control_table_trimmed)), id.vars = c('Genus'))

exp_table_trimmed_m <- melt(cbind(exp_table_trimmed, 
      Genus = rownames(exp_table_trimmed)), id.vars = c('Genus'))

# other catchall category
control_table_filtered <- control_table_trimmed
exp_table_filtered <- exp_table_trimmed
if (nrow(control_table_trimmed)>10 && nrow(exp_table_trimmed)>10) {
  control_table_filtered[,"Total"] <- rowSums(control_table_filtered)
  control_table_filtered <- control_table_filtered[ with (control_table_filtered, order(-Total)), ]
  exp_table_filtered[,"Total"] <- rowSums(exp_table_filtered)
  exp_table_filtered <- exp_table_filtered[ with (exp_table_filtered, order(-Total)), ]

  control_table_filtered["Other",] <- colSums(control_table_filtered[10:nrow(control_table_filtered),])
  control_table_filtered <- rbind(control_table_filtered[1:9,], control_table_filtered[nrow(control_table_filtered),])
  exp_table_filtered["Other",] <- colSums(exp_table_filtered[10:nrow(exp_table_filtered),])
  exp_table_filtered <- rbind(exp_table_filtered[1:9,], exp_table_filtered[nrow(exp_table_filtered),])
  control_table_filtered <- control_table_filtered[,-ncol(control_table_filtered)]
  exp_table_filtered <- exp_table_filtered[,-ncol(exp_table_filtered)]
}
  
control_table_filtered_m <- melt(cbind(control_table_filtered,
     Genus = rownames(control_table_filtered)), id.vars = c('Genus'))

exp_table_filtered_m <- melt(cbind(exp_table_filtered, 
     Genus = rownames(exp_table_filtered)), id.vars = c('Genus'))

# merging the tables
temprow <- matrix(c(rep.int(0,length(control_table_trimmed_m))),nrow=1,ncol=length(control_table_trimmed_m))
newrow <- data.frame(temprow)
colnames(newrow) <- colnames(control_table_trimmed_m)
newrow$variable <- "<blank>"
newrow$Genus <- "Filler"

full_table <- rbind(control_table_trimmed_m, newrow, exp_table_trimmed_m)
full_table_top30 <- rbind(control_table_filtered_m, newrow, exp_table_filtered_m)

# color palette
CbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
  "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",  
  "#ffffb3",  "#bebada",  "#fb8072",  "#80b1d3",  "#fdb462",  "#b3de69", "#fccde5",  
  "#d9d9d9",  "#bc80bd",  "#ccebc5",  "#ffed6f", "#e41a1c",  "#377eb8", "#4daf4a",  
  "#984ea3",  "#ff7f00", "#ffff33",  "#a65628",  "#f781bf", "#999999", "#000000", 
  "#a6cee4", "#1f78b5", "#b2df8b", "#a9a9a9", "#ffffff")

# plotting
org_relative_ggplot <- ggplot(full_table_top30, aes(x = variable, y = value, fill = Genus)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = CbPalette) +
  scale_y_continuous(labels = percent_format()) +
  theme(legend.position = "none", text=element_text(size=16), 
        axis.text.x = element_text(angle=90, vjust=1)) +
  guides(fill = guide_legend(ncol=10)) +
#  ggtitle(opt$gtitle) +
  xlab("Sample ID") + ylab("Relative activity of total sample")

org_absolute_ggplot <- ggplot(full_table_top30, aes(x = variable, y = value, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CbPalette) +
  theme(legend.position = "bottom", text=element_text(size=16), 
        axis.text.x = element_text(angle=90, vjust=1)) +
  guides(fill = guide_legend(ncol=4)) +
#  ggtitle("Total organism abundance in thresholded Macaque samples") +
  xlab("Sample ID") + ylab("Total reads per sample")

# both graphs displayed together: 
#grid.arrange(org_relative_ggplot, org_absolute_ggplot, ncol = 1)

# saving and finishing up
cat ("Saving PCA plot as ", opt$out, " now.\n")
pdf(file = "combined_graphs.pdf", width=9, height=12)
grid.arrange(org_relative_ggplot, org_absolute_ggplot, ncol = 1)
dev.off()