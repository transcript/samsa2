# make_combined_graphs.R
# Created 6/28/16
# Last updated 6/16/2017
# Run with --help flag for help.

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gtitle"), type="character", default=NULL, 
              help="title to be displayed on graph", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="combined_graph.pdf", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-D", "--directory"), type="character", default=NULL,
              help="working directory location", metavar="character"),
  make_option(c("-N", "--number"), type="integer", default=10,
              help="Number of top organisms to graph; default is 10", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ make_combined_graphs.R -D working_directory/ -O save.filename -N org_number -G title")

# check for necessary specs
if (is.null(opt$directory)) {
  print ("WARNING: No working directory specified with '-d' flag.")
  stop()
} else {
  cat ("Working directory is ", opt$directory, "\n")
  wd_location <- opt$directory
}

if (is.null(opt$out)) {
  print ("WARNING: No save name for plot specified; defaulting to 'combined_graph.pdf'.")
  save_filename <- "combined_graph.pdf"
} else {
  cat ("Saving results as ", opt$out, "\n")
  save_filename <- opt$out
}

suppressPackageStartupMessages({
  library(DESeq2)
  library("ggplot2")
  library(gridExtra)
  library(scales)
  library(reshape2)
  library(knitr)
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

# melting tables
control_table_trimmed_m <- melt(cbind(control_table_trimmed, 
      Genus = rownames(control_table_trimmed)), id.vars = c('Genus'))

exp_table_trimmed_m <- melt(cbind(exp_table_trimmed, 
      Genus = rownames(exp_table_trimmed)), id.vars = c('Genus'))

# other catchall category
control_table_filtered <- control_table_trimmed
exp_table_filtered <- exp_table_trimmed
if (nrow(control_table_trimmed)>opt$number && nrow(exp_table_trimmed)>opt$number) {
  control_table_filtered[,"Total"] <- rowSums(control_table_filtered)
  control_table_filtered <- control_table_filtered[ with (control_table_filtered, order(-Total)), ]
  exp_table_filtered[,"Total"] <- rowSums(exp_table_filtered)
  exp_table_filtered <- exp_table_filtered[ with (exp_table_filtered, order(-Total)), ]

  control_table_filtered["Other",] <- colSums(control_table_filtered[opt$number:nrow(control_table_filtered),])
  control_table_filtered <- rbind(control_table_filtered[1:(opt$number-1),], control_table_filtered[nrow(control_table_filtered),])
  exp_table_filtered["Other",] <- colSums(exp_table_filtered[opt$number:nrow(exp_table_filtered),])
  exp_table_filtered <- rbind(exp_table_filtered[1:(opt$number-1),], exp_table_filtered[nrow(exp_table_filtered),])
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
    "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462", "#b3de69", "#fccde5",  
    "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f", "#e41a1c", "#377eb8", "#4daf4a",  
    "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999", "#000000", 
    "#a6cee4", "#1f78b5", "#b2df8b", "#a9a9a9", "#ffffff")

Cb64k <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
           "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
           "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
           "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
           "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
           "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
           "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
           "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
           "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
           "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
           "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
           "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
           "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

# plotting
org_relative_ggplot <- ggplot(full_table_top30, aes(x = variable, y = value, fill = Genus)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = Cb64k) +
  scale_y_continuous(labels = percent_format()) +
  theme(legend.position = "none", text=element_text(size=16), 
        axis.text.x = element_text(angle=90, vjust=1)) +
  guides(fill = guide_legend(ncol=10)) +
  ggtitle(opt$gtitle) +
  xlab("Sample ID") + ylab("Relative activity of total sample")

org_absolute_ggplot <- ggplot(full_table_top30, aes(x = variable, y = value, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = Cb64k) +
  theme(legend.position = "bottom", text=element_text(size=16), 
        axis.text.x = element_text(angle=90, vjust=1)) +
  guides(fill = guide_legend(ncol=4)) +
  xlab("Sample ID") + ylab("Total reads per sample")

# saving and finishing up
cat ("Saving PCA plot as ", save_filename, " now.\n")
pdf(file = save_filename, width=9, height=12)
grid.arrange(org_relative_ggplot, org_absolute_ggplot, ncol = 1)
dev.off()