# volcano_plot.R
# Created 7/18/20
# Last updated 9/07/2020
# Run with --help flag for help.
#
# Need installation of EnhancedVolcano package beforehand
# install.packages('EnhancedVolcano')

suppressPackageStartupMessages({
  library(optparse)
})

option_list = list(
  make_option(c("-G", "--gtitle"), type="character", default="Control vs. Experimental",
              help="title to be displayed on graph", metavar="character"),
  make_option(c("-O", "--out"), type="character", default="volcano_plot.pdf",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-D", "--directory"), type="character", default="./",
              help="working directory location", metavar="character"),
  make_option(c("-F", "--fccutoff"), type="numeric", default=1.50,
              help="log2 Fold Change cut-off; default is 1.5", metavar="character"),
  make_option(c("-P", "--pcutoff"), type="numeric", default=0.05,
              help="P-adjusted cut-off; default is 0.05", metavar="character"),
  make_option(c("-C", "--connect"), type="logical", default=TRUE,
              help="label for significantly expressed data points; default is TRUE", metavar="character")
);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("USAGE: $ make_volcano_plot.R -D working_directory/ -O save.filename -G title -F fc_cutoff -P p_cutoff -C TRUE")

# check for necessary specs
if (is.null(opt$directory)) {
  print ("WARNING: No working directory specified with '-d' flag.")
  stop()
} else {
  cat ("Working directory is ", opt$directory, "\n")
  wd_location <- opt$directory
  setwd(wd_location)
}

if (is.null(opt$out)) {
  print ("WARNING: No save name for plot specified; defaulting to 'volcano_plot.pdf'.")
  save_filename <- "volcano_plot.pdf"
} else {
  cat ("Saving results as ", opt$out, "\n")
  save_filename <- opt$out
}

# importing DESeq_results table
DESeq_results <- list.files(
  pattern = "*DESeq_results.tab", full.names = T, recursive = FALSE)

for (x in DESeq_results) {
  res <- read.table(file = x, header = TRUE, sep = "\t")
} 

# reformatting the DESeq_results.tab
# head(res)
simpres <- res[,-8]; 
simpres <- simpres[,-7]; 
simpres <- simpres[,-6]; 
simpres <- simpres[,-4]; 
simpres <- simpres[,-3]; 
simpres <- simpres[,-2]
# str(simpres)
# summary(simpres)

# loading package
suppressPackageStartupMessages({
  library(EnhancedVolcano)
})

# plotting
volcano_plot <- EnhancedVolcano(
  toptable = simpres,
  lab = simpres[,1],
  x = 'log2FoldChange',
  y = 'padj',
  title = opt$gtitle,
  ylim = c(0, max(-log10(simpres$padj), na.rm=TRUE) + 2),
  xlim = c(-5, 5),
  pCutoff = opt$pcutoff,
  FCcutoff = opt$fccutoff,
  drawConnectors = opt$connect,
  widthConnectors = 0.5,
  lengthConnectors = unit(0.005, 'npc'),
  endsConnectors = 'last',
  colConnectors = 'black',
  legendPosition = 'right',
  legendLabels=c('NS','Log2 FC','p-value',
                 'p-value & Log2 FC')
)



# saving and finishing up
cat ("Saving volcano plot as ", save_filename, " now.\n")
pdf(file = save_filename, width=12, height=9)
volcano_plot
dev.off()
