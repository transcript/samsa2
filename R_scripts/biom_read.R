# biom_read.R
# Created 7/07/2016
# Turns BIOM format objects into something actually useful.

library(biom)
library(optparse)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Input file in BIOM format", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="parsed_biom_file.tab", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-d", "--deseq_out"), type="character", default="biom_DESeq_results",
              help="DESeq output file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)) {
  print ("WARNING: No input file specified with '-i' flag.")
  stop()
} else {
  cat ("Input file name is is ", opt$input, "\n")
  infile <- opt$input
}

if (is.null(opt$out)) {
  print ("WARNING: No save name for plot specified; defaulting to 'parsed_biom_file.tab'.") }

# handling the BIOM input file
biom_input <- read_biom(infile)
t1 <- biom_data(biom_input)
t2 <- as.data.frame(as.matrix(t1))

# reordering (NOTE: ONLY works for macaque study, will need to be edited for other purposes!!)
t2 <- t2[,c(2,4,9,10,11,12,15,17,18,21,23,24,1,3,5,6,7,8,13,14,16,19,20,22)]

# ordering by overall abundance
t2$Total <- rowSums(t2)
t2 <- t2[order(-t2$Total),]
t2$Total <- NULL

# DESeq stats time!
library(DESeq2)
t3 <- t2
condition <- data.frame(condition=factor(c(rep("control", 12), rep("ICD", 12))))
colnames(t3) <- condition$condition
dds <- DESeqDataSetFromMatrix(t2, condition, ~ condition)
dds <- DESeq(dds)

baseMeanPerLvl <- sapply( levels(dds$condition), 
    function(lvl) rowMeans( counts(dds,normalized=TRUE)[,dds$condition == lvl] ) )

res <- results(dds)
seed_results <- data.frame(res)
seed_results <- merge(seed_results, baseMeanPerLvl, by="row.names")
seed_results <- seed_results[,c(1,2,8,9,3,4,5,6,7)]
colnames(seed_results)[c(3,4)] <- c("controlMean", "experimentalMean")
seed_results <- seed_results[order(-seed_results$baseMean),]

# output upregulated and downregulated results as two different files
upregulated_results <- subset(seed_results, experimentalMean > controlMean)
downregulated_results <- subset(seed_results, experimentalMean < controlMean)

# saving DESeq results
cat ("Saving DESeq output (two tables, one up, one down) as ", opt$deseq_out, 
     "_[upregulated/downregulated].tsv now.\n", sep = "")
write.table(upregulated_results, file=paste(opt$deseq_out, "_upregulated.tsv", sep=""),
    append=FALSE, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
write.table(downregulated_results, file=paste(opt$deseq_out, "_downregulated.tsv", sep=""),
            append=FALSE, quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

# saving DESeq output
#cat ("Saving DESeq output tables as ", opt$deseq_out, " now.\n")
#write.table(seed_results, file=opt$deseq_out, append=FALSE, quote=FALSE, sep="\t", 
#        row.names=FALSE, col.names=TRUE)

# results output as a table (.tab)
#cat ("Saving BIOM output data table as ", opt$out, " now.\n")
#write.table(t2, file=opt$out, append=FALSE, quote=FALSE, sep="\t", row.names=TRUE, 
#            col.names=TRUE)
