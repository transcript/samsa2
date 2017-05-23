# Subsystems_DESeq_stats.R
# Created 3/06/2017, by Sam Westreich
# Last updated 3/06/2017

library("DESeq2")
library("data.table")

setwd("~/path/to/files")

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

# loading the control files in
y <- 0
for (x in control_files) {
  y <- y + 1
  if (y == 1) {
    control_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(control_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    control_table <- control_table[,-1]
    rownames(control_table) <- control_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    control_table <- merge(temp_table, control_table, by = "Level4")  
  }
}
control_table <- control_table[,-ncol(control_table)]

# loading the experimental files in
y <- 0
for (x in exp_files) {
  y <- y + 1
  if (y == 1) {
    exp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(exp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    exp_table <- exp_table[,-1]
    rownames(exp_table) <- exp_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level1", "Level2")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    exp_table <- merge(temp_table, exp_table, by = "Level4")  
  }
}
exp_table <- exp_table[,-ncol(exp_table)]

# At this point, the whole table is read in.  Next step (for statistical comparison) is to 
# get just the level we want to compare.

# for Level 1 comparisons:
l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level2", "Level3", "Level4")])
l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level2", "Level3", "Level4")])
# OR for level 2 comparisons
l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level3", "Level4")])
l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level3", "Level4")])
names(l1_control_table)[names(l1_control_table) == 'Level2'] <- 'Level1'
names(l1_exp_table)[names(l1_exp_table) == 'Level2'] <- 'Level1'
# OR for level 3 comparisons
l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level2", "Level4")])
l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level2", "Level4")])
names(l1_control_table)[names(l1_control_table) == 'Level3'] <- 'Level1'
names(l1_exp_table)[names(l1_exp_table) == 'Level3'] <- 'Level1'
# OR for level 4 comparisons
l1_control_table <- data.table(control_table[, !names(control_table) %in% c("Level1", "Level3", "Level2")])
l1_exp_table <- data.table(exp_table[, !names(exp_table) %in% c("Level1", "Level3", "Level2")])
names(l1_control_table)[names(l1_control_table) == 'Level4'] <- 'Level1'
names(l1_exp_table)[names(l1_exp_table) == 'Level4'] <- 'Level1'
l1_control_table <- data.table(as.data.frame(l1_control_table)[,c(2:(ncol(l1_control_table)),1)])
l1_exp_table <- data.table(as.data.frame(l1_exp_table)[,c(2:(ncol(l1_exp_table)),1)])

# remove blank spots (no hierarchy)
l1_control_table <- l1_control_table[-which(l1_control_table$Level1 == ""), ]
l1_exp_table <- l1_exp_table[-which(l1_exp_table$Level1 == ""), ]

# reducing stuff down to avoid duplicates
colnames(l1_control_table) <- c(control_names_trimmed, "Level1")
colnames(l1_exp_table) <- c(exp_names_trimmed, "Level1")
l1_control_table <- l1_control_table[, lapply(.SD, sum), by=Level1]
l1_exp_table <- l1_exp_table[, lapply(.SD, sum), by=Level1]
l1_table <- merge(l1_control_table, l1_exp_table, by="Level1", all.x = T)
rownames(l1_table) <- l1_table$Level1
l1_names <- l1_table$Level1
l1_table$Level1 <- NULL
l1_table[is.na(l1_table)] <- 0

# now the DESeq stuff
completeCondition <- data.frame(condition=factor(c(rep("experimental", length(exp_files)), 
  rep("control", length(control_files)))))
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

write.table(l1_results, file = "Subsystems_level_4_DESeq_results.tab", 
  append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
