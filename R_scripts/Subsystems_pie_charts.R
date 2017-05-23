# Subsystems_pie_charts.R
# Created 3/03/2017, by Sam Westreich
# Last updated 3/03/2017

library("ggplot2")
library("data.table")

setwd("~/Desktop/Projects/Lab Stuff/Aim 3/subsystems_results/reduced_files/")

# get list of files
files_list <- list.files( pattern = "*.reduced", full.names = T, recursive = FALSE)
file_names = ""
for (name in files_list) {
  file_names <- c(file_names, strsplit(name, split='./')[2])}
file_names <- file_names[-1]
file_names_trimmed = ""
for (name in file_names) {
  file_names_trimmed <- c(file_names_trimmed, unlist(strsplit(name, split='.', fixed=TRUE))[1])}
file_names_trimmed <- file_names_trimmed[-1]

# loading the files in
y <- 0
for (x in files_list) {
  y <- y + 1
  if (y == 1) {
    data_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(data_table) = c("DELETE", x, "Level4", "Level3", "Level2", "Level1")
    data_table <- data_table[,-1]
    rownames(data_table) <- data_table$Level4 }
  if (y > 1) {
    temp_table <- read.table(file = x, header = F, quote = "", sep = "\t", fill = TRUE)
    colnames(temp_table) = c("DELETE", x, "Level4", "Level3", "Level2", "Level1")
    rownames(temp_table) <- temp_table$Level4
    temp_table <- temp_table[,c(2,3)]
    data_table <- merge(temp_table, data_table, by = "Level4")  
    }
}
data_table <- data_table[,-ncol(data_table)]

# At this point, the whole table is read in.  Next step (for graphs) is to consolidate it
# down to averages.
data_table$Average <- rowSums(data_table[,c(1:length(files_list)+1)]) / length(files_list)
data_table <- data_table[order(-data_table$Average),]
# remove blank spots (no hierarchy)
data_table <- data_table[-which(data_table$Level1 == ""), ]

# Now, the next block run will change based on which level is wanted for the final pie chart:
# Starting with level 1
l1_table <- data.table(data_table[,c("Average","Level1")])  #NOTE: Change this for different levels
l1_table <- l1_table[, lapply(.SD, sum), by=Level1]
l1_table <- l1_table[order(-l1_table$Average)]

# Pie chart time?
CbPalette <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
      "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7",  
      "#ffffb3",  "#bebada",  "#fb8072",  "#80b1d3",  "#fdb462",  "#b3de69", 
      "#fccde5",  "#d9d9d9",  "#bc80bd",  "#ccebc5",  "#ffed6f", "#e41a1c",  
      "#377eb8", "#4daf4a",  "#984ea3",  "#ff7f00", "#ffff33",  "#a65628",  
      "#f781bf", "#999999", "#000000", "#a6cee4", "#1f78b5", "#b2df8b")

bp<- ggplot(l1_table[c(0:10),], aes(x="", y=Average, fill=Level1))+
  geom_bar(width = 1, stat = "identity") +
  scale_fill_manual(values = CbPalette) +
  ggtitle("SEED Subsystems hierarchy")

pie <- bp + coord_polar("y", start=0) +
#    geom_text(aes(x=1, label=Level1)) +
    theme(axis.text=element_blank(), legend.position = "right" ) 

cat ("Saving Subsystems pie chart as ", save_filename, " now.\n")
pdf(file = paste("Subsystems_pie_chart", ".pdf", sep = ""), width=10, height=7)
pie
dev.off()