# install_DESeq.R
# Created 9/21/2017, by Michelle Treiber

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2",lib=".")

packageurl1 <- "https://cran.r-project.org/src/contrib/optparse_1.4.4.tar.gz"
install.packages(packageurl1, repos=NULL, type="source", lib=".")

packageurl2 <- "https://cran.r-project.org/src/contrib/data.table_1.10.4.tar.gz"
install.packages(packageurl2, repos=NULL, type="source", lib=".")
