# install_DESeq.R
# Created 9/21/2017, by Michelle Treiber

## try http:// if https:// URLs are not supported
## This command works for up to R version 3.5.  If using R 3.6, comment out this section and uncomment the following block.
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2", lib=".")

## Un-comment this stuff if you are using R version 3.6
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DESeq2")

packageurl1 <- "https://cran.r-project.org/src/contrib/Archive/getopt/getopt_1.20.0.tar.gz"
install.packages(packageurl1, repos=NULL, type="source", lib=".")

packageurl2 <- "https://cran.r-project.org/src/contrib/Archive/optparse/optparse_1.4.4.tar.gz"
install.packages(packageurl2, repos=NULL, type="source", lib=".")

packageurl3 <- "https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.10.4.tar.gz"
install.packages(packageurl3, repos=NULL, type="source", lib=".")

# installing EnhancedVolcano; this is pulled from GitHub as the CRAN version isn't compatible with older bioconductor versions
devtools::install_github("kevinblighe/EnhancedVolcano")
