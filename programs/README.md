## Programs readme file
******

The following programs are used in the SAMSA2 pipeline:

1. PEAR, a flexible paired-end read merger
2. Trimmomatic, a read cleanup and adaptor removal program
3. SortMeRNA, a program to detect and remove ribosomal sequences
4. DIAMOND, a superfast BLAST-like aligner
5. Several specific R packages.

Several of these dependencies are included here in this folder, when licenses allow.  
For files that cannot be provided here, linkes are given to provide access to download the program, along with details such as version number verified to work with SAMSA2.

### PEAR

The PEAR version used for SAMSA2 is 0.9.6.  

PEAR is free for academic use under the Creative Commmons license; gzipped copies (for both x32 and x86 versions) are provided in this directory.

PEAR's homepage: [https://sco.h-its.org/exelixis/web/software/pear/](https://sco.h-its.org/exelixis/web/software/pear/)

### Trimmomatic

The Trimmomatic version used for SAMSA2 is 0.36.

Trimmomatic is available under the GNU 3.0 General License; a gzipped copy is provided in this directory.

Trimmomatic's homepage: [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

### SortMeRNA

The SortMeRNA version used for SAMSA2 is 2.1.

SortMeRNA is available under the GNU 3.0 General License; a gzipped copy of SortMeRNA is provided in this directory.

SortMeRNA's Github page: [https://github.com/biocore/sortmerna](https://github.com/biocore/sortmerna)

### DIAMOND

The DIAMOND version used for SAMSA2 is 0.8.38 (NOTE: more recent versions of DIAMOND are available, but have not yet been tested with SAMSA2.  Older versions of DIAMOND, such as 0.7.9, will **not work.**)

DIAMOND is available under the AGPL GNU General License; a gzipped copy of DIAMOND is provided in this directory.

DIAMOND's Github page: [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)

### R packages

The R packages used by SAMSA2's analysis scripts include:

* optparse, version 1.4.4, MIT license
* data.table, version 1.10.4, GNU General Public License v3.0
* DESeq2, version 1.16.1, GNU General Public License v3.0

In R, you can use the following commands to install the above-listed versions of these packages:

    install_version("optparse", version = "1.4.4", repos = "http://cran.us.r-project.org")
    install_version("data.table", version = "1.10.4", repos = "http://cran.us.r-project.org")

For DESeq2, as it is a BioConductor package, the above-listed version can be accessed at the following URL, under BioConductor version 3.5: [https://bioconductor.org/packages/3.5/bioc/html/DESeq2.html](https://bioconductor.org/packages/3.5/bioc/html/DESeq2.html)
