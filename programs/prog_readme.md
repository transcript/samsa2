# Programs readme file
******

The following programs are used in the SAMSA2 pipeline:

1. PEAR, a flexible paired-end read merger
2. Trimmomatic, a read cleanup and adaptor removal program
3. SortMeRNA, a program to detect and remove ribosomal sequences
4. DIAMOND, a superfast BLAST-like aligner
5. Several specific R packages.

Several of these dependencies are included here in this folder, when licenses allow.  
For files that cannot be provided here, linkes are given to provide access to download the program, along with details such as version number verified to work with SAMSA2.

## PEAR

The PEAR version used for SAMSA2 is 0.9.6.  

PEAR requires registration to download, and so cannot be hosted here; however, PEAR is free for academic use.

PEAR's homepage: [https://sco.h-its.org/exelixis/web/software/pear/](https://sco.h-its.org/exelixis/web/software/pear/)

PEAR's download page: [https://www.h-its.org/downloads/pear-academic/](https://www.h-its.org/downloads/pear-academic/)

## Trimmomatic

The Trimmomatic version used for SAMSA2 is 0.36.

Trimmomatic does not have any explicit license, but can be downloaded from the lab's page.

Trimmomatic's homepage: [http://www.usadellab.org/cms/?page=trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

## SortMeRNA

The SortMeRNA version used for SAMSA2 is 2.1.

SortMeRNA is available under the GNU 3.0 General License; a gzipped copy of SortMeRNA is provided in this directory.

SortMeRNA's Github page: [https://github.com/biocore/sortmerna](https://github.com/biocore/sortmerna)

## DIAMOND

The DIAMOND version used for SAMSA2 is 0.8.38 (NOTE: more recent versions of DIAMOND are available, but have not yet been tested with SAMSA2.  Older versions of DIAMOND, such as 0.7.9, will **not work.**)

DIAMOND is available under the AGPL GNU General License; a gzipped copy of DIAMOND is provided in this directory.

DIAMOND's Github page: [https://github.com/bbuchfink/diamond](https://github.com/bbuchfink/diamond)

## R packages

The R packages used by SAMSA2's analysis scripts include:

* optparse, MIT license
* DESeq, GNU General Public License v3.0
* data.table, GNU General Public License v3.0
