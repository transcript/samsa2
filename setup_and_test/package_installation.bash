#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --cpus-per-task 10
#SBATCH --mem-per-cpu 2000

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# package_installation.bash
# Created September 28, 2017 by Michelle Treiber
#
####################################################################
#
# This script was created to install all SAMSA2 package dependencies.
# 
# The packages are as follows:
# 1. R packages:
#		a. DESeq2
#		b. getopt
#		c. optparse
#		d. data.table
# 2. PEAR
# 3. Trimmomatic
# 5. DIAMOND
# 4. Sortmerna
#
####################################################################
#
# Set pathway for SAMSA to location of samsa2 GitHub download:
IGNORE_DEPS=1 source "${BASH_SOURCE%/*}/../bash_scripts/lib/common.sh"

cd $PROGRAMS

# 1.  Unpack PEAR package:
echo "Now extracting PEAR at: "; date
tar -xvzf ${PEAR_DIR}.tar.gz
echo "PEAR extraction finished at: "; date

# 2. Unzip Trimmomatic package:
echo "Now extracting Trimmomatic at: "; date
unzip ${TRIMMOMATIC_DIR}.zip
echo "Trimmomatic extraction finished at: "; date

# 3. Unpack DIAMOND package
echo "Now extracting DIAMOND at: "; date
tar -xvzf $PROGRAMS/diamond-linux64.tar.gz
echo "DIAMOND extraction finished at: "; date

# 4. Unpack and build sortmerna package:
echo "Now extracting SortMeRNA at: "; date
tar -xvzf ${SORTMERNA_DIR}.tar.gz
cd $SORTMERNA_DIR
bash ./build.sh
echo "SortMeRNA extraction finished at: "; date

cd ../

# 5. Download necessary R packages:
echo "Now extracting R packages at: "; date
mkdir $SAMSA/R_scripts/packages
cd $SAMSA/R_scripts/packages
# If using R version 3.6, change this file in the line below to install_packages_R_3.6.R .
R --save < ../install_packages.R -v
echo "R package extraction finished at: "; date

# Index silva-bac-16s-id90.fasta for sortmerna use:
echo "Indexing the SILVA rRNA database for SortMeRNA; feel free to abort at this step if this isn't needed."
$SORTMERNA_DIR/indexdb_rna --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db

echo "Finished extracting and installing all SAMSA2 package dependencies at: "; date
echo "Completed!"
exit
