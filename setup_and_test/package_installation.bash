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
# Set pathway for starting_location to location of samsa2 GitHub download:
starting_location=/home/samsa2

# 1. Download necessary R packages:
echo "Now downloading R packages at: "; date
mkdir $starting_location/R_scripts/packages
cd $starting_location/R_scripts/packages
R --save < ../install_packages.R -v
echo "R package downloads finished at: "; date

cd $starting_location/programs

# 2.  Unpack PEAR package:
echo "Now downloading PEAR at: "; date
tar -xvzf $starting_location/programs/pear-0.9.10-linux-x86_64.tar.gz
echo "PEAR download finished at: "; date

# 3. Unzip Trimmomatic package:
echo "Now downloading Trimmomatic at: "; date
unzip $starting_location/programs/Trimmomatic-0.36.zip
echo "Trimmomatic download finished at: "; date

# 4. Unpack DIAMOND package
echo "Now downloading diamond at: "; date
tar -xvzf $starting_location/programs/diamond-linux64.tar.gz
echo "DIAMOND download finished at: "; date

# 5. Unpack and build sortmerna package:
echo "Now downloading sortmerna at: "; date
tar -xvzf $starting_location/programs/sortmerna-2.1.tar.gz
cd $starting_location/programs/sortmerna-2.1
bash ./build.sh
echo "Sortmerna download finished at: "; date

# Index silva-bac-16s-id90.fasta for sortmerna use:
sortmerna_location=$starting_location/programs/sortmerna-2.1
$sortmerna_location/indexdb_rna --ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db

echo "Finished downloading and installing all SAMSA2 package dependencies at: "; date
echo "Completed!"