#!/bin/bash
#SBATCH -t 30:00:00
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 2000

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# full_database_download.bash
# Created September 28, 2017 by Michelle Treiber
#
####################################################################
#
# This script was created to download full databases for use in 
# the SAMSA2 master script.
#
# NOTE: The databases are up to 28GB and may require many hours to download.
# Users may want to consider running this download overnight.
# 
####################################################################
#
# Set pathway for starting_location to location of samsa2 GitHub download:
starting_location=/home/samsa2

# Make database directory
mkdir $starting_location/full_databases
cd $starting_location/full_databases

echo "NOTE: The databases are up to 28GB and may require many hours to download. Users may want to consider running this download overnight."
# Download NCBI RefSeq database:
echo "NOW DOWNLOADING NCBI REFSEQ DATABASE AT: "; date
wget "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/RefSeq_bac.fa"

# Download SEED Subsystems database:
echo "NOW DOWNLOADING SEED SUBSYSTEMS DATABASE AT: "; date
wget "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/subsys_db.fa"

echo "\n\nNOTE: IF USERS ARE USING A DIFFERENT VERSION OF DIAMOND OR
WOULD RATHER MAKE THEIR OWN DIAMOND COMPATIBLE DATABASES
TO SAVE TIME:
1. Run package_installation.bash located at 
		https://github.com/transcript/samsa2/tree/master/setup
		OR download and install other desired DIAMOND version 
		See https://github.com/bbuchfink/diamond for more details
2. Comment out the next 4 lines of code
3. Uncomment the last 2 lines of code and if using a different
		version of DIAMOND than the one on SAMSA2 Github, change 
		the diamond location\n\n"

# Download DIAMOND compatible RefSeq database:
echo "NOW DOWNLOADING DIAMOND COMPATIBLE REFSEQ DATABASE AT: "; date
wget "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/RefSeq_bac.dmnd"

# Download DIAMOND compatible Subsystems database:
echo "NOW DOWNLOADING DIAMOND COMPATIBLE SUBSYSTEMS DATABASE AT: "; date
wget "https://bioshare.bioinformatics.ucdavis.edu/bioshare/download/2c8s521xj9907hn/subsys_db.dmnd"

# $starting_location/programs/diamond makedb --in $starting_location/full_databases/RefSeq_bac.fa --db $starting_location/full_databases/RefSeq_bac
# $starting_location/programs/diamond makedb --in $starting_location/full_databases/subsys_db.fa --db $starting_location/full_databases/subsys_db

echo "Completed!"