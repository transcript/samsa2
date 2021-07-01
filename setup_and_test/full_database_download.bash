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
# Updated 1 July 2021 by Sam Westreich - files are now hosted on Zenodo, rather than Bioshare
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
# Set pathway for SAMSA to location of samsa2 GitHub download:
source "${BASH_SOURCE%/*}/../bash_scripts/lib/common.sh"

# Make database directory
mkdir $SAMSA/full_databases
cd $SAMSA/full_databases

echo -e "NOTE: The databases are up to 28GB and may require hours to download. Users may want to consider running this download overnight.\n"
# Download NCBI RefSeq database:
echo "NOW DOWNLOADING NCBI REFSEQ DATABASE AT: "; date
wget "https://zenodo.org/record/5022377/files/RefSeq_bac.fa.bz2" --no-check-certificate

# Download SEED Subsystems database:
echo "NOW DOWNLOADING SEED SUBSYSTEMS DATABASE AT: "; date
wget "https://zenodo.org/record/5022377/files/subsys_db.fa.bz2" --no-check-certificate

# unzipping
bunzip2 *.bz2

echo -e "\n\nNOTE: IF USERS ARE USING A DIFFERENT VERSION OF DIAMOND OR
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
wget "https://zenodo.org/record/5022377/files/RefSeq_bac.dmnd.bz2" --no-check-certificate

# Download DIAMOND compatible Subsystems database:
echo "NOW DOWNLOADING DIAMOND COMPATIBLE SUBSYSTEMS DATABASE AT: "; date
wget "https://zenodo.org/record/5022377/files/subsys_db.dmnd.bz2" --no-check-certificate

# unzipping
bunzip2 *.bz2

echo "Completed!"
exit
