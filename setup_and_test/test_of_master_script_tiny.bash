#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 2000

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# test_of_master_script_tiny.bash
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified by Michelle Treiber on September 28, 2017
#
####################################################################
#
# This script was created to test the SAMSA2 pipeline, up to the DIAMOND
# annotation step. It uses TINY databases and will not give valid results
# files. To download the full databases, run full_database_download.bash
# that can be accessed at https://github.com/transcript/samsa2/tree/master/setup
# and run samsa2_master_script.bash located at
# https://github.com/transcript/samsa2/tree/master/sample_files_paired-end.
#
# The steps in this script are:
#	1.   Merging with PEAR, if applicable
#	2.   Read cleaning with Trimmomatic
#	3.   rRNA removal with SortMeRNA
#	4.   Annotation using DIAMOND (by default against the RefSeq database)
#   4.1  Annotation using DIAMOND against the Subsystems database
#
# NOTE: BEFORE running this script, please download and run package_installation.bash
# that can be accessed at https://github.com/transcript/samsa2/tree/master/setup
# in order to set up SAMSA2 dependencies.
#
####################################################################
#

set -e # terminate script after first thing fails

echo "WARNING: THIS IS A TEST RUN."
echo "You are running an example script to test SAMSA2 using TINY databases."
echo "To run real metatranscriptomes, please download one or both of the databases by running full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup."
echo -e "NOTE: Before running this script, please download and run package_installation.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"

####################################################################
#
# VARIABLES
#
# 	0. Starting location- Set pathway to location of samsa2 GitHub download
source "${BASH_SOURCE%/*}/../bash_scripts/common.sh"

#	1. Starting files location
starting_files_location=$SAMSA/sample_files_paired-end/1_starting_files

#	2. Output files location
mkdir $SAMSA/output_test
output_location=$SAMSA/output_test

#	3. Python scripts location
python_programs=$SAMSA/python_scripts

#	4. PEAR
pear_location=$SAMSA/programs/pear-0.9.10-linux-x86_64/bin

# 	5. Trimmomatic
trimmomatic_location=$SAMSA/programs/Trimmomatic-0.36

#	6. SortMeRNA
sortmerna_location=$SAMSA/programs/sortmerna-2.1

#	7. DIAMOND
diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"
diamond_location="$SAMSA/programs"

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#		Example: control_1.R1.fastq
#				 control_1.R2.fastq
#
# Note: if using single-end sequencing, skip this step (comment out).

#for file in $starting_files_location/*.gz
#do
#	gunzip $file
#done

cd $starting_files_location
for file in $starting_files_location/*R1*
do
	file1=$file
	file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'`
	out_path=`echo $file | awk -F "_R1" '{print $1 ".merged"}'`
	out_name=`echo ${out_path##*/}`
	$pear_location/pear -f $file1 -r $file2 -o $out_name
done

mkdir $output_location/step_1_output_test/
mv $starting_files_location/*merged* $output_location/step_1_output_test/
echo -e "\nPaired-end merging step completed.\n"

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $output_location/step_1_output_test/ folder!

for file in $output_location/step_1_output_test/*merged*
do
	shortname=`echo $file | awk -F "merged" '{print $1 "cleaned.fastq"}'`

	java -jar $trimmomatic_location/trimmomatic-0.36.jar SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

mkdir $output_location/step_2_output_test/
mv $output_location/step_1_output_test/*cleaned.fastq $output_location/step_2_output_test/
echo -e "\nCleaning files with Trimmomatic completed.\n"

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [ -f $output_location/step_2_output_test/raw_counts.txt ]
then
	rm $output_location/step_2_output_test/raw_counts.txt
	touch $output_location/step_2_output_test/raw_counts.txt
else
	touch $output_location/step_2_output_test/raw_counts.txt
fi

for file in $output_location/step_2_output_test/*cleaned.fastq
do
	python $python_programs/raw_read_counter.py -I $file -O $output_location/step_2_output_test/raw_counts.txt
done
echo -e "\nCounting raw sequences completed!\n"

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $output_location/step_2_output_test/*cleaned.fastq
do
	shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`

	$sortmerna_location/sortmerna --ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db --reads $file --aligned $file.ribosomes --other $shortname --fastx --num_alignments 0 --log -v

done

mkdir $output_location/step_3_output_test/
mv $output_location/step_2_output_test/*ribodepleted* $output_location/step_3_output_test/

echo -e "\nRibosomal read removal step completed!\n"

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $output_location/step_3_output_test/*ribodepleted.fastq
do
	shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
	echo "Now starting on " $file
	echo "Converting to " $shortname

	$diamond_location/diamond blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1
	$diamond_location/diamond view --daa $file.RefSeq.daa -o $shortname -f tab
done

mkdir $output_location/step_4_output_test/
mkdir $output_location/step_4_output_test/daa_binary_files/

mv $output_location/step_3_output_test/*annotated* $output_location/step_4_output_test/
mv $output_location/step_3_output_test/*.daa $output_location/step_4_output_test/daa_binary_files/

echo -e "\nRefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $output_location/step_3_output_test/*ribodepleted.fastq
do
	shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
	echo "Now starting on Subsystems annotations for " $file

	$diamond_location/diamond blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1
	$diamond_location/diamond view --daa $file.Subsys.daa -o $shortname -f tab
done

mv $output_location/step_3_output_test/*subsys_annotated* $output_location/step_4_output_test/
mv $output_location/step_3_output_test/*.daa $output_location/step_4_output_test/daa_binary_files/

echo -e "\nDIAMOND Subsystems annotations completed at: "; date

echo "Now removing output files."
rm -r $output_location
echo "SUCCESS! test_of_master_script.bash successfully ran."
exit
##################################################################
