#!/bin/bash
#SBATCH --mem=100000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user=stwestreich@ucdavis.edu
#SBATCH --mail-type=END

####################################################################
#
# master_script.sh
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified June 16, 2017
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#	1.   Merging with PEAR, if applicable
#	2.   Read cleaning with Trimmomatic
#	3.   rRNA removal with SortMeRNA
#	4.   Annotation using DIAMOND (by default against the RefSeq database)
#	5.   Aggregation using analysis_counter.py
#   4.1  Annotation using DIAMOND against the Subsystems database
#   5.1  Aggregation using Subsystems-specific analysis counter.py
#
####################################################################
#
# VARIABLES - set these paths for each step.
#
#	0. Starting files location
starting_location=/share/milklab/sam/test_files

#	1. PEAR
pear_location=/home/swestreich/programs/pear-0.9.6

# 	2. Trimmomatic
trimmomatic_location=/home/swestreich/programs/Trimmomatic-0.33

#	3. SortMeRNA
sortmerna_location=/home/swestreich/programs/sortmerna-2.1-linux-64

#	4. DIAMOND
diamond_database="/share/milklab/sam/databases/bct_new"
diamond_subsys_db="/share/milklab/sam/databases/subsys_new"
diamond_location="/home/swestreich/programs"

#	5. Aggregation
python_programs=/share/milklab/sam/python_scripts

#	6. R scripts and paths
export R_LIBS="/share/milklab/sam/R_scripts/packages"
R_programs=/share/milklab/sam/R_scripts

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
# Note: if using single-end sequencing, skip this step (comment out).
# Note: if performing R analysis (step 6), be sure to name files with 
# 	the appropriate prefix ("control_$file" and "experimental_$file")!

for file in $starting_location/*R1*
do
	file1=$file
	file2=`echo $file1 | awk -F"R1" '{print $1 "R2" $2}'`
	out_path=`echo $file | awk -F"R1" '{print $1 "merged"}'`
	out_name=`echo ${out_path##*/}`

	$pear_location/pear-0.9.6 -f $file1 -r $file2 -o $out_name
done

mkdir $starting_location/step_1_output/
mv $starting_location/*merged* $starting_location/step_1_output/

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping FLASH, make sure that all starting files are in the
# $starting_location/step_1_output/ folder!

for file in $starting_location/step_1_output/*.assembled*
do
	shortname=`echo $file | awk -F"assembled" '{print $1 "cleaned.fastq"}'`

	java -jar $trimmomatic_location/trimmomatic-0.33.jar SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

mkdir $starting_location/step_2_output/
mv $starting_location/step_1_output/*cleaned.fastq $starting_location/step_2_output/

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $starting_location/step_2_output/*.cleaned.fastq
do
	unzip_name=`echo $file | awk -F".gz" '{print $1}'`
	shortname=`echo $file | awk -F"fna" '{print $1 "ribodepleted"}'`

	$sortmerna_location/sortmerna --ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db --reads $file --aligned $file.ribosomes --other $shortname --fastx --num_alignments 0 --log -v

done

mkdir $starting_location/step_3_output/
mv $starting_location/step_2_output/*ribodepleted $starting_location/step_3_output/

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $starting_location/step_3_output/*ribodepleted
do
	shortname=`echo $file | awk -F "ribodepleted" '{print $1 "annotated"}'`
	echo "Now starting on " $file 
	echo "Converting to " $shortname

	$diamond_location/diamond blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1 --sensitive
	$diamond_location/diamond view --daa $file.dmd -o $shortname -f tab
done

mkdir $starting_location/step_4_output/
mkdir $starting_location/step_4_output/daa_binary_files/

mv $starting_location/step_3_output/*annotated* $starting_location/step_4_output/
mv $starting_location/step_3_output/*.daa $starting_location/step_4_output/daa_binary_files/

echo "RefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $starting_location/step_4_output/*annotated*
do
	python $python_programs/DIAMOND_analysis_counter.py -I $file -D $diamond_database -O
	python $python_programs/DIAMOND_analysis_counter.py -I $file -D $diamond_database -F
done

mkdir $starting_location/step_5_output/RefSeq_results/org_results/
mkdir $starting_location/step_5_output/RefSeq_results/func_results/
mv $starting_location/step_4_output/*organism.tsv $starting_location/step_5_output/RefSeq_results/org_results/
mv $starting_location/step_4_output/*function.tsv $starting_location/step_5_output/RefSeq_results/func_results/

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $starting_location/step_3_output/*ribodepleted
do
	shortname=`echo $file | awk -F"ribodepleted" '{print $1 "subsys_annotated"}'`
	echo "Now starting on Subsystems annotations for " $file

	$diamond_location/diamond blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1 --sensitive
	$diamond_location/diamond view --daa $file.dmd -o $shortname -f tab
done

mv $starting_location/step_3_output/*subsys_annotated* $starting_location/step_4_output/
mv $starting_location/step_3_output/*.daa $starting_location/step_4_output/daa_binary_files/

echo "DIAMOND Subsystems annotations completed at: "; date

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $starting_location/step_4_output/*subsys_annotated*
do
	python $python_programs/DIAMOND_subsystems_analysis_counter.py -I $file -D $diamond_subsys_db.fa -P $file.receipt
	
	# This quick program reduces down identical hierarchy annotations
	python $python_programs/subsys_reducer.py -I $file.hierarchy
done

mkdir $starting_location/step_5_output/Subsystems_results/
mkdir $starting_location/step_5_output/Subsystems_results/receipts/
mv $starting_location/step_4_output/*.reduced $starting_location/step_5_output/Subsystems_results/
mv $starting_location/step_4_output/*.receipt $starting_location/step_5_output/Subsystems_results/receipts/
rm $starting_location/step_4_output/*.hierarchy

##################################################################
#
# At this point, all the results files are ready for analysis using R.
# This next step performs basic DESeq2 analysis of the RefSeq organism, function,
# and Subsystems annotations.
#
# More complex R analyses may be performed using specific .sh analysis scripts.
#
# STEP 6: R ANALYSIS
# Note: For R to properly identify files to compare/contrast, they must include 
# the appropriate prefix (either "control_$file" or experimental_$file")!

Rscript $R_programs/run_DESeq_stats.R -I $starting_location/step_5_output/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab
Rscript $R_programs/run_DESeq_stats.R -I $starting_location/step_5_output/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab
Rscript $R_programs/Subsystems_DESeq_stats.R -I $starting_location/step_5_output/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1

echo "Master bash script finished running at: "; date
exit
####################################################################
