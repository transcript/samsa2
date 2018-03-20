#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH --cpus-per-task 20
#SBATCH --mem-per-cpu 2000

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script_for_sample_files.bash
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified by Michelle Treiber on August 9, 2017
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
#	6.	 Running R scripts to get DESeq statistical analysis.
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set 
# up SAMSA2 dependencies and download full databases.
#
####################################################################
#
echo -e "NOTE: Before running this script, user must run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - Set pathway for starting_location to location of samsa2 GitHub download
#
# 	0. Set starting location:
starting_location=/home/samsa2

#	00. Starting files location
starting_files_location=$starting_location/sample_files_paired-end/1_starting_files

#	1. PEAR
pear_location=$starting_location/programs/pear-0.9.10-linux-x86_64/bin

# 	2. Trimmomatic
trimmomatic_location=$starting_location/programs/Trimmomatic-0.36

#	3. SortMeRNA
sortmerna_location=$starting_location/programs/sortmerna-2.1

#	4. DIAMOND
diamond_database="$starting_location/full_databases/RefSeq_bac"
diamond_subsys_db="$starting_location/full_databases/subsys_db"
diamond_location="$starting_location/programs/diamond"

#	5. Aggregation
python_programs=$starting_location/python_scripts
RefSeq_db="$starting_location/full_databases/RefSeq_bac.fa"
Subsys_db="$starting_location/full_databases/subsys_db.fa"

#	6. R scripts and paths
export R_LIBS="$starting_location/R_scripts/packages"
R_programs=$starting_location/R_scripts

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#		Example: control_1.R1.fastq
#				 control_1.R2.fastq
#
# Note: if using single-end sequencing, skip this step (comment out).
# Note: if performing R analysis (step 6), be sure to name files with 
# 	the appropriate prefix ("control_$file" and "experimental_$file")!

#for file in $starting_files_location/*.gz
#do
#	gunzip $file
#done

cd $starting_files_location
for file in $starting_files_location/*R1*
do
	file1=$file
	file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'`
	out_path=`echo $file | awk -F "R1" '{print $1 "merged"}'`
	out_name=`echo ${out_path##*/}`
	$pear_location/pear -f $file1 -r $file2 -o $out_name
done

mkdir $starting_files_location/step_1_output/
mv $starting_files_location/*merged* $starting_files_location/step_1_output/
echo "STEP 1 DONE"

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $starting_files_location/step_1_output/ folder!

for file in $starting_files_location/step_1_output/*assembled*
do
	shortname=`echo $file | awk -F "merged" '{print $1 "cleaned.fastq"}'`

	java -jar $trimmomatic_location/trimmomatic-0.36.jar SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

mkdir $starting_files_location/step_2_output/
mv $starting_files_location/step_1_output/*cleaned.fastq $starting_files_location/step_2_output/
echo "STEP 2 DONE"

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [ -f $starting_files_location/step_2_output/raw_counts.txt ]
then
	rm $starting_files_location/step_2_output/raw_counts.txt
	touch $starting_files_location/step_2_output/raw_counts.txt
else
	touch $starting_files_location/step_2_output/raw_counts.txt
fi	

for file in $starting_files_location/step_2_output/*cleaned.fastq
do
	python $python_programs/raw_read_counter.py -I $file -O $starting_files_location/step_2_output/raw_counts.txt
done
echo "STEP 2.9 DONE" 

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $starting_files_location/step_2_output/*cleaned.fastq
do
	shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`

	$sortmerna_location/sortmerna --ref $sortmerna_location/rRNA_databases/silva-bac-16s-id90.fasta,$sortmerna_location/index/silva-bac-16s-db --reads $file --aligned $file.ribosomes --other $shortname --fastx --num_alignments 0 --log -v

done

mkdir $starting_files_location/step_3_output/
mv $starting_files_location/step_2_output/*ribodepleted* $starting_files_location/step_3_output/

echo "STEP 3 DONE"

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $starting_files_location/step_3_output/*ribodepleted.fastq
do
	shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
	echo "Now starting on " $file 
	echo "Converting to " $shortname

	$diamond_location blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1 --sensitive
	$diamond_location view --daa $file.RefSeq.daa -o $shortname -f tab
done

mkdir $starting_files_location/step_4_output/
mkdir $starting_files_location/step_4_output/daa_binary_files/

mv $starting_files_location/step_3_output/*annotated* $starting_files_location/step_4_output/
mv $starting_files_location/step_3_output/*.daa $starting_files_location/step_4_output/daa_binary_files/

echo "RefSeq DIAMOND annotations completed at: "; date
echo "STEP 4 DONE"

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $starting_files_location/step_4_output/*RefSeq_annotated*
do
	python $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
	python $python_programs/standardized_DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
done

mkdir $starting_files_location/step_5_output/
mkdir $starting_files_location/step_5_output/RefSeq_results/
mkdir $starting_files_location/step_5_output/RefSeq_results/org_results/
mkdir $starting_files_location/step_5_output/RefSeq_results/func_results/
mv $starting_files_location/step_4_output/*organism.tsv $starting_files_location/step_5_output/RefSeq_results/org_results/
mv $starting_files_location/step_4_output/*function.tsv $starting_files_location/step_5_output/RefSeq_results/func_results/

echo "STEP 5 DONE"

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $starting_files_location/step_3_output/*ribodepleted.fastq
do
	shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
	echo "Now starting on Subsystems annotations for " $file

	$diamond_location blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1 --sensitive
	$diamond_location view --daa $file.Subsys.daa -o $shortname -f tab
done

mv $starting_files_location/step_3_output/*subsys_annotated* $starting_files_location/step_4_output/
mv $starting_files_location/step_3_output/*.daa $starting_files_location/step_4_output/daa_binary_files/

echo "DIAMOND Subsystems annotations completed at: "; date
echo "STEP 4.1 DONE"

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $starting_files_location/step_4_output/*subsys_annotated
do
	python $python_programs/DIAMOND_subsystems_analysis_counter.py -I $file -D $Subsys_db -O $file.hierarchy -P $file.receipt
	
	# This quick program reduces down identical hierarchy annotations
	python $python_programs/subsys_reducer.py -I $file.hierarchy
done

mkdir $starting_files_location/step_5_output/Subsystems_results/
mkdir $starting_files_location/step_5_output/Subsystems_results/receipts/
mv $starting_files_location/step_4_output/*.reduced $starting_files_location/step_5_output/Subsystems_results/
mv $starting_files_location/step_4_output/*.receipt $starting_files_location/step_5_output/Subsystems_results/receipts/
rm $starting_files_location/step_4_output/*.hierarchy

echo "STEP 5.1 DONE"

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

Rscript $R_programs/run_DESeq_stats.R -I $starting_files_location/step_5_output/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab -R $starting_files_location/step_2_output/raw_counts.txt
Rscript $R_programs/run_DESeq_stats.R -I $starting_files_location/step_5_output/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab -R $starting_files_location/step_2_output/raw_counts.txt
Rscript $R_programs/Subsystems_DESeq_stats.R -I $starting_files_location/step_5_output/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1 -R $starting_files_location/step_2_output/raw_counts.txt

echo "STEP 6 DONE"
echo "Master bash script finished running at: "; date
exit
####################################################################
