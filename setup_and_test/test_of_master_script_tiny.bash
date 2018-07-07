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
#   1.   Merging with PEAR, if applicable
#   2.   Read cleaning with Trimmomatic
#   3.   rRNA removal with SortMeRNA
#   4.   Annotation using DIAMOND (by default against the RefSeq database)
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
# Starting location- Set pathway to location of samsa2 GitHub download
source "${BASH_SOURCE%/*}/../bash_scripts/common.sh"

# Starting files location
starting_files_location=$SAMSA/sample_files_paired-end/1_starting_files

# Output files location
OUT_DIR=$SAMSA/output_test
if [[ ! -d "$OUT_DIR" ]]; then
  mkdir "$OUT_DIR"
fi

# DIAMOND
diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"

# R scripts and paths
export R_LIBS="$SAMSA/R_scripts/packages"
R_programs=$SAMSA/R_scripts

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#       Example: control_1.R1.fastq
#                control_1.R2.fastq
#
# Note: if using single-end sequencing, skip this step (comment out).

#for file in $starting_files_location/*.gz
#do
#   gunzip $file
#done

cd $starting_files_location
for file in $starting_files_location/*R1*
do
    file1=$file
    file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'`
    out_path=`echo $file | awk -F "_R1" '{print $1 ".merged"}'`
    out_name=`echo ${out_path##*/}`
    checked $PEAR -f $file1 -r $file2 -o $out_name
done

[[ ! -d "$OUT_DIR/step_1_output_test" ]] && mkdir $OUT_DIR/step_1_output_test/
mv $starting_files_location/*merged* $OUT_DIR/step_1_output_test/
echo -e "\nPaired-end merging step completed.\n"

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $OUT_DIR/step_1_output_test/ folder!

for file in $OUT_DIR/step_1_output_test/*merged*
do
    shortname=`echo $file | awk -F "merged" '{print $1 "cleaned.fastq"}'`
    checked java -jar $TRIMMOMATIC SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

[[ ! -d "$OUT_DIR/step_2_output_test" ]] && mkdir $OUT_DIR/step_2_output_test/
mv $OUT_DIR/step_1_output_test/*cleaned.fastq $OUT_DIR/step_2_output_test/
echo -e "\nCleaning files with Trimmomatic completed.\n"

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [ -f $OUT_DIR/step_2_output_test/raw_counts.txt ]
then
    rm $OUT_DIR/step_2_output_test/raw_counts.txt
    touch $OUT_DIR/step_2_output_test/raw_counts.txt
else
    touch $OUT_DIR/step_2_output_test/raw_counts.txt
fi

for file in $OUT_DIR/step_2_output_test/*cleaned.fastq
do
    checked python $PY_DIR/raw_read_counter.py -I $file -O $OUT_DIR/step_2_output_test/raw_counts.txt
done
echo -e "\nCounting raw sequences completed!\n"

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $OUT_DIR/step_2_output_test/*cleaned.fastq
do
    shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`
    checked $SORTMERNA --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db --reads $file --aligned $file.ribosomes --other $shortname --fastx --num_alignments 0 --log -v
done

[[ ! -d "$OUT_DIR/step_3_output_test" ]] && mkdir $OUT_DIR/step_3_output_test/
mv $OUT_DIR/step_2_output_test/*ribodepleted* $OUT_DIR/step_3_output_test/

echo -e "\nRibosomal read removal step completed!\n"

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $OUT_DIR/step_3_output_test/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
    echo "Now starting on " $file
    echo "Converting to " $shortname
    checked $DIAMOND blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1
    checked $DIAMOND view --daa $file.RefSeq.daa -o $shortname -f tab
done

[[ ! -d "$OUT_DIR/step_4_output_test" ]] && mkdir $OUT_DIR/step_4_output_test/
[[ ! -d "$OUT_DIR/step_4_output_test/daa_binary_files" ]] &&
mkdir $OUT_DIR/step_4_output_test/daa_binary_files/

mv $OUT_DIR/step_3_output_test/*annotated* $OUT_DIR/step_4_output_test/
mv $OUT_DIR/step_3_output_test/*.daa $OUT_DIR/step_4_output_test/daa_binary_files/

echo -e "\nRefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $OUT_DIR/step_4_output_test/*RefSeq_annotated
do
    python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $diamond_database.fa -O
    python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $diamond_database.fa -F
done

[[ ! -d "$OUT_DIR/step_5_output_test/RefSeq_results/org_results" ]] && mkdir -p $OUT_DIR/step_5_output_test/RefSeq_results/org_results/
[[ ! -d "$OUT_DIR/step_5_output_test/RefSeq_results/func_results/" ]] && mkdir $OUT_DIR/step_5_output_test/RefSeq_results/func_results/
mv $OUT_DIR/step_4_output_test/*organism.tsv $OUT_DIR/step_5_output_test/RefSeq_results/org_results/
mv $OUT_DIR/step_4_output_test/*function.tsv $OUT_DIR/step_5_output_test/RefSeq_results/func_results/

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $OUT_DIR/step_3_output_test/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
    echo "Now starting on Subsystems annotations for " $file
    checked $DIAMOND blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1
    checked $DIAMOND view --daa $file.Subsys.daa -o $shortname -f tab
done

mv $OUT_DIR/step_3_output_test/*subsys_annotated* $OUT_DIR/step_4_output_test/
mv $OUT_DIR/step_3_output_test/*.daa $OUT_DIR/step_4_output_test/daa_binary_files/

echo -e "\nDIAMOND Subsystems annotations completed at: "; date

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $OUT_DIR/step_4_output_test/*subsys_annotated*
do
    checked python $PY_DIR/DIAMOND_subsystems_analysis_counter.py -I $file -D $diamond_subsys_db.fa -O $file.hierarchy -P $file.receipt

    # This quick program reduces down identical hierarchy annotations
    checked python $PY_DIR/subsys_reducer.py -I $file.hierarchy
done

[[ ! -d "$OUT_DIR/step_5_output_test/Subsystems_results/receipts/" ]] && mkdir -p $OUT_DIR/step_5_output_test/Subsystems_results/receipts/
mv $OUT_DIR/step_4_output_test/*.reduced $OUT_DIR/step_5_output_test/Subsystems_results/
mv $OUT_DIR/step_4_output_test/*.receipt $OUT_DIR/step_5_output_test/Subsystems_results/receipts/
rm $OUT_DIR/step_4_output_test/*.hierarchy

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

[[ ! -d "$OUT_DIR/step_5_output_test/RefSeq_results/org_results" ]] && mkdir -p "$OUT_DIR/step_5_output_test/RefSeq_results/org_results"
[[ ! -d "$OUT_DIR/step_5_output_test/RefSeq_results/func_results" ]] && mkdir "$OUT_DIR/step_5_output_test/RefSeq_results/func_results"
checked Rscript $R_programs/run_DESeq_stats.R -I $OUT_DIR/step_5_output_test/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab -R $OUT_DIR/step_2_output_test/raw_counts.txt
checked Rscript $R_programs/run_DESeq_stats.R -I $OUT_DIR/step_5_output_test/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab -R $OUT_DIR/step_2_output_test/raw_counts.txt
checked Rscript $R_programs/Subsystems_DESeq_stats.R -I $OUT_DIR/step_5_output_test/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1 -R $OUT_DIR/step_2_output_test/raw_counts.txt

if [[ -z "$TEST_NO_RM" ]]; then
  echo "Now removing output files."
  rm -r $OUT_DIR
fi
echo "SUCCESS! test_of_master_script.bash successfully ran."
exit
##################################################################
