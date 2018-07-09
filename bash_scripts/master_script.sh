#!/bin/bash
#SBATCH --mem=100000
#SBATCH --time=7-0:0:0

# Lines starting with #SBATCH are for SLURM job management systems
# and may be removed if user is not submitting jobs to SLURM

####################################################################
#
# master_script.sh
# Created April 2017 by Sam Westreich, github.com/transcript
# This version modified February 21, 2018
#
####################################################################
#
# This script sets up and runs through ALL steps in the SAMSA pipeline
# before the analysis (which is done in R, likely in RStudio).  Each
# step is set up below.
#
# The steps are:
#   1.   Merging with PEAR, if applicable
#   2.   Read cleaning with Trimmomatic
#   3.   rRNA removal with SortMeRNA
#   4.   Annotation using DIAMOND (by default against the RefSeq database)
#   5.   Aggregation using analysis_counter.py
#   4.1  Annotation using DIAMOND against the Subsystems database
#   5.1  Aggregation using Subsystems-specific analysis counter.py
#   6.   Running R scripts to get DESeq statistical analysis.
#
# NOTE: BEFORE running this script, please run package_installation.bash
# and full_database_download.bash located at:
# https://github.com/transcript/samsa2/tree/master/setup in order to set
# up SAMSA2 dependencies and download full databases.
#
####################################################################
#
echo -e "NOTE: Before running this script, please run package_installation.bash and full_database_download.bash located at https://github.com/transcript/samsa2/tree/master/setup in order to set up SAMSA2 dependencies.\n\n"
#
# VARIABLES - set starting location and starting files location pathways
#
source "${BASH_SOURCE%/*}/../bash_scripts/common.sh"

INPUT_DIR=$SAMSA/input_files
OUT_DIR=$SAMSA

STEP_1="$OUT_DIR/step_1_output_test"
STEP_2="$OUT_DIR/step_2_output_test"
STEP_3="$OUT_DIR/step_3_output_test"
STEP_4="$OUT_DIR/step_4_output_test"
STEP_5="$OUT_DIR/step_5_output_test"

# Diamond databases
diamond_database="$SAMSA/full_databases/RefSeq_bac"
diamond_subsys_db="$SAMSA/full_databases/subsys_db"

# Aggregation databases
RefSeq_db="$SAMSA/full_databases/RefSeq_bac.fa"
Subsys_db="$SAMSA/full_databases/subsys_db.fa"

####################################################################
#
# STEP 1: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#       Example: control_1.R1.fastq
#                control_1.R2.fastq
#
# Note: if using single-end sequencing, skip this step (comment out).
# Note: if performing R analysis (step 6), be sure to name files with
#   the appropriate prefix ("control_$file" and "experimental_$file")!

for file in $INPUT_DIR/*.gz
do
    gunzip $file
done

for f in $INPUT_DIR/*_R1*
do
    f2=`echo $f | awk -F "R1" '{print $1 "R2" $2}'`
    out_path=`echo $f | awk -F "_R1" '{print $1 ".merged"}'`

    checked $PEAR -f $f -r $f2 -o ${out_path##*/}
done

$MKDIR $STEP_1
mv $INPUT_DIR/*merged* $STEP_1

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $STEP_1/ folder!

for file in $STEP_1/*.merged*
do
    shortname=`echo $file | awk -F "merged" '{print $1 "cleaned.fastq"}'`
    checked java -jar $TRIMMOMATIC SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

$MKDIR $STEP_2
mv $STEP_1/*cleaned.fastq $STEP_2

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [[ -f $STEP_2/raw_counts.txt ]]; then
    rm $STEP_2/raw_counts.txt
fi
touch $STEP_2/raw_counts.txt

for file in $STEP_2/*.cleaned.fastq
do
    checked python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
done

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $STEP_2/*.cleaned.fastq
do
    shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`
    checked $SORTMERNA \
      --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db \
      --reads $file --aligned $file.ribosomes --other $shortname --fastx \
      --num_alignments 0 --log -v
done

$MKDIR $STEP_3
mv $STEP_2/*ribodepleted* $STEP_3

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $STEP_3/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
    echo "Now starting on " $file
    echo "Converting to " $shortname
    checked $DIAMOND blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1
    checked $DIAMOND view --daa $file.RefSeq.daa -o $shortname -f tab
done

$MKDIR $STEP_4/daa_binary_files

mv $STEP_3/*annotated* $STEP_4
mv $STEP_3/*.daa $STEP_4/daa_binary_files

echo "RefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $STEP_4/*RefSeq_annotated
do
    checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
    checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
done

$MKDIR $STEP_5/RefSeq_results/org_results
$MKDIR $STEP_5/RefSeq_results/func_results
mv $STEP_4/*organism.tsv $STEP_5/RefSeq_results/org_results
mv $STEP_4/*function.tsv $STEP_5/RefSeq_results/func_results

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $STEP_3/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
    echo "Now starting on Subsystems annotations for " $file
    checked $DIAMOND blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1
    checked $DIAMOND view --daa $file.Subsys.daa -o $shortname -f tab
done

mv $STEP_3/*subsys_annotated* $STEP_4
mv $STEP_3/*.daa $STEP_4/daa_binary_files

echo "DIAMOND Subsystems annotations completed at: "; date

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $STEP_4/*subsys_annotated*
do
    checked python $PY_DIR/DIAMOND_subsystems_analysis_counter.py -I $file \
      -D $Subsys_db -O $file.hierarchy -P $file.receipt

    # This quick program reduces down identical hierarchy annotations
    checked python $PY_DIR/subsys_reducer.py -I $file.hierarchy
done

$MKDIR $STEP_5/Subsystems_results/receipts
mv $STEP_4/*.reduced $STEP_5/Subsystems_results
mv $STEP_4/*.receipt $STEP_5/Subsystems_results/receipts
rm $STEP_4/*.hierarchy

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

checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/org_results \
  -O RefSeq_org_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/run_DESeq_stats.R \
  -I $STEP_5/RefSeq_results/func_results \
  -O RefSeq_func_DESeq_results.tab \
  -R $STEP_2/raw_counts.txt
checked Rscript $R_DIR/Subsystems_DESeq_stats.R \
  -I $STEP_5/Subsystems_results \
  -O Subsystems_level-1_DESeq_results.tab -L 1 \
  -R $STEP_2/raw_counts.txt

echo "Master bash script finished running at: "; date
exit
####################################################################
