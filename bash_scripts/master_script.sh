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
#   0. Set starting location and starting files pathway
source "$(dirname "$0")/common.sh"

INPUT_FILES=$SAMSA/input_files

#   4. DIAMOND
diamond_database="$SAMSA/full_databases/RefSeq_bac"
diamond_subsys_db="$SAMSA/full_databases/subsys_db"

#   5. Aggregation
RefSeq_db="$SAMSA/full_databases/RefSeq_bac.fa"
Subsys_db="$SAMSA/full_databases/subsys_db.fa"

#   6. R scripts and paths
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
# Note: if performing R analysis (step 6), be sure to name files with
#   the appropriate prefix ("control_$file" and "experimental_$file")!

for file in $INPUT_FILES/*.gz
do
    gunzip $file
done

for file in $INPUT_FILES/*_R1*
do
    file1=$file
    file2=`echo $file1 | awk -F "R1" '{print $1 "R2" $2}'`
    out_path=`echo $file | awk -F "_R1" '{print $1 ".merged"}'`
    out_name=`echo ${out_path##*/}`

    $PEAR -f $file1 -r $file2 -o $out_name
done

mkdir $SAMSA/step_1_output/
mv $INPUT_FILES/*merged* $SAMSA/step_1_output/

####################################################################
#
# STEP 2: CLEANING FILES WITH TRIMMOMATIC
# Note: if skipping PEAR, make sure that all starting files are in the
# $SAMSA/step_1_output/ folder!

for file in $SAMSA/step_1_output/*.merged*
do
    shortname=`echo $file | awk -F "merged" '{print $1 "cleaned.fastq"}'`
    java -jar $TRIMMOMATIC SE -phred33 $file $shortname SLIDINGWINDOW:4:15 MINLEN:99
done

mkdir $SAMSA/step_2_output/
mv $SAMSA/step_1_output/*cleaned.fastq $SAMSA/step_2_output/

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [ -f $SAMSA/step_2_output/raw_counts.txt ]
then
    rm $SAMSA/step_2_output/raw_counts.txt
    touch $SAMSA/step_2_output/raw_counts.txt
else
    touch $SAMSA/step_2_output/raw_counts.txt
fi

for file in $SAMSA/step_2_output/*.cleaned.fastq
do
    python $PY_DIR/raw_read_counter.py -I $file -O $SAMSA/step_2_output/raw_counts.txt
done

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $SAMSA/step_2_output/*.cleaned.fastq
do
    shortname=`echo $file | awk -F "cleaned" '{print $1 "ribodepleted"}'`
    $SORTMERNA --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db --reads $file --aligned $file.ribosomes --other $shortname --fastx --num_alignments 0 --log -v

done

mkdir $SAMSA/step_3_output/
mv $SAMSA/step_2_output/*ribodepleted* $SAMSA/step_3_output/

####################################################################
#
# STEP 4: ANNOTATING WITH DIAMOND AGAINST REFSEQ
# Note: this step assumes that the DIAMOND database is already built.  If not,
# do that first before running this step.

echo "Now starting on DIAMOND org annotations at: "; date

for file in $SAMSA/step_3_output/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "RefSeq_annotated"}'`
    echo "Now starting on " $file
    echo "Converting to " $shortname
    $DIAMOND blastx --db $diamond_database -q $file -a $file.RefSeq -t ./ -k 1
    $DIAMOND view --daa $file.RefSeq.daa -o $shortname -f tab
done

mkdir -p $SAMSA/step_4_output/daa_binary_files/

mv $SAMSA/step_3_output/*annotated* $SAMSA/step_4_output/
mv $SAMSA/step_3_output/*.daa $SAMSA/step_4_output/daa_binary_files/

echo "RefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $SAMSA/step_4_output/*RefSeq_annotated
do
    python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -O
    python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
done

mkdir -p $SAMSA/step_5_output/RefSeq_results/org_results/
mkdir $SAMSA/step_5_output/RefSeq_results/func_results/
mv $SAMSA/step_4_output/*organism.tsv $SAMSA/step_5_output/RefSeq_results/org_results/
mv $SAMSA/step_4_output/*function.tsv $SAMSA/step_5_output/RefSeq_results/func_results/

####################################################################
#
# STEP 4.1: ANNOTATING WITH DIAMOND AGAINST SUBSYSTEMS

echo "Now starting on DIAMOND Subsystems annotations at: "; date

for file in $SAMSA/step_3_output/*ribodepleted.fastq
do
    shortname=`echo $file | awk -F "ribodepleted" '{print $1 "subsys_annotated"}'`
    echo "Now starting on Subsystems annotations for " $file
    $DIAMOND blastx --db $diamond_subsys_db -q $file -a $file.Subsys -t ./ -k 1
    $DIAMOND view --daa $file.Subsys.daa -o $shortname -f tab
done

mv $SAMSA/step_3_output/*subsys_annotated* $SAMSA/step_4_output/
mv $SAMSA/step_3_output/*.daa $SAMSA/step_4_output/daa_binary_files/

echo "DIAMOND Subsystems annotations completed at: "; date

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $SAMSA/step_4_output/*subsys_annotated*
do
    python $PY_DIR/DIAMOND_subsystems_analysis_counter.py -I $file -D $Subsys_db -O $file.hierarchy -P $file.receipt

    # This quick program reduces down identical hierarchy annotations
    python $PY_DIR/subsys_reducer.py -I $file.hierarchy
done

mkdir -p $SAMSA/step_5_output/Subsystems_results/receipts/
mv $SAMSA/step_4_output/*.reduced $SAMSA/step_5_output/Subsystems_results/
mv $SAMSA/step_4_output/*.receipt $SAMSA/step_5_output/Subsystems_results/receipts/
rm $SAMSA/step_4_output/*.hierarchy

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

Rscript $R_programs/run_DESeq_stats.R -I $SAMSA/step_5_output/RefSeq_results/org_results/ -O RefSeq_org_DESeq_results.tab -R $SAMSA/step_2_output/raw_counts.txt
Rscript $R_programs/run_DESeq_stats.R -I $SAMSA/step_5_output/RefSeq_results/func_results/ -O RefSeq_func_DESeq_results.tab -R $SAMSA/step_2_output/raw_counts.txt
Rscript $R_programs/Subsystems_DESeq_stats.R -I $SAMSA/step_5_output/Subsystems_results/ -O Subsystems_level-1_DESeq_results.tab -L 1 -R $SAMSA/step_2_output/raw_counts.txt

echo "Master bash script finished running at: "; date
exit
####################################################################
