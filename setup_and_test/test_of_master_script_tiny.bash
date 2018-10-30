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
#   1.   Read cleaning with Trimmomatic
#   2.   Merging with PEAR, if applicable
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
source "${BASH_SOURCE%/*}/../bash_scripts/lib/common.sh"

# Starting files location
INPUT_DIR=$SAMSA/sample_files_paired-end/1_starting_files
OUT_DIR=$SAMSA/output_test
$MKDIR "$OUT_DIR"

STEP_1="$OUT_DIR/step_1_output_test"
STEP_2="$OUT_DIR/step_2_output_test"
STEP_3="$OUT_DIR/step_3_output_test"
STEP_4="$OUT_DIR/step_4_output_test"
STEP_5="$OUT_DIR/step_5_output_test"

# DIAMOND
diamond_database="$SAMSA/setup_and_test/tiny_databases/RefSeq_bac_TINY_24MB"
diamond_subsys_db="$SAMSA/setup_and_test/tiny_databases/subsys_db_TINY_24MB"

####################################################################
#
# STEP 1: CLEANING FILES WITH TRIMMOMATIC

#for file in $INPUT_DIR/*.gz
#do
#   gunzip $file
#done

$MKDIR $STEP_1
paired=false
for f in $INPUT_DIR/*R1*q
do
    f2=`echo $f | awk -F "R1" '{print $1 "R2" $2}'`
    out_path=`echo $f | awk -F "_R1" '{print $1 ".cleaned"}'`
    if [ -f $f2 ]; then
      paired=true
      checked java -jar $TRIMMOMATIC PE -phred33 $f $f2 \
        $out_path".forward" $out_path".forward_unpaired" $out_path".reverse" $out_path".reverse_unpaired" \
        SLIDINGWINDOW:4:15 MINLEN:70
    else
      checked java -jar $TRIMMOMATIC SE -phred33 $f $out_path SLIDINGWINDOW:4:15 MINLEN:70
    fi
done

if $paired; then
  mv $INPUT_DIR/*".cleaned.forward"* $STEP_1
  mv $INPUT_DIR/*".cleaned.reverse"* $STEP_1
else
  mv $INPUT_DIR/*".cleaned" $STEP_1
fi
echo -e "\nCleaning files with Trimmomatic completed.\n"

####################################################################
#
# STEP 2: MERGING OF PAIRED-END FILES USING PEAR
# Note: paired-end files are usually named using R1 and R2 in the name.
#       Example: control_1.R1.fastq
#                control_1.R2.fastq

$MKDIR $STEP_2
if $paired; then
  for file in $STEP_1/*.cleaned.forward
  do
      f2=`echo $file | awk -F "cleaned.forward" '{print $1 "cleaned.reverse"}'`
      shortname=`echo $file | awk -F "cleaned.forward" '{print $1 "merged"}'`
      checked $PEAR -f $file -r $f2 -o $STEP_2/${shortname##*/}
  done
  echo -e "\nPaired-end merging step completed.\n"
else
  for file in $STEP_1/*.cleaned
  do
    new_name=`echo $file | awk -F "cleaned" '{print $1 "merged.assembled.fastq"}'`
    cp $file $STEP_2/${new_name##*/}
  done
fi

####################################################################
#
# STEP 2.9: GETTING RAW SEQUENCES COUNTS
# Note: These are used later for statistical analysis.

if [[ -f $STEP_2/raw_counts.txt ]]; then
    rm $STEP_2/raw_counts.txt
fi
touch $STEP_2/raw_counts.txt

if $paired; then
  for file in $STEP_1/*cleaned.forward
  do
      checked python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
else
  for file in $STEP_1/*cleaned
  do
      checked python $PY_DIR/raw_read_counter.py -I $file -O $STEP_2/raw_counts.txt
  done
fi
echo -e "\nCounting raw sequences completed!\n"

####################################################################
#
# STEP 3: REMOVING RIBOSOMAL READS WITH SORTMERNA
# Note: this step assumes that the SortMeRNA databases are indexed.  If not,
# do that first (see the SortMeRNA user manual for details).

for file in $STEP_2/*assembled.fastq
do
    shortname=`echo $file | awk -F "merged" '{print $1 "ribodepleted"}'`
    checked $SORTMERNA \
      --ref $SORTMERNA_DIR/rRNA_databases/silva-bac-16s-id90.fasta,$SORTMERNA_DIR/index/silva-bac-16s-db \
      --reads $file --aligned $file.ribosomes --other $shortname --fastx \
      --num_alignments 0 --log -v
done

$MKDIR $STEP_3
mv $STEP_2/*ribodepleted* $STEP_3

echo -e "\nRibosomal read removal step completed!\n"

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

$MKDIR $STEP_4
$MKDIR $STEP_4/daa_binary_files

mv $STEP_3/*annotated* $STEP_4
mv $STEP_3/*.daa $STEP_4/daa_binary_files

echo -e "\nRefSeq DIAMOND annotations completed at: "; date

####################################################################
#
# STEP 5: AGGREGATING WITH ANALYSIS_COUNTER

for file in $STEP_4/*RefSeq_annotated
do
    checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $diamond_database.fa -O
    checked python $PY_DIR/DIAMOND_analysis_counter.py -I $file -D $diamond_database.fa -F
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

echo -e "\nDIAMOND Subsystems annotations completed at: "; date

##################################################################
#
# STEP 5.1: PYTHON SUBSYSTEMS ANALYSIS COUNTER

for file in $STEP_4/*subsys_annotated
do
    checked python $PY_DIR/DIAMOND_subsystems_analysis_counter.py -I $file \
      -D $diamond_subsys_db.fa -O $file.hierarchy -P $file.receipt

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

# XXX: Sample data does not have control/experimental files; skip this step
#checked Rscript $R_DIR/run_DESeq_stats.R \
#  -I $STEP_5/RefSeq_results/org_results \
#  -O RefSeq_org_DESeq_results.tab \
#  -R $STEP_2/raw_counts.txt
#checked Rscript $R_DIR/run_DESeq_stats.R \
#  -I $STEP_5/RefSeq_results/func_results \
#  -O RefSeq_func_DESeq_results.tab \
#  -R $STEP_2/raw_counts.txt
#checked Rscript $R_DIR/Subsystems_DESeq_stats.R \
#  -I $STEP_5/Subsystems_results \
#  -O Subsystems_level-1_DESeq_results.tab -L 1 \
#  -R $STEP_2/raw_counts.txt

if [[ -z "$TEST_NO_RM" ]]; then
  echo "Now removing output files."
  rm -r $OUT_DIR
fi
echo "SUCCESS! test_of_master_script.bash successfully ran."
exit
##################################################################
