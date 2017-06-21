#!/bin/bash
#SBATCH --mem=32000
#SBATCH --time=7-0:0:0
#SBATCH --mail-user=stwestreich@ucdavis.edu
#SBATCH --mail-type=END

####################################################################
#
# Functions_by_organism.sh
# Created June 17, 2017 by Sam Westreich, github.com/transcript
# This version modified June 17, 2017
#
####################################################################
#
# Functional search by organism - getting DESeq functional results for
# an organism (or subset of organisms)
#
####################################################################
#
# VARIABLES

python_programs=/share/milklab/sam/python_scripts

specific_organism="Prevotella"

RefSeq_db="/share/milklab/sam/databases/bct"

RefSeq_results_location=/share/milklab/sam/test_files/output_diamond

raw_counts_file=/share/milklab/sam/test_files/raw_counts.txt
# note: generated at step 2.9 in master_script.sh

export R_LIBS="/share/milklab/sam/R_scripts/packages"
R_programs=/share/milklab/sam/R_scripts

#
####################################################################
#
# STEP 1 - run Python to subset the RefSeq results for organism of interest

for file in $RefSeq_results_location/*.RefSeq
do
	python $python_programs/DIAMOND_specific_organism_retriever.py -I $file -SO $specific_organism -D $RefSeq_db 
done

#
####################################################################
#
# STEP 2 - run Python to aggregate these results
for file in $RefSeq_results_location/*$specific_organism.tsv
do
	python $python_programs/DIAMOND_analysis_counter.py -I $file -D $RefSeq_db -F
done

#
####################################################################
#
# STEP 2.9 - move results to their own directories

mkdir $RefSeq_results_location/$specific_organism
mkdir $RefSeq_results_location/$specific_organism/aggregated

mv $RefSeq_results_location/*$specific_organism.tsv $RefSeq_results_location/$specific_organism/
mv $RefSeq_results_location/*$specific_organism*function.tsv $RefSeq_results_location/$specific_organism/aggregated

#
####################################################################
#
# STEP 3 - run R to get DESeq results for the new files
#
# NOTE: Make sure that files have the appropriate prefixes ('control_' 
#	or 'experimental_') so that R can recognize them!

Rscript $R_programs/run_DESeq_stats.R -I $RefSeq_results_location/$specific_organism/aggregated/ -O RefSeq.$specific_organism.func_DESeq_results.tab -R $raw_counts_file


